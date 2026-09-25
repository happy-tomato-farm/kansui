"""
5分ごとの計算を1日分積算して、日単位・ハウス単位の値を出す。

【この段階で初めて出る数字】
    ・1日の総蒸散量 [L/m²/日] と [L/ハウス/日]
    ・1日の総光合成量 [g-CO2/m²/日] と [g-糖/m²/日]、およびハウス単位

【積算のしかた】
センサーは5分間隔なので、各時点の速度が次の5分間続くとみなして足す。

    1日の蒸散量 = Σ（各時点の蒸散速度 × 300秒）

これは区分求積法の一種で、5分という細かさなら誤差は無視できる。

【前提条件】
1. 日の境界は 0時。積算日射のリセットに合わせている。
2. 欠測がある日は、その日の積算値が過小になる。
   データが揃っている日だけを対象にするため、
   1日288点そろっているかを確認して記録する。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from config import (
    EXTINCTION_COEFFICIENT_K,
    HOUSE_SPECS,
    LEAF_CHARACTERISTIC_LENGTH_M,
    LOCAL_PRESSURE_KPA,
    WIND_SPEED_M_PER_S,
)
from core.canopy import (
    DEFAULT_LAYERS,
    calculate_canopy_gas_exchange,
    photosynthesis_to_co2_g_per_m2_s,
    photosynthesis_to_sugar_g_per_m2_s,
)

# 1日に期待される測定点の数（5分間隔なら 24時間 × 12回 = 288点）
EXPECTED_POINTS_PER_DAY = 288

# 積算の刻み [秒]
TIME_STEP_S = 300


@dataclass(frozen=True)
class SimulationSettings:
    """
    計算の設定をひとまとめにしたもの。

    引数が多くなりすぎるのを防ぎ、どの条件で計算した結果かを
    あとから確認できるようにする。
    """

    lai: float
    house: str = "中央"
    wind_speed_m_per_s: float = WIND_SPEED_M_PER_S
    extinction_coefficient: float = EXTINCTION_COEFFICIENT_K
    n_layers: int = DEFAULT_LAYERS
    pressure_kpa: float = LOCAL_PRESSURE_KPA
    characteristic_length_m: float = LEAF_CHARACTERISTIC_LENGTH_M

    # 葉温をエネルギー収支から求めるか。
    # False なら葉温＝気温（光合成.xlsx と同じ扱い）だが、
    # 実データでは蒸散量が物理的に小さくなりすぎるため、
    # 実運用では True を使う。
    solve_energy_balance: bool = True

    def floor_area_m2(self) -> float:
        """対象ハウスの床面積 [m²]。"""
        if self.house not in HOUSE_SPECS:
            raise ValueError(
                f"ハウス名が config.HOUSE_SPECS にありません: '{self.house}'。"
                f"使えるのは {list(HOUSE_SPECS.keys())} です"
            )
        return HOUSE_SPECS[self.house]["floor_area_m2"]

    def describe(self) -> str:
        """設定内容を日本語で説明する。"""
        leaf_temp_mode = (
            "エネルギー収支で解く" if self.solve_energy_balance else "葉温＝気温"
        )
        return (
            f"ハウス={self.house}（床面積 {self.floor_area_m2():.0f} m²） / "
            f"LAI={self.lai:.2f} / 風速={self.wind_speed_m_per_s} m/s / "
            f"消散係数 k={self.extinction_coefficient} / {self.n_layers}層 / "
            f"葉温は{leaf_temp_mode}"
        )


def simulate_timeseries(
    df: pd.DataFrame,
    settings: SimulationSettings,
) -> pd.DataFrame:
    """
    5分ごとの環境データから、各時点の群落ガス交換を計算する。

    Args:
        df: loader.add_derived_values() を通した DataFrame。
            timestamp, dry_bulb_c, vapor_pressure_kpa, solar_w_per_m2,
            co2_ppm の列が必要。
        settings: 計算の設定

    Returns:
        入力に計算結果の列を加えた新しい DataFrame。
        加わる列:
            transpiration_mg_per_m2_s   … 蒸散速度（地面あたり）
            net_photosynthesis_mmol     … 純光合成速度（地面あたり）
            gross_photosynthesis_mmol   … 総光合成速度
            respiration_mmol            … 呼吸（光呼吸＋暗呼吸）
    """
    required = [
        "timestamp", "dry_bulb_c", "vapor_pressure_kpa",
        "solar_w_per_m2", "co2_ppm",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"計算に必要な列がありません: {missing}。"
            f"loader.add_derived_values() を先に通してください"
        )

    result = df.copy()

    # 1点ずつ計算する。群落を層に分けて葉モデルを呼ぶため、
    # numpy でまとめて処理することができない。
    # 1日分（288点）で約0.5秒、1か月で約15秒が目安。
    transpiration = np.empty(len(result))
    net_photo = np.empty(len(result))
    gross_photo = np.empty(len(result))
    respiration = np.empty(len(result))
    mean_leaf_temp = np.empty(len(result))
    top_leaf_temp = np.empty(len(result))

    temps = result["dry_bulb_c"].to_numpy(dtype=float)
    vapor = result["vapor_pressure_kpa"].to_numpy(dtype=float)
    solar = result["solar_w_per_m2"].to_numpy(dtype=float)
    co2 = result["co2_ppm"].to_numpy(dtype=float)

    for i in range(len(result)):
        # 欠測があれば、その時点は計算せず NaN にする。
        # 0 で埋めると積算値が過小になるうえ、欠測だったことが
        # 分からなくなるため。
        if np.isnan(temps[i]) or np.isnan(vapor[i]) or np.isnan(co2[i]):
            transpiration[i] = np.nan
            net_photo[i] = np.nan
            gross_photo[i] = np.nan
            respiration[i] = np.nan
            mean_leaf_temp[i] = np.nan
            top_leaf_temp[i] = np.nan
            continue

        canopy = calculate_canopy_gas_exchange(
            air_temp_c=temps[i],
            vapor_pressure_kpa=vapor[i],
            solar_above_w_per_m2=solar[i],
            air_co2_ppm=co2[i],
            wind_speed_m_per_s=settings.wind_speed_m_per_s,
            pressure_kpa=settings.pressure_kpa,
            characteristic_length_m=settings.characteristic_length_m,
            lai=settings.lai,
            extinction_coefficient=settings.extinction_coefficient,
            n_layers=settings.n_layers,
            solve_energy_balance=settings.solve_energy_balance,
        )

        transpiration[i] = canopy.transpiration_mg_per_m2_s
        net_photo[i] = canopy.net_photosynthesis_mmol_per_m2_s
        gross_photo[i] = canopy.gross_photosynthesis_mmol_per_m2_s
        respiration[i] = canopy.respiration_mmol_per_m2_s
        mean_leaf_temp[i] = canopy.mean_leaf_temp_c
        top_leaf_temp[i] = canopy.top_layer_leaf_temp_c

    result["transpiration_mg_per_m2_s"] = transpiration
    result["net_photosynthesis_mmol"] = net_photo
    result["gross_photosynthesis_mmol"] = gross_photo
    result["respiration_mmol"] = respiration
    result["mean_leaf_temp_c"] = mean_leaf_temp
    result["top_leaf_temp_c"] = top_leaf_temp
    result["leaf_air_diff_c"] = mean_leaf_temp - result["dry_bulb_c"]

    return result


def aggregate_daily(
    timeseries: pd.DataFrame,
    settings: SimulationSettings,
) -> pd.DataFrame:
    """
    5分ごとの計算結果を1日分積算する。

    Args:
        timeseries: simulate_timeseries() の戻り値
        settings: 計算の設定（ハウス面積の換算に使う）

    Returns:
        日別の集計。列は以下のとおり。
            date                     … 日付
            points                   … その日の測定点数（288なら完全）
            is_complete              … 288点そろっているか
            transpiration_l_per_m2   … 蒸散量 [L/m²/日]（= mm/日）
            transpiration_l_per_house… 蒸散量 [L/ハウス/日]
            co2_fixed_g_per_m2       … CO2固定量 [g-CO2/m²/日]
            sugar_g_per_m2           … 糖換算 [g-糖/m²/日]
            sugar_kg_per_house       … 糖換算 [kg-糖/ハウス/日]
            daily_radiation_mj       … その日の積算日射 [MJ/m²]
            mean_temp_c, mean_vpd_kpa, mean_co2_ppm … 日平均の環境
    """
    work = timeseries.copy()
    work["date"] = work["timestamp"].dt.date

    floor_area = settings.floor_area_m2()

    rows = []
    for date, group in work.groupby("date"):
        points = len(group)

        # --- 蒸散量 -------------------------------------------------------
        # mg/(m²·s) × 300秒 = mg/m² を積み上げる。
        # 1 mg の水は 1 mm³ = 1e-6 L。
        transpiration_mg = np.nansum(
            group["transpiration_mg_per_m2_s"].to_numpy(dtype=float)
        ) * TIME_STEP_S
        transpiration_l_per_m2 = transpiration_mg * 1.0e-6

        # --- 光合成量 -----------------------------------------------------
        # mmol/(m²·s) × 300秒 = mmol/m² を積み上げる。
        net_mmol = np.nansum(
            group["net_photosynthesis_mmol"].to_numpy(dtype=float)
        ) * TIME_STEP_S
        gross_mmol = np.nansum(
            group["gross_photosynthesis_mmol"].to_numpy(dtype=float)
        ) * TIME_STEP_S
        respiration_mmol = np.nansum(
            group["respiration_mmol"].to_numpy(dtype=float)
        ) * TIME_STEP_S

        # 速度あたりの換算関数は「1秒あたり」を前提にしているので、
        # 積算値（mmol/m²）に対しては同じ係数を直接掛ける。
        co2_fixed_g = net_mmol / 1000.0 * 44.0
        sugar_g = net_mmol / 1000.0 * (180.0 / 6.0)

        rows.append({
            "date": date,
            "points": points,
            "is_complete": points == EXPECTED_POINTS_PER_DAY,
            "transpiration_l_per_m2": transpiration_l_per_m2,
            "transpiration_l_per_house": transpiration_l_per_m2 * floor_area,
            "net_photosynthesis_mmol_per_m2": net_mmol,
            "gross_photosynthesis_mmol_per_m2": gross_mmol,
            "respiration_mmol_per_m2": respiration_mmol,
            "co2_fixed_g_per_m2": co2_fixed_g,
            "co2_fixed_kg_per_house": co2_fixed_g * floor_area / 1000.0,
            "sugar_g_per_m2": sugar_g,
            "sugar_kg_per_house": sugar_g * floor_area / 1000.0,
            "daily_radiation_mj": group["cumulative_radiation_mj"].max(),
            "mean_temp_c": group["dry_bulb_c"].mean(),
            "mean_vpd_kpa": group["vpd_kpa"].mean() if "vpd_kpa" in group else np.nan,
            "mean_co2_ppm": group["co2_ppm"].mean(),
            "daytime_mean_vpd_kpa": (
                group.loc[group["solar_w_per_m2"] > 20, "vpd_kpa"].mean()
                if "vpd_kpa" in group else np.nan
            ),
        })

    return pd.DataFrame(rows)


def simulate_and_aggregate(
    df: pd.DataFrame,
    settings: SimulationSettings,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    5分ごとの計算と日別集計をまとめて行う。

    通常はこの関数を呼べばよい。

    Returns:
        (5分ごとの結果, 日別集計)
    """
    timeseries = simulate_timeseries(df, settings)
    daily = aggregate_daily(timeseries, settings)

    return timeseries, daily


def summarize_monthly(daily: pd.DataFrame) -> pd.DataFrame:
    """
    日別集計を月別にまとめる。

    季節変化を見るときや、潅水実績と比べるときに使う。
    """
    work = daily.copy()
    work["month"] = pd.to_datetime(work["date"]).dt.strftime("%Y-%m")

    # 欠測のある日は積算値が過小になるので、完全な日だけで平均する
    complete = work[work["is_complete"]]

    return complete.groupby("month").agg(
        日数=("date", "count"),
        蒸散_L_m2_日=("transpiration_l_per_m2", "mean"),
        蒸散_L_ハウス_日=("transpiration_l_per_house", "mean"),
        CO2固定_g_m2_日=("co2_fixed_g_per_m2", "mean"),
        糖換算_g_m2_日=("sugar_g_per_m2", "mean"),
        糖換算_kg_ハウス_日=("sugar_kg_per_house", "mean"),
        日射_MJ_m2=("daily_radiation_mj", "mean"),
        平均気温=("mean_temp_c", "mean"),
        日中VPD_kPa=("daytime_mean_vpd_kpa", "mean"),
    ).reset_index()
