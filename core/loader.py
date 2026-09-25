"""
センサーデータ（4週記録CSV）を読み込み、モデルが使える形に整える。

【元データについて】
栽培記録アプリ/saibai-kiroku/data/ にある 4週記録CSV。
5分間隔・1か月1ファイル・文字コードは cp932（Shift_JIS）。

    4週記録中央_YYYYMM.csv … 11列
    4週記録東_YYYYMM.csv   … 10列

【このファイルで扱う3つの落とし穴】

1. 湿度の列は使えない
   「湿度1」「絶対湿度1」「飽差1」の各列は壊れている。
   2026年4月以降はすべて 0 で、それ以前も日中は 0 になる。
   代わりに「温度1（乾球）」と「温度3（湿球）」から湿度を計算し直す。

2. 東ハウスには湿球温度の列がない
   東のファイルは、温度以外の列（湿度・CO2・照度・日射）がすべて
   中央のコピーになっている。センサーが1台しかないため。
   東の湿度は中央の湿球温度を借りて計算する。この近似の妥当性は
   check_data_quality() が件数として報告する。

3. 日射は「1日積算値」で、毎日0時にリセットされる
   瞬時の日射 [W/m²] を得るには差分を取る必要がある。
   0時をまたぐと差分が負になるので、そこは当日の積算値をそのまま使う。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from config import (
    LOCAL_PRESSURE_KPA,
    MAX_PLAUSIBLE_SOLAR_W_PER_M2,
    PSYCHROMETRIC_GAMMA_PER_K,
    TETENS_A_KPA,
    TETENS_B,
    TETENS_C_C,
    WET_BULB_ANOMALY_THRESHOLD,
)

# 測定間隔 [秒]。積算日射から瞬時日射を求めるときに使う。
MEASUREMENT_INTERVAL_S = 300

# ハウスごとの列名。中央は「温度1」、東は「温度2」を乾球として使う。
COLUMN_MAP = {
    "中央": {
        "dry_bulb": "温度1(区画1)_平均",
        "wet_bulb": "温度3(区画1)_平均",
    },
    "東": {
        "dry_bulb": "温度2(区画1)_平均",
        "wet_bulb": None,  # 東には湿球の列がない。中央から借りる。
    },
}

# ハウス共通の列名（東のファイルにも同じ名前で入っているが中身は中央のコピー）
COMMON_COLUMNS = {
    "co2_ppm": "CO2濃度(区画1)_平均",
    "illuminance": "照度(区画1)_平均",
    "cumulative_radiation_mj": "1日積算日射量(区画1)_平均",
}


@dataclass(frozen=True)
class DataQualityReport:
    """データ品質の点検結果。"""

    total_rows: int
    missing_dry_bulb: int          # 乾球温度の欠測
    missing_wet_bulb: int          # 湿球温度の欠測
    wet_bulb_above_dry: int        # 湿球 > 乾球（物理的にありえない）
    wet_bulb_clipped: int          # 湿球を乾球に合わせて飽和扱いにした件数
    radiation_resets: int          # 積算日射のリセット回数（≒日数）
    negative_radiation: int        # リセット以外で差分が負になった回数
    max_radiation_w_per_m2: float  # 瞬時日射の最大値（クリップ後）
    solar_clipped: int             # 日射が上限を超えてクリップされた件数
    reliable_days: int             # 湿球温度が信頼できる日数
    total_days: int                # 全日数

    def describe(self) -> str:
        """点検結果を日本語の文章にする。"""
        lines = [
            f"  総行数              : {self.total_rows:,} 行"
            f"（{self.total_days} 日分）",
            f"  乾球温度の欠測      : {self.missing_dry_bulb} 件",
            f"  湿球温度の欠測      : {self.missing_wet_bulb} 件",
            f"  湿球>乾球（飽和扱い）: {self.wet_bulb_clipped} 件"
            f"（{self.wet_bulb_clipped / max(self.total_rows, 1) * 100:.2f}%）",
            f"  積算日射のリセット  : {self.radiation_resets} 回",
            f"  日射の異常な減少    : {self.negative_radiation} 件",
            f"  日射の上限クリップ  : {self.solar_clipped} 件"
            f"（{self.solar_clipped / max(self.total_rows, 1) * 100:.3f}%）",
            f"  瞬時日射の最大値    : {self.max_radiation_w_per_m2:.0f} W/m²",
            f"  湿球が信頼できる日数: {self.reliable_days} / {self.total_days} 日"
            f"（{self.reliable_days / max(self.total_days, 1) * 100:.1f}%）",
        ]
        return "\n".join(lines)


def _saturation_vapor_pressure_array(temp_c: np.ndarray) -> np.ndarray:
    """
    飽和水蒸気圧をまとめて計算する（Tetens式）。

    psychrometry.saturation_vapor_pressure_kpa と同じ式だが、
    こちらは numpy 配列を一度に処理する。数万行を1行ずつ関数呼び出しすると
    非常に遅くなるため、前処理用にベクトル版を用意している。
    （両者が一致することは tests/test_loader.py で確認する）
    """
    return TETENS_A_KPA * np.exp(TETENS_B * temp_c / (temp_c + TETENS_C_C))


def load_sensor_csv(csv_path: Path | str, house: str) -> pd.DataFrame:
    """
    4週記録CSVを1ファイル読み込み、列名を英語に正規化する。

    Args:
        csv_path: CSVファイルのパス
        house: "中央" または "東"

    Returns:
        正規化した DataFrame。列は timestamp, dry_bulb_c, wet_bulb_c,
        co2_ppm, cumulative_radiation_mj。
        東の場合 wet_bulb_c は NaN（あとで中央から補完する）。
    """
    if house not in COLUMN_MAP:
        raise ValueError(
            f"ハウス名は '中央' か '東' を指定してください: 受け取った値 '{house}'"
        )

    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSVが見つかりません: {csv_path}")

    # 元データは Excel から出力された cp932（Shift_JIS）。
    # utf-8 で読もうとすると UnicodeDecodeError になる。
    try:
        raw = pd.read_csv(csv_path, encoding="cp932")
    except UnicodeDecodeError as e:
        raise ValueError(
            f"CSVの文字コードを cp932 として読めませんでした: {csv_path}。"
            f"元のファイルが別の文字コードで保存されている可能性があります"
        ) from e

    columns = COLUMN_MAP[house]

    # 必要な列がそろっているか確認する。
    # 列名が変わっていた場合、黙って NaN で埋めると原因が分からなくなるため、
    # ここで名前を挙げて止める。
    required = [columns["dry_bulb"]] + list(COMMON_COLUMNS.values())
    missing = [c for c in required if c not in raw.columns]
    if missing:
        raise ValueError(
            f"{csv_path.name} に必要な列がありません: {missing}\n"
            f"実際の列: {list(raw.columns)}"
        )

    result = pd.DataFrame({
        "timestamp": pd.to_datetime(raw["日付"]),
        "dry_bulb_c": pd.to_numeric(raw[columns["dry_bulb"]], errors="coerce"),
        "co2_ppm": pd.to_numeric(raw[COMMON_COLUMNS["co2_ppm"]], errors="coerce"),
        "cumulative_radiation_mj": pd.to_numeric(
            raw[COMMON_COLUMNS["cumulative_radiation_mj"]], errors="coerce"
        ),
    })

    if columns["wet_bulb"] is not None:
        result["wet_bulb_c"] = pd.to_numeric(
            raw[columns["wet_bulb"]], errors="coerce"
        )
    else:
        # 東には湿球の列がない。呼び出し側で中央から補完する。
        result["wet_bulb_c"] = np.nan

    result["house"] = house

    return result


def load_house_period(
    data_dir: Path | str,
    house: str,
    months: list[str],
) -> pd.DataFrame:
    """
    指定した月のCSVをまとめて読み込む。

    東ハウスの場合、湿球温度を中央のファイルから自動的に借りてくる。

    Args:
        data_dir: CSVが置いてあるフォルダ
        house: "中央" または "東"
        months: 読み込む月のリスト（例: ["202511", "202512"]）

    Returns:
        時刻順に並べた DataFrame
    """
    data_dir = Path(data_dir)
    frames = []

    for month in months:
        csv_path = data_dir / f"4週記録{house}_{month}.csv"
        frame = load_sensor_csv(csv_path, house)

        if house == "東":
            # 東には湿球がないので、同じ月の中央のファイルから借りる。
            central_path = data_dir / f"4週記録中央_{month}.csv"
            if not central_path.exists():
                raise FileNotFoundError(
                    f"東の湿度計算には中央の湿球温度が必要ですが、"
                    f"{central_path} がありません"
                )
            central = load_sensor_csv(central_path, "中央")

            # 時刻で突き合わせる。行数が同じでも欠測でずれる可能性があるため、
            # 位置ではなく時刻をキーにする。
            frame = frame.drop(columns=["wet_bulb_c"]).merge(
                central[["timestamp", "wet_bulb_c"]],
                on="timestamp",
                how="left",
            )

        frames.append(frame)

    if not frames:
        raise ValueError("読み込む月が1つも指定されていません")

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values("timestamp").reset_index(drop=True)

    return combined


def add_derived_values(df: pd.DataFrame) -> pd.DataFrame:
    """
    センサーの生の値から、モデルの入力に必要な量を計算して列に加える。

    元の DataFrame は変更せず、列を加えた新しい DataFrame を返す。

    加える列:
        wet_bulb_used_c      … 実際に計算に使った湿球温度（飽和クリップ後）
        was_clipped          … 湿球を乾球に合わせた（飽和扱いにした）かどうか
        saturation_kpa       … 乾球温度での飽和水蒸気圧 es(Ta)
        vapor_pressure_kpa   … 実際の水蒸気圧 ea
        vpd_kpa              … 飽差 VPD [kPa]
        vpd_g_per_m3         … 飽差 [g/m³]（現場の感覚に近い単位）
        relative_humidity    … 相対湿度（0〜1）
        solar_w_per_m2       … 瞬時日射 [W/m²]

    Returns:
        列を加えた新しい DataFrame
    """
    result = df.copy()

    dry = result["dry_bulb_c"].to_numpy(dtype=float)
    wet = result["wet_bulb_c"].to_numpy(dtype=float)

    # ---- 湿球温度の飽和クリップ --------------------------------------
    # 湿球が乾球を上回るのは物理的にありえないが、実データでは
    # 飽和付近（夜間・早朝）で 0.1〜0.4℃ 程度上回ることがある。
    # 東ハウスでは中央の湿球を借りているため、なおさら起きやすい。
    # これは「飽和している」ことを意味するので、湿球=乾球にそろえて
    # 相対湿度100%として扱う。
    was_clipped = wet > dry
    wet_used = np.where(was_clipped, dry, wet)

    result["wet_bulb_used_c"] = wet_used
    result["was_clipped"] = was_clipped

    # ---- 水蒸気圧と飽差 ------------------------------------------------
    es_dry = _saturation_vapor_pressure_array(dry)
    es_wet = _saturation_vapor_pressure_array(wet_used)

    # 乾湿計の式（式3.16）
    #   ea = es(Tw) - γ * Pa * (Ta - Tw)
    ea = es_wet - PSYCHROMETRIC_GAMMA_PER_K * LOCAL_PRESSURE_KPA * (dry - wet_used)

    # 計算誤差で ea が負や es を超えることがあるので、範囲に収める。
    ea = np.clip(ea, 0.0, es_dry)

    result["saturation_kpa"] = es_dry
    result["vapor_pressure_kpa"] = ea
    result["vpd_kpa"] = es_dry - ea

    # 相対湿度。es が 0 になることは温度域上ありえないが、念のため守る。
    with np.errstate(divide="ignore", invalid="ignore"):
        rh = np.where(es_dry > 0, ea / es_dry, np.nan)
    result["relative_humidity"] = np.clip(rh, 0.0, 1.0)

    # 飽差を g/m³ にも直す。1m³の空気があと何gの水を受け取れるか。
    # トマトでは 3〜7 g/m³ が目安とされる。
    from config import GAS_CONSTANT_J_PER_MOL_K, MOLAR_MASS_WATER_G_PER_MOL
    temp_k = dry + 273.15
    moles_per_m3 = (LOCAL_PRESSURE_KPA * 1000.0) / (GAS_CONSTANT_J_PER_MOL_K * temp_k)
    result["vpd_g_per_m3"] = (
        result["vpd_kpa"] / LOCAL_PRESSURE_KPA * moles_per_m3 * MOLAR_MASS_WATER_G_PER_MOL
    )

    # ---- 積算日射から瞬時日射へ ----------------------------------------
    solar = _cumulative_to_instant_radiation(
        result["cumulative_radiation_mj"].to_numpy(dtype=float),
        result["timestamp"],
    )

    # ハウス内でありえない大きさの日射はセンサーのノイズとみなす。
    # クリップした事実を列に残し、点検で件数を数えられるようにする。
    result["solar_was_clipped"] = solar > MAX_PLAUSIBLE_SOLAR_W_PER_M2
    result["solar_w_per_m2"] = np.minimum(solar, MAX_PLAUSIBLE_SOLAR_W_PER_M2)

    # ---- 湿球温度が信頼できる日かどうかの判定 --------------------------
    result["wet_bulb_reliable"] = _flag_reliable_days(result)

    return result


def _flag_reliable_days(df: pd.DataFrame) -> pd.Series:
    """
    日ごとに湿球温度が信頼できるかを判定する。

    【なぜ日単位で判定するのか】
    2025年7〜10月は湿球温度が壊れているが、期間で機械的に切ると
    9月・10月に混ざっている健全な日まで捨てることになる。
    定植直後の大事な時期なので、日ごとに見て使える日は残す。

    【判定のしかた】
    その日のうち「湿球 > 乾球」となった点の割合がしきい値を超えたら、
    その日の湿球温度は信用しない。飽和状態では計算誤差で数%は
    正常に現れるため、5%を境界にしている。

    Returns:
        各行に対する真偽値（その行の日が信頼できるか）
    """
    daily_anomaly_rate = (
        df.assign(date=df["timestamp"].dt.date)
        .groupby("date")["was_clipped"]
        .transform("mean")
    )

    return daily_anomaly_rate <= WET_BULB_ANOMALY_THRESHOLD


def summarize_daily_quality(df: pd.DataFrame) -> pd.DataFrame:
    """
    日ごとの品質を一覧にする。

    どの日が使えてどの日が使えないかを目で確認したいときに使う。

    Returns:
        date, anomaly_rate（湿球異常の割合）, reliable（使えるか）,
        daily_radiation_mj（その日の積算日射）を持つ DataFrame
    """
    work = df.assign(date=df["timestamp"].dt.date)

    summary = work.groupby("date").agg(
        anomaly_rate=("was_clipped", "mean"),
        reliable=("wet_bulb_reliable", "first"),
        daily_radiation_mj=("cumulative_radiation_mj", "max"),
        mean_temp_c=("dry_bulb_c", "mean"),
        mean_vpd_kpa=("vpd_kpa", "mean"),
    )

    return summary.reset_index()


def filter_reliable_days(df: pd.DataFrame) -> pd.DataFrame:
    """
    湿球温度が信頼できる日だけを残す。

    計算を回す前にこれを通すことで、壊れたデータが結果に
    混ざるのを防ぐ。
    """
    if "wet_bulb_reliable" not in df.columns:
        raise ValueError(
            "先に add_derived_values() を呼んでください。"
            "wet_bulb_reliable 列がありません"
        )

    return df[df["wet_bulb_reliable"]].reset_index(drop=True)


def _cumulative_to_instant_radiation(
    cumulative_mj: np.ndarray,
    timestamps: pd.Series,
) -> np.ndarray:
    """
    1日積算日射量 [MJ/m²] から瞬時日射 [W/m²] を求める。

    【理屈】
    積算値は0時に0へリセットされ、日中に増えていく。
    5分間の増分が、その5分間に受けたエネルギー。

        瞬時日射 [W/m²] = 増分 [MJ/m²] × 1,000,000 ÷ 300秒

    【0時をまたぐ扱い】
    リセットの瞬間だけ差分が大きく負になる。その行は
    「当日の積算値そのもの」を増分として扱う（通常は0なので日射も0）。
    """
    increment = np.diff(cumulative_mj, prepend=np.nan)

    # 日付が変わった行を特定する
    dates = timestamps.dt.date.to_numpy()
    is_new_day = np.empty(len(dates), dtype=bool)
    is_new_day[0] = True
    is_new_day[1:] = dates[1:] != dates[:-1]

    # 日付が変わった行は、前日からの差ではなく当日の積算値を使う
    increment = np.where(is_new_day, cumulative_mj, increment)

    # それでも負になる場合はセンサーの異常。0 として扱う。
    # （日射が負になることは物理的にありえない）
    increment = np.where(increment < 0, 0.0, increment)

    # MJ/m²（5分間） → W/m²（＝J/m²/s）
    return increment * 1.0e6 / MEASUREMENT_INTERVAL_S


def check_data_quality(df: pd.DataFrame) -> DataQualityReport:
    """
    読み込んだデータの品質を点検する。

    計算を回す前に、想定外のデータが混ざっていないかを確認するためのもの。
    ここで異常が多ければ、その期間は結果を信用しない判断ができる。
    """
    dry = df["dry_bulb_c"].to_numpy(dtype=float)
    wet = df["wet_bulb_c"].to_numpy(dtype=float)
    cumulative = df["cumulative_radiation_mj"].to_numpy(dtype=float)

    # 積算日射のリセット（差分が負）の回数
    increment = np.diff(cumulative, prepend=cumulative[0])
    dates = df["timestamp"].dt.date.to_numpy()
    is_new_day = np.empty(len(dates), dtype=bool)
    is_new_day[0] = True
    is_new_day[1:] = dates[1:] != dates[:-1]

    resets = int((increment < 0).sum())
    # 日付が変わったわけでもないのに減った件数（センサー異常の疑い）
    negative_not_reset = int(((increment < 0) & (~is_new_day)).sum())

    solar = _cumulative_to_instant_radiation(cumulative, df["timestamp"])
    solar_clipped = int((solar > MAX_PLAUSIBLE_SOLAR_W_PER_M2).sum())

    daily = df.assign(date=df["timestamp"].dt.date).groupby("date")
    total_days = daily.ngroups
    if "wet_bulb_reliable" in df.columns:
        reliable_days = int(daily["wet_bulb_reliable"].first().sum())
    else:
        reliable_days = 0

    return DataQualityReport(
        total_rows=len(df),
        missing_dry_bulb=int(np.isnan(dry).sum()),
        missing_wet_bulb=int(np.isnan(wet).sum()),
        wet_bulb_above_dry=int(np.nansum(wet > dry)),
        wet_bulb_clipped=int(df["was_clipped"].sum()) if "was_clipped" in df else 0,
        radiation_resets=resets,
        negative_radiation=negative_not_reset,
        max_radiation_w_per_m2=float(np.nanmin([
            np.nanmax(solar), MAX_PLAUSIBLE_SOLAR_W_PER_M2
        ])),
        solar_clipped=solar_clipped,
        reliable_days=reliable_days,
        total_days=total_days,
    )


def load_and_prepare(
    data_dir: Path | str,
    house: str,
    months: list[str],
) -> tuple[pd.DataFrame, DataQualityReport]:
    """
    読み込みから派生量の計算・品質点検までをまとめて行う。

    通常はこの関数を呼べばよい。

    Args:
        data_dir: CSVが置いてあるフォルダ
        house: "中央" または "東"
        months: 読み込む月（例: ["202511", "202512"]）

    Returns:
        (整形済みDataFrame, 品質レポート)
    """
    raw = load_house_period(data_dir, house, months)
    prepared = add_derived_values(raw)
    report = check_data_quality(prepared)

    return prepared, report
