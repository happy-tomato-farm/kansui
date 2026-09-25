"""
根群域の水収支を1日ずつ追いかける。

【何を出すためのものか】
潅水したうち、どれだけが
  ・蒸散として植物に使われたか
  ・根群域より下へ抜けて（地下へ流亡して）しまったか
  ・土に残ったか
を分けて示す。

さらに、土壌水分が蒸散を抑え始めているか（教科書 第9章 式9.19 の下の
  蒸散 = min(可能蒸散 E_pmax, 可能吸水 U_p)
という切り替え）を毎日判定する。

【収支の式】
根群域を1つの容器とみなす（タンクモデル）。

    貯留量の変化 ＝ 潅水 − 蒸散 − 排水

  貯留量 S [mm] ＝ 体積含水率 θ × 根群域の深さ d [mm]
  排水          ＝ その時の含水率での透水係数 K(θ)（単位勾配の仮定）

排水を「圃場容水量を超えた分は即座に落とす」とはしない。
K(θ) を使えば、どれくらいの速さで抜けるかが時間とともに追える。
実際この土では、飽和近くで 1日 90 mm 近く抜けるのに、
圃場容水量まで下がると 0.5 mm/日 まで落ちる。

【点滴チューブの扱い】
点滴は床全面を濡らさない。帯状に濡れた部分だけに水が入り、
根もそこに集まり、蒸散もそこから起きる。
そこで「濡れ面積率 f」を入れ、床面積あたりの量を濡れた部分あたりに
割り戻して計算する。

    濡れた部分での潅水 [mm] ＝ 床面積あたりの潅水 [mm] ÷ f

f が小さいほど、同じ潅水量でも濡れた部分は深く濡れるので、
排水（流亡）が増える。この効き方を確かめられるようにしてある。

【時間の刻み】
排水速度は含水率にきわめて敏感で、θ が 0.50 から 0.44 に下がるだけで
90 mm/日 → 1 mm/日 と100分の1になる。1日まとめて計算すると
潅水直後の速い排水が1日続いたことになってしまうので、5分刻みで進める。

潅水と蒸散は日中（既定 6〜17時）に配分し、夜間は排水だけが続く形にする。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from config import (
    DAYTIME_WINDOW_HOUR,
    ROOT_SYSTEM_MAX_UPTAKE_MM_PER_DAY,
    ROOT_ZONE_DEPTH_M,
    DRIP_WETTED_FRACTION,
    MIN_AIR_FILLED_POROSITY,
    SOIL_EVAPORATION_MM_PER_DAY,
    SOIL_RETENTION,
    WATER_BALANCE_SUBSTEPS_PER_DAY,
)
from core.soil import (
    FIELD_CAPACITY_WATER_CONTENT,
    air_filled_porosity,
    describe_soil_state,
    drainage_flux_mm_per_day,
    uptake_limit_ratio_from_water_content,
)

MM_PER_M = 1000.0


# =============================================================================
# 1. 設定と結果の入れ物
# =============================================================================

@dataclass(frozen=True)
class WaterBalanceSettings:
    """水収支の計算条件。

    根群域の深さと濡れ面積率は実測がないので、ここを振って
    感度解析をする（結論が前提に左右されないかを確かめる）。
    """

    root_zone_depth_m: float = ROOT_ZONE_DEPTH_M
    wetted_fraction: float = DRIP_WETTED_FRACTION
    max_uptake_mm_per_day: float = ROOT_SYSTEM_MAX_UPTAKE_MM_PER_DAY
    soil_evaporation_mm_per_day: float = SOIL_EVAPORATION_MM_PER_DAY
    substeps_per_day: int = WATER_BALANCE_SUBSTEPS_PER_DAY
    daytime_window_hour: tuple[float, float] = DAYTIME_WINDOW_HOUR

    def __post_init__(self) -> None:
        if not 0.0 < self.wetted_fraction <= 1.0:
            raise ValueError(
                f"濡れ面積率は 0 より大きく 1 以下でなければならない。"
                f"渡された値: {self.wetted_fraction}"
            )
        if self.root_zone_depth_m <= 0.0:
            raise ValueError(
                f"根群域の深さは正の値でなければならない。"
                f"渡された値: {self.root_zone_depth_m} m"
            )
        start, end = self.daytime_window_hour
        if not 0.0 <= start < end <= 24.0:
            raise ValueError(
                f"日中の時間帯は 0 ≦ 開始 < 終了 ≦ 24 でなければならない。"
                f"渡された値: {self.daytime_window_hour}"
            )

    @property
    def root_zone_depth_mm(self) -> float:
        return self.root_zone_depth_m * MM_PER_M


@dataclass(frozen=True)
class DailyWaterInput:
    """1日ぶんの入力。日付・潅水量・可能蒸散量。

    潅水量と可能蒸散量はどちらも **床面積あたり** [mm/日]（= L/m²/日）。
    """

    date: str
    irrigation_mm: float
    potential_transpiration_mm: float


@dataclass(frozen=True)
class DailyWaterBalance:
    """1日ぶんの結果。すべて床面積あたり [mm]。

    「濡れた部分あたり」ではなく「床面積あたり」に戻してあるので、
    潅水量・蒸散量・排水量をそのまま足し引きして確かめられる。
    """

    date: str
    irrigation_mm: float                  # 入れた水
    potential_transpiration_mm: float     # 気象から決まる蒸散要求 E_pmax
    actual_transpiration_mm: float        # 実際の蒸散 = min(E_pmax, U_p)
    soil_evaporation_mm: float            # 土壌面蒸発（白マルチ100%なら0）
    drainage_mm: float                    # 根群域より下へ抜けた水
    storage_change_mm: float              # 土に残った水の増減
    start_water_content: float            # 日の初めの体積含水率
    end_water_content: float              # 日の終わりの体積含水率
    min_water_content: float
    max_water_content: float
    start_pf: float
    end_pf: float
    min_air_filled_porosity: float        # その日いちばん苦しかったときの空気率
    hours_above_field_capacity: float     # 圃場容水量を超えていた時間
    hours_anoxic: float                   # 空気率が下限を割っていた時間
    uptake_limited_hours: float           # 土壌水分が蒸散を抑えていた時間

    @property
    def drainage_fraction(self) -> float:
        """潅水のうち地下へ流亡した割合。潅水がゼロの日は 0 を返す。"""
        if self.irrigation_mm <= 0.0:
            return 0.0
        return self.drainage_mm / self.irrigation_mm

    @property
    def water_stress(self) -> bool:
        """土壌水分が蒸散を抑えた日かどうか。"""
        return self.actual_transpiration_mm < self.potential_transpiration_mm - 1e-9

    @property
    def balance_residual_mm(self) -> float:
        """収支の閉じ具合の点検用。0 に近くなければ計算が壊れている。

        潅水 − 蒸散 − 土壌面蒸発 − 排水 − 貯留変化 = 0 のはずだ。
        """
        return (
            self.irrigation_mm
            - self.actual_transpiration_mm
            - self.soil_evaporation_mm
            - self.drainage_mm
            - self.storage_change_mm
        )


# =============================================================================
# 2. 1日の中の配分
# =============================================================================

def _daytime_weights(settings: WaterBalanceSettings) -> list[float]:
    """1日の各ステップに、潅水・蒸散をどの割合で割り振るかを返す。合計は1。

    日中の時間帯のあいだけ、正弦の半周期の形で配分する。
    朝夕は少なく、正午前後に山が来る形で、蒸散の日変化に近い。
    潅水も同じ形にする（日射比例で潅水しているため実態に近い）。
    """
    import math

    start_hour, end_hour = settings.daytime_window_hour
    substeps = settings.substeps_per_day
    hours_per_step = 24.0 / substeps

    weights = []
    for index in range(substeps):
        # そのステップの中心時刻
        hour = (index + 0.5) * hours_per_step
        if start_hour <= hour < end_hour:
            phase = (hour - start_hour) / (end_hour - start_hour)
            weights.append(math.sin(math.pi * phase))
        else:
            weights.append(0.0)

    total = sum(weights)
    if total <= 0.0:
        raise ValueError(
            f"日中の時間帯 {settings.daytime_window_hour} に1ステップも入らなかった。"
            f"substeps_per_day={settings.substeps_per_day} が小さすぎる可能性がある。"
        )
    return [w / total for w in weights]


# =============================================================================
# 3. 1日ぶんを進める
# =============================================================================

def step_one_day(
    water_content: float,
    day: DailyWaterInput,
    settings: WaterBalanceSettings,
) -> DailyWaterBalance:
    """1日ぶんの水収支を解いて、その日の結果を返す。

    引数の water_content は書き換えない（新しい値を結果の中に返す）。

    【min(E_pmax, U_p) を日単位で比べる理由】
    教科書の可能吸水速度 U_p は「1日にどれだけ吸えるか」を表す量で、
    根系の能力を日単位でとらえたものだ。これを5分刻みの瞬間の蒸散速度と
    比べてはいけない。

    蒸散は日中に集中するので、正午前後の瞬間の速度は日平均の 2.5 倍ほどに
    なる。日単位の U_p をそのまま瞬間値と比べると、水が十分あっても
    毎日ピーク時に頭を打ち、蒸散量が不当に小さく出る
    （実際、可能蒸散 3.5 mm/日・最大吸水 6 mm/日 の条件で 2.28 mm/日 に
    抑えられてしまった）。

    そこで、
      1. その日の可能蒸散の合計と、その日の可能吸水の合計を比べて
         実際の蒸散量（日量）を決める
      2. 決まった日量を、日中の時間帯に配分して土から引く
    という順で計算する。

    可能吸水は日の初めの含水率で評価する。1日のあいだの Aw の変化は
    小さく（潅水していれば 0.05 程度）、Up* はその範囲でほとんど動かない。

    【濡れ面積の扱い】
    計算はすべて「濡れた部分あたり [mm]」で行い、最後に床面積あたりへ
    戻す。潅水も蒸散も濡れた部分に集中するため、
      濡れた部分あたり ＝ 床面積あたり ÷ 濡れ面積率
    で割り戻してから計算に入れる。
    """
    fraction = settings.wetted_fraction
    depth_mm = settings.root_zone_depth_mm
    substeps = settings.substeps_per_day
    days_per_step = 1.0 / substeps
    weights = _daytime_weights(settings)

    # --- 床面積あたり → 濡れた部分あたり に割り戻す ---
    irrigation_wet_mm = day.irrigation_mm / fraction
    potential_transpiration_wet_mm = day.potential_transpiration_mm / fraction
    max_uptake_wet_mm_per_day = settings.max_uptake_mm_per_day / fraction
    evaporation_wet_mm = settings.soil_evaporation_mm_per_day / fraction

    theta = water_content
    theta_start = theta
    theta_min = theta
    theta_max = theta

    # --- その日の蒸散量を先に決める（式9.19の下: 蒸散 = min(E_pmax, U_p)）---
    uptake_capacity_wet_mm = (
        max_uptake_wet_mm_per_day * uptake_limit_ratio_from_water_content(theta_start)
    )
    actual_transpiration_wet_mm = min(
        potential_transpiration_wet_mm, uptake_capacity_wet_mm
    )
    is_uptake_limited = (
        uptake_capacity_wet_mm < potential_transpiration_wet_mm - 1e-12
    )

    drained_wet_mm = 0.0
    transpired_wet_mm = 0.0
    evaporated_wet_mm = 0.0

    hours_above_fc = 0.0
    hours_anoxic = 0.0
    hours_uptake_limited = 0.0
    hours_per_step = 24.0 / substeps

    theta_residual = SOIL_RETENTION["THETA_R"]

    for weight in weights:
        # --- (1) 潅水を入れる ---
        theta += irrigation_wet_mm * weight / depth_mm

        # 飽和を超えたぶんは地表にたまらず即座に抜ける（＝地表流亡）
        theta_saturated = SOIL_RETENTION["THETA_S"]
        if theta > theta_saturated:
            drained_wet_mm += (theta - theta_saturated) * depth_mm
            theta = theta_saturated

        # --- (2) 蒸散: 日量として決めた分を、日中の形に沿って引く ---
        actual_mm = actual_transpiration_wet_mm * weight
        if is_uptake_limited and weight > 0.0:
            hours_uptake_limited += hours_per_step

        # --- (3) 土壌面蒸発（白マルチ100%なら0）---
        evaporation_mm = evaporation_wet_mm * days_per_step

        # --- (4) 排水: 単位勾配の仮定で K(θ) をそのまま流束とする ---
        drainage_mm = drainage_flux_mm_per_day(theta) * days_per_step

        # --- (5) 出ていく水の合計が、残留含水率を割らないように抑える ---
        available_mm = (theta - theta_residual) * depth_mm
        outflow_mm = actual_mm + evaporation_mm + drainage_mm
        if outflow_mm > available_mm:
            # 出せる量まで比例配分で縮める。無理に引くと含水率が
            # 物理的にありえない値になり、水分特性曲線が解けなくなる。
            scale = available_mm / outflow_mm if outflow_mm > 0.0 else 0.0
            actual_mm *= scale
            evaporation_mm *= scale
            drainage_mm *= scale
            outflow_mm = available_mm

        theta -= outflow_mm / depth_mm

        transpired_wet_mm += actual_mm
        evaporated_wet_mm += evaporation_mm
        drained_wet_mm += drainage_mm

        # --- (6) その時点の土壌の状態を記録する ---
        theta_min = min(theta_min, theta)
        theta_max = max(theta_max, theta)
        if theta > FIELD_CAPACITY_WATER_CONTENT:
            hours_above_fc += hours_per_step
        if air_filled_porosity(theta) < MIN_AIR_FILLED_POROSITY:
            hours_anoxic += hours_per_step

    # --- 濡れた部分あたり → 床面積あたり に戻す ---
    start_state = describe_soil_state(theta_start)
    end_state = describe_soil_state(theta)

    return DailyWaterBalance(
        date=day.date,
        irrigation_mm=day.irrigation_mm,
        potential_transpiration_mm=day.potential_transpiration_mm,
        actual_transpiration_mm=transpired_wet_mm * fraction,
        soil_evaporation_mm=evaporated_wet_mm * fraction,
        drainage_mm=drained_wet_mm * fraction,
        storage_change_mm=(theta - theta_start) * depth_mm * fraction,
        start_water_content=theta_start,
        end_water_content=theta,
        min_water_content=theta_min,
        max_water_content=theta_max,
        start_pf=start_state.pf,
        end_pf=end_state.pf,
        min_air_filled_porosity=air_filled_porosity(theta_max),
        hours_above_field_capacity=hours_above_fc,
        hours_anoxic=hours_anoxic,
        uptake_limited_hours=hours_uptake_limited,
    )


def simulate_water_balance(
    days: Sequence[DailyWaterInput],
    settings: WaterBalanceSettings | None = None,
    initial_water_content: float | None = None,
) -> list[DailyWaterBalance]:
    """複数日ぶんの水収支を通して計算する。

    initial_water_content を省いたら、圃場容水量から始める。
    """
    if settings is None:
        settings = WaterBalanceSettings()
    theta = (
        FIELD_CAPACITY_WATER_CONTENT
        if initial_water_content is None
        else initial_water_content
    )

    results: list[DailyWaterBalance] = []
    for day in days:
        result = step_one_day(theta, day, settings)
        results.append(result)
        theta = result.end_water_content
    return results


@dataclass(frozen=True)
class SteadyState:
    """同じ天気がつづいたときの落ち着き先。"""

    #: 落ち着いたあとの1日ぶんの収支
    balance: DailyWaterBalance

    #: 落ち着くまでにかかった日数
    days: int

    #: 落ち着いたかどうか。False なら日数の上限まで動き続けている
    #: （潅水が足りずに乾きつづける、など）
    settled: bool


def steady_state(
    irrigation_mm: float,
    potential_transpiration_mm: float,
    settings: WaterBalanceSettings | None = None,
    initial_water_content: float | None = None,
    max_days: int = 120,
    tolerance: float = 1.0e-5,
) -> SteadyState:
    """同じ潅水・同じ蒸散を毎日くり返したときの、落ち着き先を求める。

    Args:
        irrigation_mm: 1日の潅水量 [mm]
        potential_transpiration_mm: 1日の可能蒸散量 [mm]
        settings: 土壌側の設定
        initial_water_content: 初日の朝の体積含水率。省くと圃場容水量から。
        max_days: 何日まで回すか
        tolerance: 朝の含水率の変化がこれを下回ったら落ち着いたとみなす

    Returns:
        SteadyState。settled が False なら落ち着いていない
        （潅水が蒸散に足りず、乾きつづけている場合など）。

    【なぜこれが要るか】
    1日だけ計算すると、潅水のうち蒸散に使われなかったぶんは
    「土に溜まった」ことになって流亡に出てこない。
    実際には翌日も翌々日も同じだけ入るので、土はやがて満杯になり、
    入れた水がそのまま抜けるようになる。この落ち着き先が実態に近い。

    たとえば5月に 3 L/m²/10MJ で潅水すると、1日だけの計算では流亡率 3% だが、
    3日目から抜け始めて定常では 53% になる。
    実測から逆算した流亡率 55〜57% と合うのはこちらのほう。
    """
    if irrigation_mm < 0.0:
        raise ValueError(f"潅水量が負: {irrigation_mm} mm")
    if potential_transpiration_mm < 0.0:
        raise ValueError(f"可能蒸散量が負: {potential_transpiration_mm} mm")
    if max_days < 1:
        raise ValueError(f"max_days は1以上で指定する。渡された値: {max_days}")

    if settings is None:
        settings = WaterBalanceSettings()
    theta = (
        FIELD_CAPACITY_WATER_CONTENT
        if initial_water_content is None
        else initial_water_content
    )

    day_input = DailyWaterInput(
        date="定常",
        irrigation_mm=irrigation_mm,
        potential_transpiration_mm=potential_transpiration_mm,
    )

    result = step_one_day(theta, day_input, settings)
    for day in range(1, max_days + 1):
        result = step_one_day(theta, day_input, settings)
        if abs(result.end_water_content - theta) < tolerance:
            return SteadyState(balance=result, days=day, settled=True)
        theta = result.end_water_content
    return SteadyState(balance=result, days=max_days, settled=False)


# =============================================================================
# 4. まとめ
# =============================================================================

@dataclass(frozen=True)
class WaterBalanceSummary:
    """期間をまとめた結果。すべて床面積あたり [mm]。"""

    days: int
    irrigation_mm: float
    potential_transpiration_mm: float
    actual_transpiration_mm: float
    soil_evaporation_mm: float
    drainage_mm: float
    storage_change_mm: float
    stressed_days: int
    anoxic_hours: float
    hours_above_field_capacity: float

    @property
    def drainage_fraction(self) -> float:
        """潅水のうち地下へ流亡した割合。"""
        if self.irrigation_mm <= 0.0:
            return 0.0
        return self.drainage_mm / self.irrigation_mm

    @property
    def transpiration_fraction(self) -> float:
        """潅水のうち蒸散に使われた割合。

        ※1を超えることがある。土に元からあった水を使った場合や、
        期間の初めが湿っていた場合で、おかしくはない。
        """
        if self.irrigation_mm <= 0.0:
            return 0.0
        return self.actual_transpiration_mm / self.irrigation_mm

    @property
    def stress_ratio(self) -> float:
        """可能蒸散に対して実際の蒸散がどれだけ届いたか（1なら無制限）。"""
        if self.potential_transpiration_mm <= 0.0:
            return 1.0
        return self.actual_transpiration_mm / self.potential_transpiration_mm


def summarize(results: Iterable[DailyWaterBalance]) -> WaterBalanceSummary:
    """日ごとの結果を期間でまとめる。"""
    rows = list(results)
    if not rows:
        raise ValueError("まとめる日がひとつもない。days が空でないか確かめることだ。")

    return WaterBalanceSummary(
        days=len(rows),
        irrigation_mm=sum(r.irrigation_mm for r in rows),
        potential_transpiration_mm=sum(r.potential_transpiration_mm for r in rows),
        actual_transpiration_mm=sum(r.actual_transpiration_mm for r in rows),
        soil_evaporation_mm=sum(r.soil_evaporation_mm for r in rows),
        drainage_mm=sum(r.drainage_mm for r in rows),
        storage_change_mm=sum(r.storage_change_mm for r in rows),
        stressed_days=sum(1 for r in rows if r.water_stress),
        anoxic_hours=sum(r.hours_anoxic for r in rows),
        hours_above_field_capacity=sum(r.hours_above_field_capacity for r in rows),
    )
