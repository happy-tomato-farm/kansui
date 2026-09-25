"""
葉のエネルギー収支を解いて、葉温を求める。

【なぜ必要か】
これまでのモデルは「葉温＝気温」と置いていた。これは光合成.xlsx が
演習用に葉温を入力として与える設計だったためで、その枠内では正しい。

しかし実データを流し込んで1日の蒸散量を出すと、日射のうち蒸散に
使われるエネルギーが18%にしかならず、物理的にありえない結果になった
（よく茂った群落なら50〜80%）。

原因は、葉温を固定すると自己調節が働かないこと。本来なら
「蒸散が足りない → 葉が冷えない → 葉温が上がる → 葉内の飽和水蒸気圧が
上がる → 蒸散が増える」という釣り合いが働く。葉温を気温に固定すると
この経路が切れてしまう。

【解く式】
葉が受け取るエネルギーと、捨てるエネルギーが釣り合う葉温を探す。

    吸収する放射 ＝ 放出する長波放射 ＋ 顕熱 ＋ 潜熱

      吸収する放射   : 日射の吸収 ＋ 周囲からの長波放射
      放出する長波   : ε σ Tl⁴ （葉の温度の4乗に比例）
      顕熱 H         : cp · gHa · (Tl − Ta)  … 葉と空気の温度差で決まる
      潜熱 λE        : λ · gv · (Cvs(Tl) − Cva) … 蒸散に伴う気化熱

葉温 Tl が上がると、右辺の3項すべてが増える（単調増加）。
したがって釣り合う点はただ1つに決まり、二分法で確実に求まる。

これは Penman-Monteith 式と同じ考え方だが、蒸散量を陽に解くのではなく
葉温を数値的に探す形をとっている。気孔コンダクタンスが日射とCO2の
経験式で与えられているため、こちらのほうが素直に書ける。

【前提条件】
1. 葉の周囲の長波放射は気温と同じ温度の面から来るとみなす。
   ハウス内は被覆・地面・他の葉に囲まれているため、この近似は妥当。
   屋外の晴天夜間のように空が冷たい場合は成り立たない。
2. 葉の熱容量を無視している（定常状態を仮定）。
   5分間隔の計算では、葉は数十秒で新しい平衡温度に達するため問題ない。
3. 気孔コンダクタンスは日射とCO2だけで決まり、葉温には依存しないとする。
   実際には高温で気孔が閉じることがあるが、経験式に含まれていない。
"""

from __future__ import annotations

from dataclasses import dataclass

import math

from config import (
    BOUNDARY_LAYER_COEF,
    HEAT_BOUNDARY_LAYER_COEF,
    LEAF_EMISSIVITY,
    LEAF_SIDES_FOR_HEAT,
    LEAF_SOLAR_ABSORPTIVITY,
    MOLAR_HEAT_CAPACITY_AIR_J_PER_MOL_K,
    MOLAR_LATENT_HEAT_WATER_J_PER_MOL,
    STEFAN_BOLTZMANN,
    STOMATA_SCALE_FACTOR,
)
from core.leaf_model import combine_in_series, stomatal_conductance_water
from core.psychrometry import saturation_vapor_pressure_kpa

# 二分法で葉温を探す範囲（気温からの差 [℃]）。
# 蒸散が盛んなときは気温より下がり、日射が強く気孔が閉じていると上がる。
SEARCH_RANGE_BELOW_AIR_C = 15.0
SEARCH_RANGE_ABOVE_AIR_C = 30.0

# 二分法の収束条件 [℃]。0.001℃まで追い込めば十分。
TEMPERATURE_TOLERANCE_C = 0.001

# 二分法の最大繰り返し回数。範囲45℃を0.001℃まで狭めるには
# log2(45/0.001) ≒ 16 回で足りるが、余裕をみる。
MAX_ITERATIONS = 60


@dataclass(frozen=True)
class EnergyBalance:
    """葉のエネルギー収支の内訳。単位はすべて W/m²（葉面積あたり）。"""

    leaf_temp_c: float          # 求まった葉温 [℃]
    air_temp_c: float           # 気温 [℃]
    absorbed_solar: float       # 吸収した日射
    absorbed_longwave: float    # 周囲から受けた長波放射
    emitted_longwave: float     # 葉が放出した長波放射
    sensible_heat: float        # 顕熱（正なら葉から空気へ）
    latent_heat: float          # 潜熱（蒸散によって奪われる熱）
    net_radiation: float        # 正味放射（吸収 − 放出）
    residual: float             # 収支の残差（0に近いほど正しく解けている）

    @property
    def leaf_air_temp_difference_c(self) -> float:
        """葉温と気温の差。正なら葉のほうが暖かい。"""
        return self.leaf_temp_c - self.air_temp_c

    @property
    def latent_heat_fraction(self) -> float:
        """正味放射のうち潜熱（蒸散）に使われた割合。"""
        if self.net_radiation <= 0:
            return 0.0
        return self.latent_heat / self.net_radiation


def heat_boundary_layer_conductance(
    wind_speed_m_per_s: float,
    characteristic_length_m: float,
) -> float:
    """
    熱の境界層コンダクタンスを求める。

        gHa = 0.135 * sqrt(u / d) * （両面ぶん）

    水蒸気の係数（0.147）と違うのは、熱と水蒸気で拡散のしかたが
    異なるため。また顕熱は葉の両面から出入りするので2倍する。

    Returns:
        熱の境界層コンダクタンス [mol/(m²·s)]
    """
    if characteristic_length_m <= 0:
        raise ValueError(
            f"葉の特性長が 0 以下です: {characteristic_length_m} m"
        )
    if wind_speed_m_per_s < 0:
        raise ValueError(f"風速が負です: {wind_speed_m_per_s} m/s")

    one_side = HEAT_BOUNDARY_LAYER_COEF * math.sqrt(
        wind_speed_m_per_s / characteristic_length_m
    )

    return one_side * LEAF_SIDES_FOR_HEAT


def _energy_residual(
    leaf_temp_c: float,
    air_temp_c: float,
    vapor_pressure_kpa: float,
    solar_radiation_w_per_m2: float,
    vapor_conductance: float,
    heat_conductance: float,
    pressure_kpa: float,
    longwave_view_factor: float,
    solar_absorptivity: float,
) -> tuple[float, dict]:
    """
    ある葉温を仮定したときの、エネルギー収支の残差を求める。

    残差 = 吸収 − 放出 − 顕熱 − 潜熱

    残差が 0 になる葉温が答え。正なら葉が受け取りすぎ（もっと熱くなる）、
    負なら捨てすぎ（もっと冷える）。

    Returns:
        (残差 [W/m²], 内訳の辞書)
    """
    air_temp_k = air_temp_c + 273.15
    leaf_temp_k = leaf_temp_c + 273.15

    # --- 吸収する放射 --------------------------------------------------
    absorbed_solar = solar_absorptivity * solar_radiation_w_per_m2

    # --- 長波放射の正味交換 ---------------------------------------------
    # 葉は周囲と長波放射をやりとりする。周囲には2種類ある。
    #
    #   (1) 外界（ハウスの被覆・地面・空）… 温度は気温とみなす
    #   (2) 周囲の葉                      … 温度はその葉自身とほぼ同じ
    #
    # 群落の奥にある葉は、周りを同じような温度の葉に囲まれているため、
    # 長波をやりとりしても差し引きゼロに近い。正味の損失が生じるのは
    # 外界が見える群落の表面付近だけである。
    #
    # そこで「外界がどれだけ見えるか」を表す視野率を掛ける。
    # 視野率は群落の表面で1に近く、奥へ行くほど0に近づく。
    #
    # 層ごとに損失を積み上げていた以前の実装は、この点を無視していたため
    # 長波の損失を大きく見積もりすぎていた（LAI倍になっていた）。
    absorbed_longwave = (
        LEAF_SIDES_FOR_HEAT * LEAF_EMISSIVITY * STEFAN_BOLTZMANN * air_temp_k ** 4
    )
    emitted_longwave = (
        LEAF_SIDES_FOR_HEAT * LEAF_EMISSIVITY * STEFAN_BOLTZMANN * leaf_temp_k ** 4
    )

    # 外界が見えない分は、周囲の葉と交換して差し引きゼロになる。
    net_longwave = longwave_view_factor * (absorbed_longwave - emitted_longwave)

    # --- 顕熱（空気への熱伝達）------------------------------------------
    sensible_heat = (
        MOLAR_HEAT_CAPACITY_AIR_J_PER_MOL_K
        * heat_conductance
        * (leaf_temp_c - air_temp_c)
    )

    # --- 潜熱（蒸散による気化熱）----------------------------------------
    # 葉の内部は飽和しているとみなす。葉温が高いほど葉内の水蒸気圧が
    # 上がり、蒸散が増える。ここが自己調節の要になる。
    leaf_saturation_kpa = saturation_vapor_pressure_kpa(leaf_temp_c)
    vapor_gradient = (leaf_saturation_kpa - vapor_pressure_kpa) / pressure_kpa
    transpiration_mol = vapor_conductance * vapor_gradient
    latent_heat = MOLAR_LATENT_HEAT_WATER_J_PER_MOL * transpiration_mol

    residual = absorbed_solar + net_longwave - sensible_heat - latent_heat

    breakdown = {
        "absorbed_solar": absorbed_solar,
        "absorbed_longwave": longwave_view_factor * absorbed_longwave,
        "emitted_longwave": longwave_view_factor * emitted_longwave,
        "net_longwave": net_longwave,
        "sensible_heat": sensible_heat,
        "latent_heat": latent_heat,
        "transpiration_mol": transpiration_mol,
    }

    return residual, breakdown


def solve_leaf_temperature(
    air_temp_c: float,
    vapor_pressure_kpa: float,
    solar_radiation_w_per_m2: float,
    air_co2_ppm: float,
    wind_speed_m_per_s: float,
    pressure_kpa: float,
    characteristic_length_m: float,
    longwave_view_factor: float = 1.0,
    solar_absorptivity: float = LEAF_SOLAR_ABSORPTIVITY,
    stomata_scale_factor: float = STOMATA_SCALE_FACTOR,
) -> EnergyBalance:
    """
    エネルギー収支が釣り合う葉温を求める。

    【解き方】
    葉温を上げると「放出する長波＋顕熱＋潜熱」はすべて増える。
    つまり残差は葉温について単調減少する。
    したがって残差が正になる温度と負になる温度で挟み、
    二分法で挟み撃ちにすれば必ず解が求まる。

    Args:
        air_temp_c: 気温 [℃]
        vapor_pressure_kpa: 大気の水蒸気圧 [kPa]
        solar_radiation_w_per_m2: その葉が受ける日射 [W/m²]
        air_co2_ppm: CO2濃度 [ppm]（気孔コンダクタンスの計算に使う）
        wind_speed_m_per_s: 風速 [m/s]
        pressure_kpa: 大気圧 [kPa]
        characteristic_length_m: 葉の特性長 [m]
        longwave_view_factor: 外界がどれだけ見えるか（0〜1）。
            1 なら孤立した葉。群落の中では小さくなる。
            canopy.py が層ごとに計算して渡す。
        solar_absorptivity: 葉面積あたりの日射吸収率。
            群落の中では、Beer則の消散係数から決まる値を使う。

    Returns:
        EnergyBalance（葉温と収支の内訳）
    """
    if not 0.0 <= longwave_view_factor <= 1.0:
        raise ValueError(
            f"長波の視野率は 0〜1 で指定してください: {longwave_view_factor}"
        )
    # 気孔コンダクタンスは日射とCO2で決まり、葉温には依存しない。
    # ループの外で1回だけ計算すればよい。
    par_w_per_m2 = solar_radiation_w_per_m2 / 2.0
    gvs = stomatal_conductance_water(par_w_per_m2, air_co2_ppm, stomata_scale_factor)
    gva = BOUNDARY_LAYER_COEF * math.sqrt(
        wind_speed_m_per_s / characteristic_length_m
    )
    vapor_conductance = combine_in_series(gvs, gva)

    heat_conductance = heat_boundary_layer_conductance(
        wind_speed_m_per_s, characteristic_length_m
    )

    def residual_at(temp_c: float) -> float:
        value, _ = _energy_residual(
            temp_c, air_temp_c, vapor_pressure_kpa, solar_radiation_w_per_m2,
            vapor_conductance, heat_conductance, pressure_kpa,
            longwave_view_factor, solar_absorptivity,
        )
        return value

    low = air_temp_c - SEARCH_RANGE_BELOW_AIR_C
    high = air_temp_c + SEARCH_RANGE_ABOVE_AIR_C

    residual_low = residual_at(low)
    residual_high = residual_at(high)

    # 探索範囲の内側に解があることを確かめる。
    # 範囲外なら、そのまま二分法を回しても誤った答えになるため、
    # 黙って続けず端の値を返して知らせる。
    if residual_low < 0:
        # 下限でも「捨てすぎ」＝もっと冷える方向。ありえないほど蒸散が大きい。
        raise ValueError(
            f"葉温が探索範囲の下限を下回りました（気温 {air_temp_c}℃ の"
            f"{SEARCH_RANGE_BELOW_AIR_C}℃ 下）。"
            f"水蒸気圧 {vapor_pressure_kpa} kPa、日射 {solar_radiation_w_per_m2} W/m² を"
            f"確認してください"
        )
    if residual_high > 0:
        raise ValueError(
            f"葉温が探索範囲の上限を超えました（気温 {air_temp_c}℃ の"
            f"{SEARCH_RANGE_ABOVE_AIR_C}℃ 上）。"
            f"日射 {solar_radiation_w_per_m2} W/m² が大きすぎる可能性があります"
        )

    # 二分法
    for _ in range(MAX_ITERATIONS):
        middle = (low + high) / 2.0
        if high - low < TEMPERATURE_TOLERANCE_C:
            break

        if residual_at(middle) > 0:
            # まだ受け取りすぎ → もっと熱い側に解がある
            low = middle
        else:
            high = middle

    leaf_temp_c = (low + high) / 2.0
    residual, breakdown = _energy_residual(
        leaf_temp_c, air_temp_c, vapor_pressure_kpa, solar_radiation_w_per_m2,
        vapor_conductance, heat_conductance, pressure_kpa,
        longwave_view_factor, solar_absorptivity,
    )

    return EnergyBalance(
        leaf_temp_c=leaf_temp_c,
        air_temp_c=air_temp_c,
        absorbed_solar=breakdown["absorbed_solar"],
        absorbed_longwave=breakdown["absorbed_longwave"],
        emitted_longwave=breakdown["emitted_longwave"],
        sensible_heat=breakdown["sensible_heat"],
        latent_heat=breakdown["latent_heat"],
        net_radiation=breakdown["absorbed_solar"] + breakdown["net_longwave"],
        residual=residual,
    )


def longwave_view_factor_in_canopy(
    cumulative_lai_above: float,
    cumulative_lai_below: float,
    extinction_coefficient: float,
) -> float:
    """
    群落内のある層から、外界がどれだけ見えるかを求める。

    【考え方】
    葉は上向きと下向きの両方に長波放射を出す。
    上向きには、その上にある葉を通り抜けた分だけ外界（被覆・空）が見える。
    下向きには、その下にある葉を通り抜けた分だけ外界（地面）が見える。

        上向きの視野率 = exp(-k × 上にある葉面積)
        下向きの視野率 = exp(-k × 下にある葉面積)

    両面の平均をとって、その層の視野率とする。
    群落の表面なら1に近く、奥に行くほど0に近づく。

    外界が見えない分は、周囲の葉（自分とほぼ同じ温度）と長波を交換するので
    差し引きゼロになる。

    Args:
        cumulative_lai_above: その層より上にある葉面積指数
        cumulative_lai_below: その層より下にある葉面積指数
        extinction_coefficient: 消散係数 k

    Returns:
        視野率（0〜1）
    """
    upward = math.exp(-extinction_coefficient * cumulative_lai_above)
    downward = math.exp(-extinction_coefficient * cumulative_lai_below)

    return (upward + downward) / 2.0
