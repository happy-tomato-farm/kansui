"""
葉1枚のモデルを群落全体（地面1m²あたり）に広げる。

【なぜスケーリングが必要か】
leaf_model.py が計算するのは「葉面積1m²あたり」の値。
知りたいのは「地面1m²あたり」「ハウス1棟あたり」の値なので、
葉がどれだけあるか（LAI）を掛ける必要がある。

ただし単純に LAI を掛けるだけでは間違いになる。
群落の上のほうの葉は日射を強く受けるが、下のほうの葉は上の葉に
遮られて暗い。暗い葉は光合成も蒸散も少ない。
そこで群落を層に分け、層ごとに日射を減衰させて計算し、足し合わせる。

【Beer則（ビールの法則）】
群落内の日射は、上からの累積葉面積に対して指数関数的に減る。

    I(L) = I0 * exp(-k * L)

        I0 : 群落の上の日射 [W/m²]
        L  : その高さまでの累積葉面積指数 [m²/m²]
        k  : 消散係数（葉の傾きや配置で決まる。トマトでは 0.6〜0.8）

LAI=2.2、k=0.7 なら、最下層に届く光は exp(-0.7*2.2) = 21%。
上の葉と下の葉で5倍近い差がある。

【前提条件（崩れると結果がずれる）】
1. 群落内の気温・湿度・CO2・風速は一様と仮定している。
   実際には上下で勾配があるが、実測がないため一様として扱う。
   センサーの設置高さが分かれば、その高さの値として解釈することになる。
2. 直達光と散乱光を分けていない（sunlit/shaded を区別しない）。
   葉の隙間から差し込む光の斑を無視するため、光合成をやや過小に
   見積もる傾向がある。ただしハウス内の日射は被覆を通った時点で
   散乱成分が多いため、この影響は屋外ほど大きくない。
3. すべての葉が同じ性能を持つと仮定している。
   実際には下位の老化葉は光合成能力が落ちる。
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass

import numpy as np

from config import (
    EXTINCTION_COEFFICIENT_K,
    LEAF_AREA_PER_LEAF_M2,
    LEAF_COUNT_OUT_OF_SEASON_PER_M2,
    LEAF_COUNT_PLAN_PER_M2,
    LEAF_SOLAR_ABSORPTIVITY,
    MOLAR_MASS_CO2_G_PER_MOL,
    MOLAR_MASS_WATER_G_PER_MOL,
)
from core.energy_balance import (
    heat_boundary_layer_conductance,
    longwave_view_factor_in_canopy,
)
from core.leaf_model import calculate_leaf_gas_exchange

# 群落を何層に分けて計算するか。
# 層を増やすほど正確になるが計算時間が増える。
# 20層あれば実用上じゅうぶん収束する（tests/test_canopy.py で確認）。
DEFAULT_LAYERS = 20


@dataclass(frozen=True)
class CanopyGasExchange:
    """
    群落（地面1m²あたり）のガス交換の計算結果。

    単位に注意。leaf_model は「葉面積1m²あたり」、
    こちらは「地面1m²あたり」。
    """

    lai: float                              # 葉面積指数 [m²/m²]
    transpiration_mg_per_m2_s: float        # 蒸散速度 [mg/(m²·s)] 地面あたり
    net_photosynthesis_mmol_per_m2_s: float # 純光合成速度 [mmol/(m²·s)] 地面あたり
    gross_photosynthesis_mmol_per_m2_s: float  # 総光合成速度 [mmol/(m²·s)]
    respiration_mmol_per_m2_s: float        # 呼吸（光呼吸＋暗呼吸）[mmol/(m²·s)]
    sunlit_fraction: float                  # 上層が受ける日射に対する群落平均の割合
    top_layer_solar_w_per_m2: float         # 最上層に届く日射 [W/m²]
    bottom_layer_solar_w_per_m2: float      # 最下層に届く日射 [W/m²]
    mean_leaf_temp_c: float                 # 葉面積で重みづけした平均葉温 [℃]
    top_layer_leaf_temp_c: float            # 最上層の葉温 [℃]（最も高温になる）

    # --- 群落を「1枚の大きな葉」とみなしたときの量 ---
    # Penman-Monteith 式（教科書 式14.12）と突き合わせるために必要。
    # 各層の値を葉面積で重みづけして足し上げたもの。
    canopy_conductance_water: float         # 群落の水蒸気コンダクタンス gv [mol/(m²·s)]
    canopy_conductance_heat: float          # 群落の熱コンダクタンス gHa [mol/(m²·s)]
    absorbed_radiation_w_per_m2: float      # 群落が吸収した日射 [W/m²]
    below_canopy_solar_w_per_m2: float      # 群落を通り抜けた日射 [W/m²]


def leaf_area_index(
    leaves_per_m2: float,
    leaf_area_per_leaf_m2: float = LEAF_AREA_PER_LEAF_M2,
) -> float:
    """
    葉枚数から葉面積指数（LAI）を求める。

    【LAIとは】
    Leaf Area Index。地面1m²の上に葉が何m²あるかという値。
    LAI=2 なら、地面を2枚重ねで覆うだけの葉があるということ。

        LAI = 葉枚数 [枚/m²] × 葉1枚の面積 [m²/枚]

    葉枚数は季節で大きく変わる。冬 16 枚/m² から春 40 枚/m² まで倍以上。
    日付から引くには leaf_count_per_m2() を使う。
    葉1枚の面積は現在は暫定値で、実測後に
    config.LEAF_AREA_PER_LEAF_M2 を書き換える。

    Args:
        leaves_per_m2: 地面1m²あたりの葉枚数 [枚/m²]
        leaf_area_per_leaf_m2: 葉1枚の面積 [m²/枚]

    Returns:
        葉面積指数 LAI [m²/m²]
    """
    if leaves_per_m2 < 0:
        raise ValueError(f"葉枚数が負です: {leaves_per_m2} 枚/m²")
    if leaf_area_per_leaf_m2 <= 0:
        raise ValueError(
            f"葉1枚の面積が 0 以下です: {leaf_area_per_leaf_m2} m²/枚。"
            f"config.LEAF_AREA_PER_LEAF_M2 を確認してください"
        )

    return leaves_per_m2 * leaf_area_per_leaf_m2


def leaf_count_per_m2(date: dt.date) -> float:
    """作業計画（週ごとの目標葉枚数）から、その日の葉枚数を求める [枚/m²]。

    【なぜ日付から引くのか】
    葉枚数は季節で倍以上変わる。冬 16 枚/m²、春 40 枚/m²。
    年間固定にすると春の蒸散を4割ほど過小に見積もる。

    【週の切れ目をなめらかにつなぐ】
    計画は週単位だが、そのまま使うと週が変わる日に値が飛ぶ。
    週の中心どうしを直線で結んで、日ごとに少しずつ変わるようにしている。

    【栽培していない時期】
    週27〜37（7月上旬〜9月中旬）は作期の外なので計画値がない。
    その期間は config.LEAF_COUNT_OUT_OF_SEASON_PER_M2 を返す。

    Args:
        date: 日付

    Returns:
        葉枚数 [枚/m²]
    """
    year, week, weekday = date.isocalendar()
    if week not in LEAF_COUNT_PLAN_PER_M2:
        return LEAF_COUNT_OUT_OF_SEASON_PER_M2

    current = LEAF_COUNT_PLAN_PER_M2[week]

    # 週の前半は前の週へ、後半は次の週へ向かって直線で補間する。
    # weekday は月曜1〜日曜7。週の中心を木曜（4）とみなす。
    offset = (weekday - 4) / 7.0
    neighbour_week = week + (1 if offset > 0 else -1)
    # 週53の次は週1、週1の前は週52（年またぎ）
    if neighbour_week > 53:
        neighbour_week = 1
    elif neighbour_week < 1:
        neighbour_week = 52

    neighbour = LEAF_COUNT_PLAN_PER_M2.get(neighbour_week)
    if neighbour is None:
        return current
    return current + (neighbour - current) * abs(offset)


def leaf_area_index_for_date(
    date: dt.date,
    leaf_area_per_leaf_m2: float = LEAF_AREA_PER_LEAF_M2,
) -> float:
    """日付から LAI を求める。作業計画の葉枚数 × 葉1枚の面積。"""
    return leaf_area_index(leaf_count_per_m2(date), leaf_area_per_leaf_m2)


def solar_at_depth(
    solar_above_w_per_m2: float,
    cumulative_lai: float,
    extinction_coefficient: float = EXTINCTION_COEFFICIENT_K,
) -> float:
    """
    群落内のある深さに届く日射を Beer則で求める。

        I(L) = I0 * exp(-k * L)

    Args:
        solar_above_w_per_m2: 群落上部の日射 I0 [W/m²]
        cumulative_lai: その深さまでの累積葉面積指数 L [m²/m²]
        extinction_coefficient: 消散係数 k

    Returns:
        その深さの日射 [W/m²]
    """
    if cumulative_lai < 0:
        raise ValueError(f"累積葉面積指数が負です: {cumulative_lai}")

    return solar_above_w_per_m2 * np.exp(-extinction_coefficient * cumulative_lai)


def calculate_canopy_gas_exchange(
    air_temp_c: float,
    vapor_pressure_kpa: float,
    solar_above_w_per_m2: float,
    air_co2_ppm: float,
    wind_speed_m_per_s: float,
    pressure_kpa: float,
    characteristic_length_m: float,
    lai: float,
    extinction_coefficient: float = EXTINCTION_COEFFICIENT_K,
    n_layers: int = DEFAULT_LAYERS,
    leaf_temp_c: float | None = None,
    solve_energy_balance: bool = False,
) -> CanopyGasExchange:
    """
    群落全体（地面1m²あたり）の蒸散速度と光合成速度を求める。

    【計算の流れ】
    1. 群落を n_layers 層に等分する（各層の葉面積は LAI/n）
    2. 各層の中央に届く日射を Beer則で求める
    3. その日射で葉1枚あたりのガス交換を計算する
    4. 各層の値に層の葉面積 LAI/n を掛けて足し合わせる

    Args:
        air_temp_c: 気温 [℃]
        vapor_pressure_kpa: 水蒸気圧 [kPa]
        solar_above_w_per_m2: 群落上部の日射 [W/m²]
        air_co2_ppm: CO2濃度 [ppm]
        wind_speed_m_per_s: 風速 [m/s]
        pressure_kpa: 大気圧 [kPa]
        characteristic_length_m: 葉の特性長 [m]
        lai: 葉面積指数 [m²/m²]
        extinction_coefficient: 消散係数 k
        n_layers: 群落を分割する層の数
        leaf_temp_c: 葉温 [℃]。指定するとすべての層でその値を使う。
        solve_energy_balance: True なら層ごとに葉温をエネルギー収支から求める。
            上層は日射を強く受けるので葉温が高く、下層は低くなる。

    Returns:
        CanopyGasExchange（地面1m²あたりの値）
    """
    if lai < 0:
        raise ValueError(f"葉面積指数が負です: {lai}")
    if n_layers < 1:
        raise ValueError(f"層の数は1以上にしてください: {n_layers}")

    # 葉がなければガス交換も起きない。
    if lai == 0.0:
        return CanopyGasExchange(
            lai=0.0,
            transpiration_mg_per_m2_s=0.0,
            net_photosynthesis_mmol_per_m2_s=0.0,
            gross_photosynthesis_mmol_per_m2_s=0.0,
            respiration_mmol_per_m2_s=0.0,
            sunlit_fraction=0.0,
            top_layer_solar_w_per_m2=solar_above_w_per_m2,
            bottom_layer_solar_w_per_m2=solar_above_w_per_m2,
            mean_leaf_temp_c=air_temp_c,
            top_layer_leaf_temp_c=air_temp_c,
            canopy_conductance_water=0.0,
            canopy_conductance_heat=0.0,
            absorbed_radiation_w_per_m2=0.0,
            below_canopy_solar_w_per_m2=solar_above_w_per_m2,
        )

    layer_lai = lai / n_layers  # 1層あたりの葉面積指数

    total_transpiration = 0.0
    total_net = 0.0
    total_gross = 0.0
    total_respiration = 0.0
    weighted_solar = 0.0
    weighted_leaf_temp = 0.0
    top_leaf_temp = air_temp_c
    total_conductance_water = 0.0
    total_conductance_heat = 0.0
    total_absorbed_radiation = 0.0

    for i in range(n_layers):
        # 層の中央の深さ（累積LAI）で代表させる。
        # 層の上端を使うと日射を過大に、下端を使うと過小に見積もるため。
        cumulative_lai = (i + 0.5) * layer_lai

        layer_solar = solar_at_depth(
            solar_above_w_per_m2, cumulative_lai, extinction_coefficient
        )

        # この層から外界がどれだけ見えるか。
        # 群落の奥ほど周囲を葉に囲まれ、長波の正味交換が小さくなる。
        view_factor = longwave_view_factor_in_canopy(
            cumulative_lai_above=cumulative_lai,
            cumulative_lai_below=lai - cumulative_lai,
            extinction_coefficient=extinction_coefficient,
        )

        leaf = calculate_leaf_gas_exchange(
            air_temp_c=air_temp_c,
            vapor_pressure_kpa=vapor_pressure_kpa,
            solar_radiation_w_per_m2=layer_solar,
            air_co2_ppm=air_co2_ppm,
            wind_speed_m_per_s=wind_speed_m_per_s,
            pressure_kpa=pressure_kpa,
            characteristic_length_m=characteristic_length_m,
            leaf_temp_c=leaf_temp_c,
            solve_energy_balance=solve_energy_balance,
            longwave_view_factor=view_factor,
        )

        # 葉面積1m²あたりの値に、この層の葉面積を掛けて地面あたりにする
        total_transpiration += leaf.transpiration_mg_per_m2_s * layer_lai
        total_net += leaf.net_photosynthesis_mmol * layer_lai
        total_gross += leaf.gross_photosynthesis_mmol * layer_lai
        total_respiration += (
            leaf.photorespiration_mmol + leaf.dark_respiration_mmol
        ) * layer_lai
        weighted_solar += layer_solar * layer_lai
        weighted_leaf_temp += leaf.leaf_temp_c * layer_lai

        # 層は互いに並列につながっているので、コンダクタンスは足し算になる。
        # （直列なら抵抗の和だが、葉は横に並んでいるので並列）
        total_conductance_water += leaf.total_conductance_water * layer_lai
        total_conductance_heat += (
            heat_boundary_layer_conductance(
                wind_speed_m_per_s, characteristic_length_m
            ) * layer_lai
        )
        # この層の葉が吸収した日射（葉面積あたり × その層の葉面積）
        total_absorbed_radiation += (
            LEAF_SOLAR_ABSORPTIVITY * layer_solar * layer_lai
        )

        if i == 0:
            top_leaf_temp = leaf.leaf_temp_c

    # 群落平均の受光量が、上部の日射の何割にあたるか
    if solar_above_w_per_m2 > 0:
        sunlit_fraction = (weighted_solar / lai) / solar_above_w_per_m2
    else:
        sunlit_fraction = 0.0

    return CanopyGasExchange(
        lai=lai,
        transpiration_mg_per_m2_s=total_transpiration,
        net_photosynthesis_mmol_per_m2_s=total_net,
        gross_photosynthesis_mmol_per_m2_s=total_gross,
        respiration_mmol_per_m2_s=total_respiration,
        sunlit_fraction=sunlit_fraction,
        top_layer_solar_w_per_m2=solar_at_depth(
            solar_above_w_per_m2, 0.5 * layer_lai, extinction_coefficient
        ),
        bottom_layer_solar_w_per_m2=solar_at_depth(
            solar_above_w_per_m2, (n_layers - 0.5) * layer_lai, extinction_coefficient
        ),
        mean_leaf_temp_c=weighted_leaf_temp / lai,
        top_layer_leaf_temp_c=top_leaf_temp,
        canopy_conductance_water=total_conductance_water,
        canopy_conductance_heat=total_conductance_heat,
        absorbed_radiation_w_per_m2=total_absorbed_radiation,
        below_canopy_solar_w_per_m2=solar_at_depth(
            solar_above_w_per_m2, lai, extinction_coefficient
        ),
    )


def transpiration_to_mm_per_s(transpiration_mg_per_m2_s: float) -> float:
    """
    蒸散速度の単位を mg/(m²·s) から mm/s（＝L/(m²·s)）に直す。

    水1mgは1mm³。地面1m²に1mm³の水は 1e-6 mm の深さ。
    つまり mg/m² を mm に直すには 1e-6 を掛ける。

    農業では蒸散量を「1日に何mm」「1m²あたり何L」で表すことが多く、
    1 mm = 1 L/m² の関係にある。
    """
    return transpiration_mg_per_m2_s * 1.0e-6


def photosynthesis_to_co2_g_per_m2_s(net_photosynthesis_mmol_per_m2_s: float) -> float:
    """
    純光合成速度を mmol/(m²·s) から CO2 の g/(m²·s) に直す。

    CO2 の分子量は 44 g/mol。
    """
    return net_photosynthesis_mmol_per_m2_s / 1000.0 * MOLAR_MASS_CO2_G_PER_MOL


def photosynthesis_to_sugar_g_per_m2_s(
    net_photosynthesis_mmol_per_m2_s: float,
) -> float:
    """
    純光合成速度を、固定した炭水化物（ブドウ糖換算）の重さに直す。

    【換算の根拠】
    光合成の全体反応は
        6 CO2 + 6 H2O → C6H12O6 + 6 O2
    なので、CO2 6分子からブドウ糖（分子量180）が1分子できる。
    CO2 1 mol あたり 180/6 = 30 g のブドウ糖に相当する。

    実際の乾物にはタンパク質や有機酸も含まれるため、これは
    「炭水化物に換算するとどれだけか」という目安の値。
    """
    GLUCOSE_G_PER_CO2_MOL = 180.0 / 6.0

    return net_photosynthesis_mmol_per_m2_s / 1000.0 * GLUCOSE_G_PER_CO2_MOL
