"""
葉1枚（正確には葉面積1m²あたり）の蒸散速度と光合成速度を求めるモデル。

【出典】
生物環境物理学/光合成.xlsx をそのまま移植したものだ。
式番号は「生物環境物理学の基礎」に対応する。

【このモデルの位置づけ】
ここで計算するのは「葉面積1m²あたり」の値であって、
「地面1m²あたり」ではない。地面あたりに直すには LAI（葉面積指数）を
掛ける必要があり、それは core/canopy.py が担当する。

【モデルの骨格】
蒸散も光合成も「濃度差 ÷ 抵抗」という同じ形をしている。

    蒸散   = 総コンダクタンス(水蒸気) × (葉内の水蒸気モル分率 - 大気の水蒸気モル分率)
    光合成 = 総コンダクタンス(CO2)   × (大気のCO2濃度 - 葉緑体のCO2濃度)

コンダクタンス（通しやすさ）は抵抗の逆数。気孔と境界層という
2つの関門が直列につながっているので、電気回路の直列抵抗と同じく
逆数の和の逆数で合成する。

【前提条件（崩れると結果がずれる）】
1. 葉温 = 気温 と仮定している。日射が強いとき実際の葉温は気温より
   数℃高くなるため、蒸散はやや過小評価になる。
2. 気孔コンダクタンスは日射とCO2濃度だけで決まる経験式を使っている。
   土壌水分の不足による気孔閉鎖は考慮していない。
3. 風速は一定値を使う。実際には天窓・サイドの開閉で変わる。
"""

import math
from dataclasses import dataclass

from config import (
    LIGHT_RESPONSE_MODEL,
    NON_RECTANGULAR_LIGHT_PARAMS,
    PAR_UMOL_PER_JOULE,
    BOUNDARY_LAYER_COEF,
    CO2_CONDUCTANCE_RATIO,
    MOLAR_MASS_CO2_G_PER_MOL,
    MOLAR_MASS_WATER_G_PER_MOL,
    PHOTOSYNTHESIS_PARAMS,
    STOMATA_PARAMS,
    STOMATA_SCALE_FACTOR,
)
from core.psychrometry import mole_fraction, saturation_vapor_pressure_kpa


@dataclass(frozen=True)
class LeafGasExchange:
    """
    葉のガス交換の計算結果。

    中間値もすべて保持しているのは、計算がおかしいときに
    どの段階でずれたかを追えるようにするためだ。
    frozen=True にして、作った後から書き換えられないようにしている
    （計算結果が知らないうちに変わる事故を防ぐ）。
    """

    # --- 水蒸気側 ---
    stomatal_conductance_water: float      # 気孔コンダクタンス gvs [mol/(m²·s)]
    boundary_conductance_water: float      # 境界層コンダクタンス gva [mol/(m²·s)]
    total_conductance_water: float         # 総コンダクタンス gv [mol/(m²·s)]
    air_vapor_mole_fraction: float         # 大気の水蒸気モル分率 Cva [mol/mol]
    leaf_vapor_mole_fraction: float        # 葉内の水蒸気モル分率 Cvs [mol/mol]
    transpiration_mol_per_m2_s: float      # 蒸散速度 Fv [mol/(m²·s)]
    transpiration_mg_per_m2_s: float       # 蒸散速度 Fv [mg/(m²·s)]

    # --- CO2側 ---
    stomatal_conductance_co2: float        # 気孔コンダクタンス gcs [mol/(m²·s)]
    boundary_conductance_co2: float        # 境界層コンダクタンス gca [mol/(m²·s)]
    total_conductance_co2: float           # 総コンダクタンス gc [mol/(m²·s)]

    # --- 光合成の各項 ---
    light_function: float                  # 光合成光関数 f(Stp) [0-1]
    temp_function: float                   # 光合成温度関数 h(Tl) [0-1]
    max_photosynthesis_mol_per_m2_s: float # CO2最大時光合成速度 Pm [mol/(m²·s)]
    gross_photosynthesis_mmol: float       # 総光合成速度 P [mmol/(m²·s)]
    photorespiration_mmol: float           # 光呼吸 Rp [mmol/(m²·s)]
    dark_respiration_mmol: float           # 暗呼吸 Rd [mmol/(m²·s)]
    net_photosynthesis_mmol: float         # 純光合成速度 Pn [mmol/(m²·s)]
    net_photosynthesis_mg: float           # 純光合成速度 Pn [mg/(m²·s)]

    # --- 葉温 ---
    leaf_temp_c: float                     # 計算に使った葉温 [℃]
    air_temp_c: float                      # 気温 [℃]

    @property
    def leaf_air_temp_difference_c(self) -> float:
        """葉温と気温の差。正なら葉のほうが暖かい。"""
        return self.leaf_temp_c - self.air_temp_c


def stomatal_conductance_water(
    par_w_per_m2: float,
    co2_ppm: float,
    scale_factor: float = STOMATA_SCALE_FACTOR,
) -> float:
    """
    気孔コンダクタンス（水蒸気）を求める。

    【理屈】
    気孔は光が当たると開き、CO2が濃いと閉じる。
    植物にとって気孔を開くのは「CO2を取り込むため」だが、
    同時に水も失う。だからCO2が十分あるなら気孔を絞る。

        gvs = gvs'（日射による開き） - gvs''（CO2による閉じ）

    【水準の補正について】
    光合成.xlsx の式は演習用の簡略式で、上限が 0.207 mol/(m²·s) と
    実測（日中 0.3〜0.6）より2〜3倍低い。式の形はそのままに、
    scale_factor で水準だけを実測に合わせる。
    根拠と出典は config.STOMATA_SCALE_FACTOR のコメントを参照。

    Args:
        par_w_per_m2: 光合成有効日射 Stp [W/m²]
        co2_ppm: 大気CO2濃度 Cca [ppm]
        scale_factor: 水準の補正倍率。1.0 なら光合成.xlsx のまま。

    Returns:
        気孔コンダクタンス gvs [mol/(m²·s)]
    """
    if scale_factor <= 0:
        raise ValueError(
            f"気孔コンダクタンスの補正倍率が 0 以下です: {scale_factor}"
        )

    p = STOMATA_PARAMS

    # 日射による開き。弱光では直線、強光では頭打ちになる折れ線近似。
    #   min(max(A*Stp + B, C*Stp), 上限)
    light_term = min(
        max(
            p["LIGHT_SLOPE_A"] * par_w_per_m2 + p["LIGHT_INTERCEPT_B"],
            p["LIGHT_SLOPE_C"] * par_w_per_m2,
        ),
        p["MAX_CONDUCTANCE"],
    )

    # CO2による閉じ。一定濃度を超えてから効き始める。
    co2_term = max(p["CO2_SLOPE_D"] * co2_ppm - p["CO2_INTERCEPT_E"], 0.0)

    conductance = (light_term - co2_term) * scale_factor

    # CO2が極端に濃いと差が負になりうる。気孔が「逆に開く」ことはないので
    # 0 で下限を切る（完全閉鎖）。
    return max(conductance, 0.0)


def boundary_layer_conductance_water(
    wind_speed_m_per_s: float,
    characteristic_length_m: float,
) -> float:
    """
    境界層コンダクタンス（水蒸気）を求める。（生物環境物理学の基礎 p.113）

    【境界層とは】
    葉の表面にはりつく、動きの遅い薄い空気の層。
    ここを水蒸気やCO2が拡散で通り抜ける必要があるため、抵抗になる。
    風が強いほど境界層は薄くなり、通しやすくなる。

        gva = 0.147 * sqrt(u / d)

    風速の平方根に比例するので、風速を 0.2 → 0.5 m/s にしても
    コンダクタンスは約1.6倍にしかならない。

    Args:
        wind_speed_m_per_s: 風速 u [m/s]
        characteristic_length_m: 葉の特性長 d [m]

    Returns:
        境界層コンダクタンス gva [mol/(m²·s)]
    """
    if wind_speed_m_per_s < 0:
        raise ValueError(f"風速が負です: {wind_speed_m_per_s} m/s")
    if characteristic_length_m <= 0:
        raise ValueError(
            f"葉の特性長が 0 以下です: {characteristic_length_m} m。"
            f"葉幅の設定を確認してください"
        )

    return BOUNDARY_LAYER_COEF * math.sqrt(wind_speed_m_per_s / characteristic_length_m)


def combine_in_series(conductance_a: float, conductance_b: float) -> float:
    """
    直列につながった2つのコンダクタンスを合成する。

    気孔と境界層は直列（水蒸気は必ず両方を通る）なので、
    抵抗の和 = 1/a + 1/b、その逆数が合成コンダクタンスになる。

        g = 1 / (1/a + 1/b)

    どちらかが 0 なら全体も 0（通れない）。
    """
    if conductance_a <= 0 or conductance_b <= 0:
        return 0.0

    return 1.0 / (1.0 / conductance_a + 1.0 / conductance_b)


def photosynthesis_light_function(par_w_per_m2: float, kl_w_per_m2: float) -> float:
    """
    光合成の光応答関数 f(Stp)。

        f = Stp / (Stp + KL)

    直角双曲線。KL は「最大の半分の光合成をする光量」。
    弱光では光に比例し、強光では 1 に漸近する（光飽和）。

    Returns:
        0〜1 の係数
    """
    denominator = par_w_per_m2 + kl_w_per_m2
    if denominator <= 0:
        return 0.0

    return par_w_per_m2 / denominator


def photosynthesis_light_response_mol(
    par_w_per_m2: float,
    quantum_yield_mol_per_mol: float,
    pmax_mol_per_m2_s: float,
    curvature_theta: float,
    par_umol_per_joule: float,
) -> float:
    """
    非直角双曲線による光応答。光だけで決まる総光合成速度 [mol/(m²·s)] を返す。

        P = [ φI + Pmax − √((φI + Pmax)² − 4θ·φI·Pmax) ] / (2θ)

    【各項の意味】
        φI   : 弱光のときの光合成速度。光子が来た分だけ固定できる領域。
        Pmax : 強光のときの上限。酵素の処理能力で頭打ちになる領域。
        θ    : この2つの領域をどれくらい角ばってつなぐか。
               θ→0 なら直角双曲線（ゆるやかに飽和）、
               θ→1 なら2直線を折れ線でつないだ形（急に飽和）。

    【直角双曲線との違い】
    直角双曲線では初期勾配と飽和点が KL ひとつで結ばれてしまうため、
    片方を合わせるともう片方が外れる。θ を入れると独立に決められる。

    Args:
        par_w_per_m2: 光合成有効放射 [W/m²]
        quantum_yield_mol_per_mol: 量子収率 φ [mol CO2/mol光子]
        pmax_mol_per_m2_s: 光飽和時の総光合成速度 Pmax [mol/(m²·s)]
        curvature_theta: 曲率 θ（0 < θ ≦ 1）
        par_umol_per_joule: PAR を光子数に直す係数 [μmol光子/J]

    Returns:
        光だけで決まる総光合成速度 [mol/(m²·s)]
    """
    if not 0.0 < curvature_theta <= 1.0:
        raise ValueError(
            f"曲率 θ は 0 より大きく 1 以下でなければならない。"
            f"渡された値: {curvature_theta}"
        )
    if par_w_per_m2 <= 0.0 or pmax_mol_per_m2_s <= 0.0:
        return 0.0

    # PAR [W/m²] → 光子束 [mol光子/(m²·s)]
    photon_flux_mol = par_w_per_m2 * par_umol_per_joule * 1e-6
    light_limited = quantum_yield_mol_per_mol * photon_flux_mol

    total = light_limited + pmax_mol_per_m2_s
    discriminant = total ** 2 - 4.0 * curvature_theta * light_limited * pmax_mol_per_m2_s
    # 丸め誤差で判別式がわずかに負になることがあるので 0 で止める
    root = math.sqrt(max(discriminant, 0.0))

    return (total - root) / (2.0 * curvature_theta)


def photosynthesis_temp_function(
    leaf_temp_c: float,
    coef_a_c: float,
    optimal_temp_c: float,
) -> float:
    """
    光合成の温度応答関数 h(Tl)。

        h = (2*(Tl+a)² * (Tm+a)² - (Tl+a)⁴) / (Tm+a)⁴

    最適温度 Tm で 1 になり、そこから離れると下がる釣り鐘型。

    【注意】
    高温側では式が負の値を返す。Tm=28℃, a=5 の場合、
    葉温が約41.7℃を超えると負になる。光合成が負になることは
    ないので 0 で下限を切る（暗呼吸は別途 Rd で差し引く）。

    Returns:
        0〜1 の係数
    """
    leaf_term = (leaf_temp_c + coef_a_c) ** 2
    optimal_term = (optimal_temp_c + coef_a_c) ** 2

    if optimal_term <= 0:
        raise ValueError(
            f"光合成温度関数の分母が 0 以下です: "
            f"最適温度 {optimal_temp_c} ℃, 温度係数 {coef_a_c} ℃"
        )

    value = (2.0 * leaf_term * optimal_term - leaf_term ** 2) / (optimal_term ** 2)

    return max(value, 0.0)


def gross_photosynthesis_mol(
    total_conductance_co2: float,
    air_co2_ppm: float,
    max_photosynthesis_mol: float,
    kc_ppm: float,
) -> float:
    """
    総光合成速度 P を求める。

    【理屈】
    光合成速度は次の2つを同時に満たす必要がある。

      (1) 輸送の式  : P = gc * (大気CO2 - 葉緑体CO2)
          … 大気から葉緑体までCO2が運ばれる速さ
      (2) 生化学の式: P = Pm * g(Ccc)  ただし g(Ccc) = Ccc / (Ccc + KC)
          … 葉緑体でCO2を固定する速さ

    未知数は P と 葉緑体CO2濃度 Ccc の2つ。連立させると Ccc が消去でき、
    P についての2次方程式になる。その小さいほうの解が答え。

        P = gc * [ (Ca + KC + Pm/gc) - sqrt((Ca + KC + Pm/gc)² - 4*Ca*Pm/gc) ] / 2

    （光合成.xlsx では葉緑体CO2濃度を手で調整して合わせていたが、
      ここでは解析解を使うので調整は不要。）

    Args:
        total_conductance_co2: CO2の総コンダクタンス gc [mol/(m²·s)]
        air_co2_ppm: 大気CO2濃度 Cca [ppm]
        max_photosynthesis_mol: CO2最大時光合成速度 Pm [mol/(m²·s)]
        kc_ppm: 光合成CO2係数 KC [ppm]

    Returns:
        総光合成速度 P [mol/(m²·s)]
    """
    if total_conductance_co2 <= 0:
        # 気孔が完全に閉じていればCO2は入らない。
        return 0.0

    # ppm（μmol/mol）を mol/mol に直す。
    ca = air_co2_ppm / 1.0e6
    kc = kc_ppm / 1.0e6

    ratio = max_photosynthesis_mol / total_conductance_co2
    b = ca + kc + ratio

    discriminant = b ** 2 - 4.0 * ca * ratio
    if discriminant < 0:
        # 数学的には常に 0 以上になるが、丸め誤差でわずかに負になりうる。
        # その場合は重解として扱う。
        discriminant = 0.0

    return total_conductance_co2 * (b - math.sqrt(discriminant)) / 2.0


def photorespiration_mmol(
    gross_photosynthesis_mmol_value: float,
    air_co2_ppm: float,
    rp_max_ratio: float,
    co2_compensation_ppm: float,
) -> float:
    """
    光呼吸 Rp を求める。

    【光呼吸とは】
    ルビスコ（CO2を固定する酵素）が、CO2の代わりに酸素と反応してしまう
    無駄な反応。CO2濃度が低いほど起きやすい。
    CO2施用が効くのは、この無駄を減らせるからでもある。

        Rp = Rpmax * (1 - Cca/Ccp) * P

    CO2濃度が Ccp（補償点）に達すると光呼吸はゼロになる。

    Returns:
        光呼吸速度 Rp [mmol/(m²·s)]
    """
    if co2_compensation_ppm <= 0:
        raise ValueError(
            f"光呼吸がなくなるCO2濃度が 0 以下です: {co2_compensation_ppm} ppm"
        )

    ratio = rp_max_ratio * (1.0 - air_co2_ppm / co2_compensation_ppm)

    # CO2が補償点を超えると式は負になるが、光呼吸が負（＝CO2を増やす）
    # ことはないので 0 で下限を切る。
    return max(ratio, 0.0) * gross_photosynthesis_mmol_value


def dark_respiration_mmol(
    leaf_temp_c: float,
    rd20_mmol: float,
    q10: float,
) -> float:
    """
    暗呼吸 Rd を求める。

    【暗呼吸とは】
    植物が生きるために常時行っている呼吸。光の有無に関係なく起きる。
    温度が10℃上がるとおよそ2倍になる（Q10 = 2）。

        Rd = Rd20 * Q10^((Tl - 20) / 10)

    夜温を下げると呼吸による消耗が減るのは、この関係による。

    Returns:
        暗呼吸速度 Rd [mmol/(m²·s)]
    """
    return rd20_mmol * q10 ** ((leaf_temp_c - 20.0) / 10.0)


def calculate_leaf_gas_exchange(
    air_temp_c: float,
    vapor_pressure_kpa: float,
    solar_radiation_w_per_m2: float,
    air_co2_ppm: float,
    wind_speed_m_per_s: float,
    pressure_kpa: float,
    characteristic_length_m: float,
    leaf_temp_c: float | None = None,
    solve_energy_balance: bool = False,
    longwave_view_factor: float = 1.0,
    stomata_scale_factor: float = STOMATA_SCALE_FACTOR,
    light_response_model: str = LIGHT_RESPONSE_MODEL,
) -> LeafGasExchange:
    """
    葉面積1m²あたりの蒸散速度と純光合成速度を求める。

    Args:
        air_temp_c: 気温 Ta [℃]
        vapor_pressure_kpa: 大気の水蒸気圧 ea [kPa]
            （乾湿球から求める場合は psychrometry.vapor_pressure_from_wet_bulb_kpa）
        solar_radiation_w_per_m2: 日射 St [W/m²]
            光合成有効日射 Stp はこの半分として扱う（光合成.xlsx の扱いに合わせる）
        air_co2_ppm: 大気CO2濃度 Cca [ppm]
        wind_speed_m_per_s: 風速 u [m/s]
        pressure_kpa: 大気圧 Pa [kPa]
        characteristic_length_m: 葉の特性長 d [m]
        leaf_temp_c: 葉温 Tl [℃]。指定するとその値を使う。
        solve_energy_balance: True なら葉温をエネルギー収支から求める。
            False（既定）なら葉温＝気温とする（光合成.xlsx と同じ扱い）。

    【葉温の決め方について】
    既定では葉温＝気温としている。これは光合成.xlsx の扱いに合わせたもので、
    xlsx との照合テストを通すために必要。
    ただし実データを流して1日の蒸散量を出す用途では、この仮定により
    蒸散量が物理的にありえないほど小さくなる（日射の18%しか潜熱に回らない）。
    実運用では solve_energy_balance=True を使う。

    Returns:
        LeafGasExchange（蒸散速度・光合成速度と、その中間値）
    """
    if leaf_temp_c is None:
        if solve_energy_balance:
            # エネルギー収支が釣り合う葉温を数値的に求める。
            # 循環importを避けるため、ここで読み込む。
            from core.energy_balance import solve_leaf_temperature

            balance = solve_leaf_temperature(
                air_temp_c=air_temp_c,
                vapor_pressure_kpa=vapor_pressure_kpa,
                solar_radiation_w_per_m2=solar_radiation_w_per_m2,
                air_co2_ppm=air_co2_ppm,
                wind_speed_m_per_s=wind_speed_m_per_s,
                pressure_kpa=pressure_kpa,
                characteristic_length_m=characteristic_length_m,
                longwave_view_factor=longwave_view_factor,
                stomata_scale_factor=stomata_scale_factor,
            )
            leaf_temp_c = balance.leaf_temp_c
        else:
            # 光合成.xlsx と同じく葉温＝気温と仮定する。
            leaf_temp_c = air_temp_c

    params = PHOTOSYNTHESIS_PARAMS

    # 光合成有効日射は全日射の半分とする。
    # （太陽光のエネルギーのうち、光合成に使える波長域が約半分のため）
    par_w_per_m2 = solar_radiation_w_per_m2 / 2.0

    # ---- 水蒸気の輸送 ----------------------------------------------------
    # 葉の内部は常に飽和していると仮定する（気孔腔は水で満たされている）。
    leaf_saturation_kpa = saturation_vapor_pressure_kpa(leaf_temp_c)

    air_vapor_fraction = mole_fraction(vapor_pressure_kpa, pressure_kpa)
    leaf_vapor_fraction = mole_fraction(leaf_saturation_kpa, pressure_kpa)

    gvs = stomatal_conductance_water(par_w_per_m2, air_co2_ppm, stomata_scale_factor)
    gva = boundary_layer_conductance_water(wind_speed_m_per_s, characteristic_length_m)
    gv = combine_in_series(gvs, gva)

    # 蒸散速度（式6.7）: 濃度差 × コンダクタンス
    transpiration_mol = gv * (leaf_vapor_fraction - air_vapor_fraction)
    transpiration_mg = transpiration_mol * MOLAR_MASS_WATER_G_PER_MOL * 1000.0

    # ---- CO2の輸送 -------------------------------------------------------
    # CO2は水蒸気より分子が大きく拡散が遅いため、コンダクタンスを割り引く。
    gcs = gvs * CO2_CONDUCTANCE_RATIO["STOMATA"]
    gca = gva * CO2_CONDUCTANCE_RATIO["BOUNDARY"]
    gc = combine_in_series(gcs, gca)

    # ---- 光合成 ----------------------------------------------------------
    h_temp = photosynthesis_temp_function(
        leaf_temp_c, params["TEMP_COEF_A_C"], params["OPTIMAL_TEMP_C"]
    )

    # 光応答。既定は非直角双曲線、"rectangular" を選ぶと光合成.xlsx と同じ式。
    if light_response_model == "non_rectangular":
        pmax_mol = NON_RECTANGULAR_LIGHT_PARAMS["PMAX_MMOL_PER_M2_S"] / 1000.0
        light_limited_mol = photosynthesis_light_response_mol(
            par_w_per_m2=par_w_per_m2,
            quantum_yield_mol_per_mol=NON_RECTANGULAR_LIGHT_PARAMS[
                "QUANTUM_YIELD_MOL_PER_MOL"
            ],
            pmax_mol_per_m2_s=pmax_mol,
            curvature_theta=NON_RECTANGULAR_LIGHT_PARAMS["CURVATURE_THETA"],
            par_umol_per_joule=PAR_UMOL_PER_JOULE,
        )
        # 報告用の「Pmax の何割に届いたか」。直角双曲線の f(Stp) と同じ意味。
        f_light = light_limited_mol / pmax_mol if pmax_mol > 0.0 else 0.0
    elif light_response_model == "rectangular":
        f_light = photosynthesis_light_function(par_w_per_m2, params["KL_W_PER_M2"])
        light_limited_mol = params["PMAX_MMOL_PER_M2_S"] * f_light / 1000.0
    else:
        raise ValueError(
            f"光応答の式は 'non_rectangular' か 'rectangular' のどちらか。"
            f"渡された値: '{light_response_model}'"
        )

    # CO2が十分にあるときの光合成速度（光と温度だけで決まる上限）
    pm_mol = light_limited_mol * h_temp

    gross_mol = gross_photosynthesis_mol(
        gc, air_co2_ppm, pm_mol, params["KC_PPM"]
    )
    gross_mmol = gross_mol * 1000.0

    rp_mmol = photorespiration_mmol(
        gross_mmol,
        air_co2_ppm,
        params["RP_MAX_RATIO"],
        params["CO2_COMPENSATION_PPM"],
    )
    rd_mmol = dark_respiration_mmol(
        leaf_temp_c, params["RD20_MMOL_PER_M2_S"], params["Q10"]
    )

    # 純光合成 = 総光合成 - 光呼吸 - 暗呼吸
    # 夜間は総光合成が 0 なので、この値は負（呼吸による放出）になる。
    net_mmol = gross_mmol - rp_mmol - rd_mmol
    net_mg = net_mmol * MOLAR_MASS_CO2_G_PER_MOL

    return LeafGasExchange(
        stomatal_conductance_water=gvs,
        boundary_conductance_water=gva,
        total_conductance_water=gv,
        air_vapor_mole_fraction=air_vapor_fraction,
        leaf_vapor_mole_fraction=leaf_vapor_fraction,
        transpiration_mol_per_m2_s=transpiration_mol,
        transpiration_mg_per_m2_s=transpiration_mg,
        stomatal_conductance_co2=gcs,
        boundary_conductance_co2=gca,
        total_conductance_co2=gc,
        light_function=f_light,
        temp_function=h_temp,
        max_photosynthesis_mol_per_m2_s=pm_mol,
        gross_photosynthesis_mmol=gross_mmol,
        photorespiration_mmol=rp_mmol,
        dark_respiration_mmol=rd_mmol,
        net_photosynthesis_mmol=net_mmol,
        net_photosynthesis_mg=net_mg,
        leaf_temp_c=leaf_temp_c,
        air_temp_c=air_temp_c,
    )
