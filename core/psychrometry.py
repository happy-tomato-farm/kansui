"""
湿り空気の計算（psychrometry = 湿度計測学）。

乾球温度と湿球温度から、水蒸気圧・相対湿度・飽差（VPD）を求める。

【なぜ湿球温度から計算するのか】
このハウスのセンサーが出力する「湿度」「絶対湿度」の値は壊れており、
2026年4月以降は 0 が並び、それ以前も日中は 0 になる。
一方で「温度3」が湿球温度を測っているため、乾球（温度1）との組み合わせで
湿度を復元できる。復元値は夜明け前 89%・日中 55% と生理的に妥当であり、
乾湿球の逆転（湿球 > 乾球）も 0% であることを確認済み。

【出典】
生物環境物理学の基礎 式3.5, 3.7, 3.8, 3.12, 3.16
および 2027作業記録.xlsx「湿度計算」シート
"""

import math

from config import (
    TETENS_A_KPA,
    TETENS_B,
    TETENS_C_C,
    PSYCHROMETRIC_GAMMA_PER_K,
    STANDARD_PRESSURE_KPA,
)


def saturation_vapor_pressure_kpa(temp_c: float) -> float:
    """
    気温から飽和水蒸気圧を求める（Tetens式、式3.8）。

    飽和水蒸気圧とは「その温度の空気が抱えられる水蒸気の上限」。
    気温が上がるとこの上限は急激に大きくなる（指数関数）。
    同じ相対湿度でも気温が違えば飽差が変わるのは、この非線形性のため。

    Args:
        temp_c: 気温 [℃]

    Returns:
        飽和水蒸気圧 [kPa]
    """
    # 分母が 0 になるのは -240.97℃ のときで、農業の温度域ではありえない。
    # それでも黙って inf を返すと後段で原因不明の異常値になるため、明示的に弾く。
    denominator = temp_c + TETENS_C_C
    if denominator <= 0:
        raise ValueError(
            f"飽和水蒸気圧を計算できません: 気温 {temp_c} ℃ は "
            f"Tetens式の適用範囲外です（-{TETENS_C_C} ℃ 以下）"
        )

    return TETENS_A_KPA * math.exp(TETENS_B * temp_c / denominator)


def air_pressure_kpa(altitude_m: float) -> float:
    """
    標高から現地の大気圧を求める（式3.7）。

    標高が上がると大気圧は下がる。大気圧は乾湿計の式と
    モル分率の計算の両方に効くため、海面標準の 101.3 kPa を
    そのまま使うと系統的な誤差になる。

    Args:
        altitude_m: 標高 [m]

    Returns:
        大気圧 [kPa]
    """
    return STANDARD_PRESSURE_KPA * math.exp(-altitude_m / 8200.0)


def vapor_pressure_from_wet_bulb_kpa(
    dry_bulb_c: float,
    wet_bulb_c: float,
    pressure_kpa: float,
) -> float:
    """
    乾球温度と湿球温度から実際の水蒸気圧を求める（乾湿計の式、式3.16）。

    【理屈】
    湿球温度計は、球部を濡らしたガーゼで包んだ温度計。
    水が蒸発するときに気化熱を奪うので、乾球より低い温度を示す。
    空気が乾いているほど蒸発が盛んになり、温度差が大きくなる。
    この温度差から、逆に空気中の水蒸気量が分かる。

        ea = es(Tw) - γ * Pa * (Ta - Tw)
             ↑          ↑
             湿球温度での  乾湿球の温度差に比例して差し引く
             飽和水蒸気圧  （温度差が大きい＝乾いている＝差し引き大）

    Args:
        dry_bulb_c: 乾球温度（普通の気温）[℃]
        wet_bulb_c: 湿球温度 [℃]
        pressure_kpa: 大気圧 [kPa]

    Returns:
        実際の水蒸気圧 ea [kPa]
    """
    if wet_bulb_c > dry_bulb_c:
        # 物理的には湿球が乾球を上回ることはない。
        # 上回ったらセンサーの異常か、データの列を取り違えている。
        # 黙って計算を続けると水蒸気圧が過大になり、飽差が負になって
        # 「蒸散が逆流する」という異常な結果になるため、ここで止める。
        raise ValueError(
            f"湿球温度が乾球温度を上回っています: "
            f"乾球 {dry_bulb_c} ℃ < 湿球 {wet_bulb_c} ℃。"
            f"センサー異常か、列の取り違えの可能性があります"
        )

    es_wet = saturation_vapor_pressure_kpa(wet_bulb_c)
    depression = dry_bulb_c - wet_bulb_c  # 乾湿球差 [K]

    vapor_pressure = es_wet - PSYCHROMETRIC_GAMMA_PER_K * pressure_kpa * depression

    # 乾湿球差が極端に大きいと、この式は負の値を返す。
    # 水蒸気圧が負になることは物理的にありえないので 0 で下限を切る。
    # （例: 乾球5℃・湿球-5℃ では ea = -0.243 kPa となり、そのまま使うと
    #   飽差が実際より大きく出て蒸散を過大評価してしまう）
    # 上限も飽和水蒸気圧で切る。過飽和は結露として扱われるべきで、
    # 気相の水蒸気圧が飽和を超えて持続することはない。
    es_dry = saturation_vapor_pressure_kpa(dry_bulb_c)

    return min(max(vapor_pressure, 0.0), es_dry)


def vapor_pressure_from_rh_kpa(temp_c: float, relative_humidity: float) -> float:
    """
    気温と相対湿度から水蒸気圧を求める（式3.11）。

    湿球温度がない場合（例: 東ハウス）や、相対湿度を直接入力する
    シミュレーション用。

    Args:
        temp_c: 気温 [℃]
        relative_humidity: 相対湿度（0〜1 の小数。60% なら 0.6）

    Returns:
        水蒸気圧 ea [kPa]
    """
    if not 0.0 <= relative_humidity <= 1.0:
        raise ValueError(
            f"相対湿度は 0〜1 の小数で指定してください（60% なら 0.6）: "
            f"受け取った値 {relative_humidity}"
        )

    return saturation_vapor_pressure_kpa(temp_c) * relative_humidity


def mole_fraction(vapor_pressure_kpa_value: float, pressure_kpa: float) -> float:
    """
    水蒸気圧をモル分率に変換する（式3.5）。

    モル分率とは「空気1molのうち何molが水蒸気か」という割合。
    蒸散の計算では、圧力ではなくモル分率で扱うと式が簡単になる
    （コンダクタンスが mol/(m²·s) の単位で統一できるため）。

        Cv = ea / Pa

    Args:
        vapor_pressure_kpa_value: 水蒸気圧 [kPa]
        pressure_kpa: 大気圧 [kPa]

    Returns:
        水蒸気のモル分率 [mol/mol]
    """
    if pressure_kpa <= 0:
        raise ValueError(f"大気圧が 0 以下です: {pressure_kpa} kPa")

    return vapor_pressure_kpa_value / pressure_kpa


def relative_humidity_from_wet_bulb(
    dry_bulb_c: float,
    wet_bulb_c: float,
    pressure_kpa: float,
) -> float:
    """
    乾球温度と湿球温度から相対湿度を求める。

    Args:
        dry_bulb_c: 乾球温度 [℃]
        wet_bulb_c: 湿球温度 [℃]
        pressure_kpa: 大気圧 [kPa]

    Returns:
        相対湿度（0〜1 の小数）。飽和付近の計算誤差で 1 をわずかに
        超えることがあるため、1.0 で頭打ちにする。
    """
    ea = vapor_pressure_from_wet_bulb_kpa(dry_bulb_c, wet_bulb_c, pressure_kpa)
    es = saturation_vapor_pressure_kpa(dry_bulb_c)

    return min(ea / es, 1.0)


def vapor_pressure_deficit_kpa(
    dry_bulb_c: float,
    wet_bulb_c: float,
    pressure_kpa: float,
) -> float:
    """
    乾球温度と湿球温度から飽差（VPD）を求める（式3.12）。

    【VPD（飽差）とは】
    Vapor Pressure Deficit。空気があとどれだけ水蒸気を受け取れるかの余力。
    飽和水蒸気圧から実際の水蒸気圧を引いた差。

        VPD = es(Ta) - ea

    これが蒸散の駆動力になる。VPD が大きい＝空気が乾いている＝
    葉から水が出ていきやすい。ただし大きすぎると気孔が閉じるため、
    トマトでは日中 0.5〜1.5 kPa 程度が適正とされる。

    Args:
        dry_bulb_c: 乾球温度 [℃]
        wet_bulb_c: 湿球温度 [℃]
        pressure_kpa: 大気圧 [kPa]

    Returns:
        飽差 VPD [kPa]
    """
    ea = vapor_pressure_from_wet_bulb_kpa(dry_bulb_c, wet_bulb_c, pressure_kpa)
    es = saturation_vapor_pressure_kpa(dry_bulb_c)

    # 飽和状態では計算誤差で ea がわずかに es を上回ることがある。
    # 負の飽差は物理的に意味がない（結露状態）ので 0 で下限を切る。
    return max(es - ea, 0.0)


def vpd_to_g_per_m3(vpd_kpa: float, temp_c: float, pressure_kpa: float) -> float:
    """
    飽差の単位を kPa から g/m³ に変換する。

    【なぜ2つの単位があるのか】
    kPa は物理計算に使いやすい単位。
    g/m³ は「1立方メートルの空気があと何グラムの水を受け取れるか」で、
    現場での感覚に近い。トマトでは 3〜7 g/m³ が目安とされる。

    Args:
        vpd_kpa: 飽差 [kPa]
        temp_c: 気温 [℃]
        pressure_kpa: 大気圧 [kPa]

    Returns:
        飽差 [g/m³]
    """
    from config import GAS_CONSTANT_J_PER_MOL_K, MOLAR_MASS_WATER_G_PER_MOL

    # 理想気体の状態方程式から、空気1m³あたりのモル数を求める。
    #   n/V = P / (R * T)      P は Pa、T は K
    temp_k = temp_c + 273.15
    moles_per_m3 = (pressure_kpa * 1000.0) / (GAS_CONSTANT_J_PER_MOL_K * temp_k)

    # 飽差をモル分率にしてから、モル数と分子量を掛けて g/m³ にする。
    vpd_mole_fraction = vpd_kpa / pressure_kpa

    return vpd_mole_fraction * moles_per_m3 * MOLAR_MASS_WATER_G_PER_MOL
