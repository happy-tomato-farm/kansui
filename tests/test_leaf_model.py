"""
葉モデルが 光合成.xlsx を正しく再現できているかを検証する。

【なぜこのテストが必要か】
このあと群落へのスケーリング・日積算・較正と計算を積み上げていく。
土台である葉モデルの移植にミスがあると、最終結果がずれたときに
「モデルが悪いのか、移植をミスしたのか」を切り分けられなくなる。

光合成.xlsx には計算例が1つ入っている。その入力と出力を使えば、
移植の正しさを数値で確認できる。ここが通ってから先に進む。

【実行のしかた】
プロジェクトのフォルダで、PowerShell から次を実行する。

    cd "c:\\Users\\kimij\\OneDrive\\Desktop\\claude作業場\\蒸散光合成モデル"
    python tests\\test_leaf_model.py

すべて一致すれば「すべて一致しました」と表示される。
"""

import sys
from pathlib import Path

# このファイルは tests/ の中にあるので、親フォルダ（プロジェクトのルート）を
# Python の検索パスに加える。こうしないと config や core を読み込めない。
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.leaf_model import calculate_leaf_gas_exchange  # noqa: E402
from core.psychrometry import (  # noqa: E402
    saturation_vapor_pressure_kpa,
    vapor_pressure_from_rh_kpa,
)


# =============================================================================
# 光合成.xlsx の計算例（10行目）
# =============================================================================
# 入力条件
XLSX_INPUT = {
    "air_temp_c": 17.0,             # 気温 Ta
    "relative_humidity": 0.6,       # 相対湿度 hr
    "solar_radiation_w_per_m2": 200.0,  # 日射 St（光合成有効日射 Stp はこの半分）
    "air_co2_ppm": 600.0,           # 大気CO2濃度 Cca
    "wind_speed_m_per_s": 0.2,      # 風速 u
    "characteristic_length_m": 0.0504,  # 特性長 d
}

# 光合成.xlsx が使っている大気圧。
# ※標高補正した 99.76 kPa ではなく、海面標準の 101.3 kPa をそのまま使っている。
#   xlsx の再現が目的なのでここでは 101.3 を使う。
#   実運用では config.LOCAL_PRESSURE_KPA（99.76 kPa）を使う。
XLSX_PRESSURE_KPA = 101.3

# 期待される出力（xlsx のセルの値）
XLSX_EXPECTED = {
    "es(Ta) 飽和水蒸気圧 [kPa]":        ("es_ta",   1.9361633314854947),
    "ea 大気水蒸気圧 [kPa]":            ("ea",      1.1616979988912968),
    "Cva 大気水蒸気モル分率 [mol/mol]": ("air_vapor_mole_fraction",    0.011467897323704805),
    "Cvs 葉内水蒸気モル分率 [mol/mol]": ("leaf_vapor_mole_fraction",   0.019113162206174678),
    "gvs 気孔C(水蒸気) [mol/m2s]":      ("stomatal_conductance_water", 0.08349999999999999),
    "gva 境界層C(水蒸気) [mol/m2s]":    ("boundary_conductance_water", 0.2928310092869264),
    "gv 総C(水蒸気) [mol/m2s]":         ("total_conductance_water",    0.06497309196440906),
    "Fv 蒸散速度 [mol/m2s]":            ("transpiration_mol_per_m2_s", 0.0004967364983009821),
    "Fv 蒸散速度 [mg/m2s]":             ("transpiration_mg_per_m2_s",  8.941256969417678),
    "gcs 気孔C(CO2) [mol/m2s]":         ("stomatal_conductance_co2",   0.05511),
    "gca 境界層C(CO2) [mol/m2s]":       ("boundary_conductance_co2",   0.2196232569651948),
    "gc 総C(CO2) [mol/m2s]":            ("total_conductance_co2",      0.04405523315615639),
    "f(Stp) 光関数 [0-1]":              ("light_function",             0.4),
    "h(Tl) 温度関数 [0-1]":             ("temp_function",              0.691358024691358),
    "Pm CO2最大時光合成 [mol/m2s]":     ("max_photosynthesis_mol_per_m2_s", 1.2568888888888887e-05),
    "P 総光合成速度 [mmol/m2s]":        ("gross_photosynthesis_mmol",  0.007885979650700354),
    "Rp 光呼吸 [mmol/m2s]":             ("photorespiration_mmol",      0.0011828969476050531),
    "Rd 暗呼吸 [mmol/m2s]":             ("dark_respiration_mmol",      0.0012021335466072284),
    "Pn 純光合成速度 [mmol/m2s]":       ("net_photosynthesis_mmol",    0.005500949156488073),
    "Pn 純光合成速度 [mg/m2s]":         ("net_photosynthesis_mg",      0.24204176288547521),
}

# 許容する相対誤差。浮動小数点の丸め誤差だけを許す厳しさにしてある。
TOLERANCE = 1e-9


def run_xlsx_comparison() -> bool:
    """
    光合成.xlsx の計算例と Python 版の結果を突き合わせる。

    Returns:
        すべて一致したら True
    """
    # xlsx は相対湿度から水蒸気圧を求めている（式3.11）
    ea = vapor_pressure_from_rh_kpa(
        XLSX_INPUT["air_temp_c"], XLSX_INPUT["relative_humidity"]
    )
    es_ta = saturation_vapor_pressure_kpa(XLSX_INPUT["air_temp_c"])

    result = calculate_leaf_gas_exchange(
        air_temp_c=XLSX_INPUT["air_temp_c"],
        vapor_pressure_kpa=ea,
        solar_radiation_w_per_m2=XLSX_INPUT["solar_radiation_w_per_m2"],
        air_co2_ppm=XLSX_INPUT["air_co2_ppm"],
        wind_speed_m_per_s=XLSX_INPUT["wind_speed_m_per_s"],
        pressure_kpa=XLSX_PRESSURE_KPA,
        characteristic_length_m=XLSX_INPUT["characteristic_length_m"],
        # xlsx の再現なので、気孔コンダクタンスの水準補正は掛けない。
        # 実運用では config.STOMATA_SCALE_FACTOR（文献に合わせた倍率）を使う。
        stomata_scale_factor=1.0,
        # 光応答も xlsx と同じ直角双曲線に戻す。
        # 実運用の既定は非直角双曲線（config.LIGHT_RESPONSE_MODEL）。
        light_response_model="rectangular",
    )

    # 中間値のうち、モデルの外で計算したものを取り込む
    computed = {"es_ta": es_ta, "ea": ea}
    for field in result.__dataclass_fields__:
        computed[field] = getattr(result, field)

    print("=" * 78)
    print("光合成.xlsx との照合")
    print("=" * 78)
    print(
        f"入力: 気温 {XLSX_INPUT['air_temp_c']} ℃ / "
        f"相対湿度 {XLSX_INPUT['relative_humidity'] * 100:.0f}% / "
        f"日射 {XLSX_INPUT['solar_radiation_w_per_m2']:.0f} W/m² / "
        f"CO2 {XLSX_INPUT['air_co2_ppm']:.0f} ppm / "
        f"風速 {XLSX_INPUT['wind_speed_m_per_s']} m/s"
    )
    print("-" * 78)
    print(f"{'項目':<34}{'xlsx':>14}{'Python':>14}{'相対誤差':>12}  判定")
    print("-" * 78)

    all_passed = True
    for label, (key, expected) in XLSX_EXPECTED.items():
        actual = computed[key]

        # 相対誤差で比較する。期待値が 0 のときだけ絶対誤差を使う。
        if expected == 0.0:
            relative_error = abs(actual)
        else:
            relative_error = abs(actual - expected) / abs(expected)

        passed = relative_error < TOLERANCE
        all_passed = all_passed and passed

        print(
            f"{label:<34}{expected:>14.6g}{actual:>14.6g}"
            f"{relative_error:>12.2e}  {'OK' if passed else '不一致'}"
        )

    print("-" * 78)
    if all_passed:
        print("すべて一致しました。葉モデルの移植は正しく行われています。")
    else:
        print("一致しない項目があります。上の「不一致」の行を確認してください。")
    print()

    return all_passed


def run_sanity_checks() -> bool:
    """
    xlsx にはない条件で、モデルが物理的におかしな値を返さないか確かめる。

    xlsx の計算例は1点しかないため、その1点だけが合っていても
    式の一部を写し間違えている可能性が残る。極端な条件を入れて
    「ありえない答え」が出ないことを確認しておく。
    """
    print("=" * 78)
    print("物理的な妥当性のチェック")
    print("=" * 78)

    checks = []

    # --- 夜間（日射ゼロ）------------------------------------------------
    # 光合成は 0 になり、呼吸のぶん純光合成は負になるはず。
    night = calculate_leaf_gas_exchange(
        air_temp_c=15.0,
        vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=0.0,
        air_co2_ppm=900.0,
        wind_speed_m_per_s=0.2,
        pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    checks.append((
        "夜間は総光合成が 0",
        night.gross_photosynthesis_mmol == 0.0,
        f"P = {night.gross_photosynthesis_mmol:.6f} mmol/m²s",
    ))
    checks.append((
        "夜間の純光合成は負（呼吸による放出）",
        night.net_photosynthesis_mmol < 0.0,
        f"Pn = {night.net_photosynthesis_mmol:.6f} mmol/m²s",
    ))

    # --- 飽和状態（湿度100%）--------------------------------------------
    # 葉温＝気温なら葉内外の水蒸気モル分率が等しくなり、蒸散は 0 になるはず。
    saturated_ea = saturation_vapor_pressure_kpa(20.0)
    saturated = calculate_leaf_gas_exchange(
        air_temp_c=20.0,
        vapor_pressure_kpa=saturated_ea,
        solar_radiation_w_per_m2=300.0,
        air_co2_ppm=400.0,
        wind_speed_m_per_s=0.2,
        pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    checks.append((
        "湿度100%では蒸散が 0",
        abs(saturated.transpiration_mol_per_m2_s) < 1e-12,
        f"Fv = {saturated.transpiration_mol_per_m2_s:.3e} mol/m²s",
    ))

    # --- 高温（45℃）------------------------------------------------------
    # 温度関数が負にならず 0 で止まること。
    hot = calculate_leaf_gas_exchange(
        air_temp_c=45.0,
        vapor_pressure_kpa=3.0,
        solar_radiation_w_per_m2=800.0,
        air_co2_ppm=400.0,
        wind_speed_m_per_s=0.2,
        pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    checks.append((
        "45℃でも温度関数が負にならない",
        hot.temp_function >= 0.0,
        f"h(Tl) = {hot.temp_function:.6f}",
    ))

    # --- 日射を増やすと光合成が増える ------------------------------------
    low_light = calculate_leaf_gas_exchange(
        air_temp_c=25.0, vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=100.0, air_co2_ppm=400.0,
        wind_speed_m_per_s=0.2, pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    high_light = calculate_leaf_gas_exchange(
        air_temp_c=25.0, vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=600.0, air_co2_ppm=400.0,
        wind_speed_m_per_s=0.2, pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    checks.append((
        "日射が増えると純光合成が増える",
        high_light.net_photosynthesis_mmol > low_light.net_photosynthesis_mmol,
        f"100 W/m²: {low_light.net_photosynthesis_mmol:.5f} → "
        f"600 W/m²: {high_light.net_photosynthesis_mmol:.5f} mmol/m²s",
    ))

    # --- CO2を増やすと光合成が増える -------------------------------------
    low_co2 = calculate_leaf_gas_exchange(
        air_temp_c=25.0, vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=400.0, air_co2_ppm=400.0,
        wind_speed_m_per_s=0.2, pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    high_co2 = calculate_leaf_gas_exchange(
        air_temp_c=25.0, vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=400.0, air_co2_ppm=800.0,
        wind_speed_m_per_s=0.2, pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    checks.append((
        "CO2が増えると純光合成が増える",
        high_co2.net_photosynthesis_mmol > low_co2.net_photosynthesis_mmol,
        f"400 ppm: {low_co2.net_photosynthesis_mmol:.5f} → "
        f"800 ppm: {high_co2.net_photosynthesis_mmol:.5f} mmol/m²s",
    ))

    # --- 風速を上げると蒸散が増える --------------------------------------
    calm = calculate_leaf_gas_exchange(
        air_temp_c=25.0, vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=400.0, air_co2_ppm=400.0,
        wind_speed_m_per_s=0.2, pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    windy = calculate_leaf_gas_exchange(
        air_temp_c=25.0, vapor_pressure_kpa=1.5,
        solar_radiation_w_per_m2=400.0, air_co2_ppm=400.0,
        wind_speed_m_per_s=0.5, pressure_kpa=99.76,
        characteristic_length_m=0.0504,
        light_response_model="rectangular",
    )
    checks.append((
        "風速が上がると蒸散が増える",
        windy.transpiration_mg_per_m2_s > calm.transpiration_mg_per_m2_s,
        f"0.2 m/s: {calm.transpiration_mg_per_m2_s:.3f} → "
        f"0.5 m/s: {windy.transpiration_mg_per_m2_s:.3f} mg/m²s "
        f"（{windy.transpiration_mg_per_m2_s / calm.transpiration_mg_per_m2_s:.2f}倍）",
    ))

    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '失敗'}  {label}")
        print(f"        {detail}")

    print("-" * 78)
    if all_passed:
        print("物理的な妥当性のチェックもすべて通りました。")
    else:
        print("妥当性チェックに失敗があります。")
    print()

    return all_passed


def main() -> int:
    xlsx_ok = run_xlsx_comparison()
    sanity_ok = run_sanity_checks()

    if xlsx_ok and sanity_ok:
        print("=" * 78)
        print("第1段階（葉モデルの移植）は完了です。")
        print("=" * 78)
        return 0

    print("=" * 78)
    print("失敗した項目があります。先に進む前に修正が必要です。")
    print("=" * 78)
    return 1


if __name__ == "__main__":
    sys.exit(main())
