"""
Penman-Monteith 式（教科書 式14.12）の検証。

【確認すること】
1. 見かけの乾湿計定数 γ* が教科書P251の式と一致するか
2. 飽和水蒸気圧曲線の傾き s が妥当か
3. 多層モデル（core/canopy.py）と結果が一致するか  ← 最も重要
4. 物理的に妥当な振る舞いをするか

【3が重要な理由】
多層モデルは群落を20層に分け、層ごとに葉温を数値的に探して解く。
Penman-Monteith 式は群落を1枚の葉とみなし、葉温を式の上で消去してある。
まったく別の道筋なので、両者が一致すれば実装が正しいと確かめられる。

【実行のしかた】
    cd "c:\\Users\\kimij\\OneDrive\\Desktop\\claude作業場\\蒸散光合成モデル"
    python tests\\test_penman_monteith.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (  # noqa: E402
    LEAF_CHARACTERISTIC_LENGTH_M,
    LOCAL_PRESSURE_KPA,
    PSYCHROMETRIC_GAMMA_PER_K,
    WIND_SPEED_M_PER_S,
)
from core.canopy import calculate_canopy_gas_exchange  # noqa: E402
from core.penman_monteith import (  # noqa: E402
    apparent_psychrometric_constant,
    calculate_from_canopy,
    saturation_slope_per_k,
)

# 多層モデルとPM式の差をどこまで許すか
# 日中は1%以内に収まるが、移流項が支配的な夕方はやや広がる
TOLERANCE_DAYTIME = 0.02
TOLERANCE_OVERALL = 0.05

# 実データから取った正午前後（11〜14時）の月平均条件
# 気温[℃], 水蒸気圧[kPa], 日射[W/m²], CO2[ppm]
REAL_CONDITIONS = [
    ("2月  正午前後", 23.7, 1.95, 339.0, 538.0, True),
    ("11月 正午前後", 24.7, 2.00, 278.0, 512.0, True),
    ("5月  正午前後", 26.1, 1.83, 489.0, 422.0, True),
    ("冬の朝 8時",    14.0, 1.30, 150.0, 900.0, False),
    ("夕方",          18.0, 1.60,  80.0, 600.0, False),
]


def test_gamma_star_matches_textbook() -> bool:
    """
    見かけの乾湿計定数 γ* が教科書P251の式と一致するかを確かめる。

    教科書は基準蒸発散（高さ12cmの草本群落）の例で、
    群落コンダクタンス 0.6 mol/(m²·s)、境界層コンダクタンス gHa = 0.2u
    のとき

        γ* = 6.67e-4 (1 + u/3)  [1/℃]

    となることを示している。この式を再現できれば、
    γ* = γ · gHa / gv という実装が正しいと確かめられる。
    """
    print("=" * 78)
    print("見かけの乾湿計定数 γ* の検証（教科書P251の基準蒸発散）")
    print("=" * 78)
    print()
    print("  基準群落: 群落コンダクタンス 0.6 mol/(m²·s)、gHa = gva = 0.2u")
    print("  教科書の式: γ* = 6.67e-4 (1 + u/3)")
    print()
    print(f"{'風速 u':>8}{'gv':>12}{'gHa':>10}{'実装のγ*':>14}{'教科書のγ*':>14}{'判定':>8}")
    print("-" * 78)

    all_passed = True
    for u in [0.5, 1.0, 2.0, 3.0, 5.0]:
        # 教科書の基準群落の設定
        canopy_conductance = 0.6
        boundary_conductance = 0.2 * u
        # 気孔と境界層の直列合成
        gv = (canopy_conductance * boundary_conductance) / (
            canopy_conductance + boundary_conductance
        )
        gha = boundary_conductance

        implemented = apparent_psychrometric_constant(gha, gv)
        textbook = PSYCHROMETRIC_GAMMA_PER_K * (1 + u / 3)

        error = abs(implemented - textbook) / textbook
        passed = error < 1e-12
        all_passed = all_passed and passed

        print(
            f"{u:>8.1f}{gv:>12.5f}{gha:>10.3f}"
            f"{implemented:>14.6f}{textbook:>14.6f}{'OK' if passed else '不一致':>8}"
        )

    print("-" * 78)
    if all_passed:
        print("  γ* = γ · gHa / gv という実装が教科書の式と完全に一致しました。")
    else:
        print("  教科書の式と一致しません。")
    print()

    return all_passed


def test_saturation_slope() -> bool:
    """
    飽和水蒸気圧曲線の傾き s が妥当かを確かめる。

    s は気温が上がるほど大きくなるはず（指数関数の微分なので）。
    数値微分と比べて、解析微分が正しいことも確認する。
    """
    print("=" * 78)
    print("飽和水蒸気圧曲線の傾き s の検証")
    print("=" * 78)

    import math

    from config import TETENS_A_KPA, TETENS_B, TETENS_C_C

    def es(t):
        return TETENS_A_KPA * math.exp(TETENS_B * t / (t + TETENS_C_C))

    print()
    print(f"{'気温℃':>8}{'解析微分 s':>14}{'数値微分':>14}{'相対誤差':>12}{'判定':>8}")
    print("-" * 78)

    all_passed = True
    previous = None
    monotonic = True

    for temp in [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]:
        analytic = saturation_slope_per_k(temp, LOCAL_PRESSURE_KPA)

        # 中心差分による数値微分
        delta = 0.001
        numeric = (es(temp + delta) - es(temp - delta)) / (2 * delta) / LOCAL_PRESSURE_KPA

        error = abs(analytic - numeric) / numeric
        passed = error < 1e-6
        all_passed = all_passed and passed

        if previous is not None and analytic <= previous:
            monotonic = False
        previous = analytic

        print(
            f"{temp:>8.1f}{analytic:>14.6f}{numeric:>14.6f}"
            f"{error:>12.2e}{'OK' if passed else '不一致':>8}"
        )

    print("-" * 78)
    print(f"  気温とともに単調に増加: {'OK' if monotonic else '失敗'}")
    print()

    return all_passed and monotonic


def test_agreement_with_layered_model() -> bool:
    """
    多層モデルと Penman-Monteith 式が一致するかを確かめる。

    これが本テストの中心。
    """
    print("=" * 78)
    print("多層モデルとの一致（最も重要な検証）")
    print("=" * 78)
    print()
    print("  多層モデル      : 群落を20層に分け、層ごとに葉温を数値的に解く")
    print("  Penman-Monteith : 群落を1枚の葉とみなし、葉温を式の上で消去")
    print()
    print(f"{'場面':<16}{'多層':>10}{'PM式':>10}{'比':>8}"
          f"{'放射項':>9}{'判定':>8}")
    print(f"{'':<16}{'mg/m²s':>10}{'mg/m²s':>10}")
    print("-" * 78)

    all_passed = True
    for label, temp, ea, solar, co2, is_daytime in REAL_CONDITIONS:
        canopy = calculate_canopy_gas_exchange(
            air_temp_c=temp,
            vapor_pressure_kpa=ea,
            solar_above_w_per_m2=solar,
            air_co2_ppm=co2,
            wind_speed_m_per_s=WIND_SPEED_M_PER_S,
            pressure_kpa=LOCAL_PRESSURE_KPA,
            characteristic_length_m=LEAF_CHARACTERISTIC_LENGTH_M,
            lai=2.16,
            solve_energy_balance=True,
        )
        pm = calculate_from_canopy(canopy, temp, ea, LOCAL_PRESSURE_KPA)

        layered = canopy.transpiration_mg_per_m2_s
        ratio = pm.transpiration_mg_per_m2_s / layered if layered else float("nan")

        tolerance = TOLERANCE_DAYTIME if is_daytime else TOLERANCE_OVERALL
        passed = abs(ratio - 1.0) < tolerance
        all_passed = all_passed and passed

        print(
            f"{label:<16}{layered:>10.2f}{pm.transpiration_mg_per_m2_s:>10.2f}"
            f"{ratio:>8.3f}{pm.radiation_fraction * 100:>8.0f}%"
            f"{'OK' if passed else '不一致':>8}"
        )

    print("-" * 78)
    print(f"  許容範囲: 日中 ±{TOLERANCE_DAYTIME * 100:.0f}% / "
          f"その他 ±{TOLERANCE_OVERALL * 100:.0f}%")
    print()
    print("  ※ 夕方にずれが大きくなるのは、移流項が支配的になるため。")
    print("     多層モデルは層ごとに葉温が違うが、PM式は1枚の葉として")
    print("     平均化するので、非線形な項の扱いで差が出る。")
    print()

    return all_passed


def test_physical_behavior() -> bool:
    """Penman-Monteith 式が物理的に妥当な振る舞いをするか確かめる。"""
    print("=" * 78)
    print("物理的な妥当性")
    print("=" * 78)

    def run(temp, ea, solar, co2, wind=WIND_SPEED_M_PER_S, lai=2.16):
        c = calculate_canopy_gas_exchange(
            air_temp_c=temp, vapor_pressure_kpa=ea, solar_above_w_per_m2=solar,
            air_co2_ppm=co2, wind_speed_m_per_s=wind, pressure_kpa=LOCAL_PRESSURE_KPA,
            characteristic_length_m=LEAF_CHARACTERISTIC_LENGTH_M,
            lai=lai, solve_energy_balance=True,
        )
        return calculate_from_canopy(c, temp, ea, LOCAL_PRESSURE_KPA)

    checks = []

    # 日射を増やすと蒸散が増える
    low = run(25.0, 2.0, 200.0, 500.0)
    high = run(25.0, 2.0, 500.0, 500.0)
    checks.append((
        "日射が増えると蒸散が増える",
        high.transpiration_mg_per_m2_s > low.transpiration_mg_per_m2_s,
        f"200 W/m²: {low.transpiration_mg_per_m2_s:.1f} → "
        f"500 W/m²: {high.transpiration_mg_per_m2_s:.1f} mg/m²s",
    ))

    # 飽差を大きくすると蒸散が増える
    humid = run(25.0, 2.8, 300.0, 500.0)
    dry = run(25.0, 1.5, 300.0, 500.0)
    checks.append((
        "空気が乾くと蒸散が増える",
        dry.transpiration_mg_per_m2_s > humid.transpiration_mg_per_m2_s,
        f"ea=2.8 kPa: {humid.transpiration_mg_per_m2_s:.1f} → "
        f"ea=1.5 kPa: {dry.transpiration_mg_per_m2_s:.1f} mg/m²s",
    ))

    # 飽和状態では移流項がゼロになる
    import math
    from config import TETENS_A_KPA, TETENS_B, TETENS_C_C
    es_20 = TETENS_A_KPA * math.exp(TETENS_B * 20.0 / (20.0 + TETENS_C_C))
    saturated = run(20.0, es_20, 300.0, 500.0)
    checks.append((
        "湿度100%では移流項がゼロ",
        abs(saturated.advection_term_w_per_m2) < 1e-9,
        f"移流項 = {saturated.advection_term_w_per_m2:.3e} W/m²",
    ))

    # 土壌熱フラックスは正味放射より十分小さい（白マルチ被覆率100%）
    noon = run(25.0, 2.0, 400.0, 500.0)
    checks.append((
        "土壌熱フラックスが正味放射に対して小さい（白マルチ100%）",
        0 <= noon.soil_heat_flux_w_per_m2 < 0.1 * noon.net_radiation_w_per_m2,
        f"G = {noon.soil_heat_flux_w_per_m2:.1f} W/m² / "
        f"正味放射 {noon.net_radiation_w_per_m2:.1f} W/m² "
        f"（{noon.soil_heat_flux_w_per_m2 / noon.net_radiation_w_per_m2 * 100:.1f}%）",
    ))

    # 気温が高いほど傾き s が大きく、放射項の効率が上がる
    cool = run(15.0, 1.2, 300.0, 500.0)
    warm = run(30.0, 2.5, 300.0, 500.0)
    checks.append((
        "気温が高いほど傾き s が大きい",
        warm.slope_s > cool.slope_s,
        f"15℃: s={cool.slope_s:.5f} → 30℃: s={warm.slope_s:.5f} /K",
    ))

    print()
    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '失敗'}  {label}")
        print(f"        {detail}")
    print()

    return all_passed


def main() -> int:
    results = [
        test_gamma_star_matches_textbook(),
        test_saturation_slope(),
        test_agreement_with_layered_model(),
        test_physical_behavior(),
    ]

    print("=" * 78)
    if all(results):
        print("第5段階（Penman-Monteith 式）は完了です。")
        print("独立した2つの定式化が一致したので、蒸散モデルの信頼性が")
        print("さらに一段上がりました。")
        print("=" * 78)
        return 0

    print("失敗した項目があります。")
    print("=" * 78)
    return 1


if __name__ == "__main__":
    sys.exit(main())
