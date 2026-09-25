"""
群落スケーリングの検証。

【確認すること】
1. 層の数を増やしたときに結果が収束するか（数値計算として妥当か）
2. Beer則が正しく効いているか（下層ほど暗い）
3. LAIを増やしたときの応答が生理的に妥当か（頭打ちになる）
4. 実際の環境条件を入れたときに、値がありうる範囲に収まるか

【実行のしかた】
    cd "c:\\Users\\kimij\\OneDrive\\Desktop\\claude作業場\\蒸散光合成モデル"
    python tests\\test_canopy.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (  # noqa: E402
    EXTINCTION_COEFFICIENT_K,
    LEAF_AREA_PER_LEAF_M2,
    LEAF_CHARACTERISTIC_LENGTH_M,
    LOCAL_PRESSURE_KPA,
    WIND_SPEED_M_PER_S,
)
from core.canopy import (  # noqa: E402
    calculate_canopy_gas_exchange,
    leaf_area_index,
    photosynthesis_to_sugar_g_per_m2_s,
    solar_at_depth,
    transpiration_to_mm_per_s,
)

# 検証に使う標準的な日中の条件（2025年11月・中央ハウスの正午の実測値に近い値）
NOON_CONDITION = {
    "air_temp_c": 24.7,
    "vapor_pressure_kpa": 2.0,   # VPD 約1.1 kPa に相当
    "solar_above_w_per_m2": 288.0,
    "air_co2_ppm": 511.0,
    "wind_speed_m_per_s": WIND_SPEED_M_PER_S,
    "pressure_kpa": LOCAL_PRESSURE_KPA,
    "characteristic_length_m": LEAF_CHARACTERISTIC_LENGTH_M,
}


def test_layer_convergence() -> bool:
    """
    層の数を増やしたときに結果が収束するかを確かめる。

    群落を何層に分けるかは計算上の都合であって、物理的な意味はない。
    層を細かくしていけば真の値に近づくはずで、ある程度の層数から先は
    ほとんど変わらなくなる。そうならなければ実装がおかしい。
    """
    print("=" * 78)
    print("層の数による収束の確認")
    print("=" * 78)

    lai = 2.2
    results = {}

    print(f"\n条件: LAI={lai}, 日射={NOON_CONDITION['solar_above_w_per_m2']:.0f} W/m², "
          f"気温={NOON_CONDITION['air_temp_c']}℃, CO2={NOON_CONDITION['air_co2_ppm']:.0f} ppm")
    print("-" * 78)
    print(f"{'層数':>6}{'蒸散 mg/m²s':>16}{'純光合成 mmol/m²s':>20}{'前回との差':>14}")
    print("-" * 78)

    previous_net = None
    for n in [1, 2, 5, 10, 20, 50, 100, 200]:
        result = calculate_canopy_gas_exchange(**NOON_CONDITION, lai=lai, n_layers=n)
        results[n] = result

        if previous_net is None:
            diff_text = "—"
        else:
            diff = abs(result.net_photosynthesis_mmol_per_m2_s - previous_net)
            relative = diff / abs(previous_net) * 100 if previous_net else 0
            diff_text = f"{relative:.3f}%"

        print(
            f"{n:>6}{result.transpiration_mg_per_m2_s:>16.4f}"
            f"{result.net_photosynthesis_mmol_per_m2_s:>20.6f}{diff_text:>14}"
        )
        previous_net = result.net_photosynthesis_mmol_per_m2_s

    # 20層と200層の差が0.5%以内なら、20層で実用上じゅうぶん
    diff_20_vs_200 = abs(
        results[20].net_photosynthesis_mmol_per_m2_s
        - results[200].net_photosynthesis_mmol_per_m2_s
    ) / abs(results[200].net_photosynthesis_mmol_per_m2_s)

    passed = diff_20_vs_200 < 0.005
    print("-" * 78)
    print(
        f"  20層と200層の差: {diff_20_vs_200 * 100:.3f}% "
        f"→ {'OK（20層で十分）' if passed else '収束していない'}"
    )
    print()

    return passed


def test_beer_law() -> bool:
    """Beer則が正しく効いているかを確かめる。"""
    print("=" * 78)
    print("群落内の光の減衰（Beer則）")
    print("=" * 78)

    solar_top = 288.0
    lai = 2.2
    k = EXTINCTION_COEFFICIENT_K

    print(f"\n群落上部の日射 {solar_top:.0f} W/m², LAI={lai}, 消散係数 k={k}")
    print("-" * 78)
    print(f"{'累積LAI':>10}{'日射 W/m²':>14}{'上部に対する割合':>18}")
    print("-" * 78)

    for cumulative in [0.0, 0.5, 1.0, 1.5, 2.0, 2.2]:
        value = solar_at_depth(solar_top, cumulative, k)
        print(f"{cumulative:>10.1f}{value:>14.1f}{value / solar_top * 100:>17.1f}%")

    print("-" * 78)

    checks = []

    # 深いほど暗くなる
    values = [solar_at_depth(solar_top, x, k) for x in [0.0, 1.0, 2.0]]
    checks.append((
        "深いほど日射が減る",
        values[0] > values[1] > values[2],
        f"{values[0]:.1f} > {values[1]:.1f} > {values[2]:.1f} W/m²",
    ))

    # 累積LAI=0 では減衰なし
    checks.append((
        "群落最上部では減衰しない",
        abs(solar_at_depth(solar_top, 0.0, k) - solar_top) < 1e-9,
        f"{solar_at_depth(solar_top, 0.0, k):.1f} W/m²",
    ))

    # 群落の中で実際に何割の光が使われているか
    result = calculate_canopy_gas_exchange(**NOON_CONDITION, lai=lai)
    checks.append((
        "最下層は最上層より暗い",
        result.bottom_layer_solar_w_per_m2 < result.top_layer_solar_w_per_m2,
        f"最上層 {result.top_layer_solar_w_per_m2:.1f} → "
        f"最下層 {result.bottom_layer_solar_w_per_m2:.1f} W/m² "
        f"（{result.bottom_layer_solar_w_per_m2 / result.top_layer_solar_w_per_m2 * 100:.0f}%）",
    ))

    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '失敗'}  {label}")
        print(f"        {detail}")
    print()

    return all_passed


def test_lai_response() -> bool:
    """
    LAIを変えたときの応答が生理的に妥当かを確かめる。

    LAIを増やせば蒸散も光合成も増えるが、下層の葉は暗いので
    増え方は頭打ちになる（収穫逓減）。さらに増やすと、暗い葉の
    呼吸だけがかさんで純光合成が下がることもありうる。
    """
    print("=" * 78)
    print("LAIに対する応答")
    print("=" * 78)

    print(f"\n条件: 日射={NOON_CONDITION['solar_above_w_per_m2']:.0f} W/m²（日中）")
    print("-" * 78)
    print(f"{'LAI':>6}{'蒸散 mg/m²s':>14}{'蒸散 L/m²/h':>14}"
          f"{'純光合成 mmol/m²s':>20}{'LAI1あたり':>14}")
    print("-" * 78)

    results = {}
    for lai in [0.0, 0.5, 1.0, 1.5, 2.0, 2.2, 3.0, 4.0, 5.0]:
        r = calculate_canopy_gas_exchange(**NOON_CONDITION, lai=lai)
        results[lai] = r

        per_lai = (
            r.net_photosynthesis_mmol_per_m2_s / lai if lai > 0 else 0.0
        )
        # mg/(m²·s) → L/(m²·h): mg を g にして 1000倍/時間、水1g=1mL
        liters_per_hour = r.transpiration_mg_per_m2_s * 3600 / 1e6

        print(
            f"{lai:>6.1f}{r.transpiration_mg_per_m2_s:>14.3f}"
            f"{liters_per_hour:>14.4f}"
            f"{r.net_photosynthesis_mmol_per_m2_s:>20.6f}{per_lai:>14.6f}"
        )

    print("-" * 78)

    checks = []

    checks.append((
        "LAI=0 では蒸散も光合成も 0",
        results[0.0].transpiration_mg_per_m2_s == 0.0
        and results[0.0].net_photosynthesis_mmol_per_m2_s == 0.0,
        "両方とも 0",
    ))

    checks.append((
        "LAIが増えると蒸散が増える",
        results[1.0].transpiration_mg_per_m2_s
        < results[2.2].transpiration_mg_per_m2_s
        < results[4.0].transpiration_mg_per_m2_s,
        f"LAI1.0: {results[1.0].transpiration_mg_per_m2_s:.2f} → "
        f"2.2: {results[2.2].transpiration_mg_per_m2_s:.2f} → "
        f"4.0: {results[4.0].transpiration_mg_per_m2_s:.2f} mg/m²s",
    ))

    # 収穫逓減：LAI1あたりの光合成が、LAIが大きいほど小さくなる
    per_lai_1 = results[1.0].net_photosynthesis_mmol_per_m2_s / 1.0
    per_lai_4 = results[4.0].net_photosynthesis_mmol_per_m2_s / 4.0
    checks.append((
        "LAIが大きいほど葉1単位あたりの光合成は下がる（収穫逓減）",
        per_lai_4 < per_lai_1,
        f"LAI=1: {per_lai_1:.6f} → LAI=4: {per_lai_4:.6f} mmol/m²s",
    ))

    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '失敗'}  {label}")
        print(f"        {detail}")
    print()

    return all_passed


def test_realistic_magnitude() -> bool:
    """
    実際の環境条件を入れたとき、値がありうる大きさに収まるかを確かめる。

    ここが最も大事な確認。式が正しくても桁が合っていなければ
    どこかで単位を取り違えている。
    """
    print("=" * 78)
    print("実際の条件での値の大きさ")
    print("=" * 78)

    lai = leaf_area_index(18.0)  # 目標葉枚数 18枚/m²
    print(f"\nLAI = 葉枚数18枚/m² × 葉1枚{LEAF_AREA_PER_LEAF_M2} m² = {lai:.2f}")

    scenarios = [
        ("冬の朝  (2月 8時)",  14.0, 1.0, 150.0, 900.0),
        ("冬の日中(2月12時)",  22.0, 1.5, 350.0, 500.0),
        ("秋の日中(11月12時)", 24.7, 2.0, 288.0, 511.0),
        ("春の日中(5月12時)",  28.0, 2.2, 500.0, 420.0),
        ("夜間    (11月 2時)", 13.0, 1.4,   0.0, 800.0),
    ]

    print("-" * 78)
    print(f"{'場面':<20}{'蒸散':>12}{'蒸散':>14}{'純光合成':>14}{'糖換算':>14}")
    print(f"{'':<20}{'mg/m²s':>12}{'L/m²/h':>14}{'μmol/m²s':>14}{'g/m²/h':>14}")
    print("-" * 78)

    checks = []
    for label, temp, ea, solar, co2 in scenarios:
        r = calculate_canopy_gas_exchange(
            air_temp_c=temp,
            vapor_pressure_kpa=ea,
            solar_above_w_per_m2=solar,
            air_co2_ppm=co2,
            wind_speed_m_per_s=WIND_SPEED_M_PER_S,
            pressure_kpa=LOCAL_PRESSURE_KPA,
            characteristic_length_m=LEAF_CHARACTERISTIC_LENGTH_M,
            lai=lai,
        )

        liters_per_hour = r.transpiration_mg_per_m2_s * 3600 / 1e6
        micromol = r.net_photosynthesis_mmol_per_m2_s * 1000
        sugar_per_hour = (
            photosynthesis_to_sugar_g_per_m2_s(r.net_photosynthesis_mmol_per_m2_s)
            * 3600
        )

        print(
            f"{label:<20}{r.transpiration_mg_per_m2_s:>12.2f}"
            f"{liters_per_hour:>14.4f}{micromol:>14.2f}{sugar_per_hour:>14.3f}"
        )

        if solar > 0:
            # 日中の群落光合成は、施設トマトで 5〜35 μmol/(m²·s) 程度
            checks.append((
                f"{label} の光合成がありうる範囲（0〜40 μmol/m²s）",
                0.0 < micromol < 40.0,
                f"{micromol:.2f} μmol/m²s",
            ))
        else:
            # 夜間は呼吸のみなので負になる
            checks.append((
                f"{label} は呼吸のみで負",
                micromol < 0.0,
                f"{micromol:.2f} μmol/m²s",
            ))

    print("-" * 78)
    print()

    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '失敗'}  {label}")
        print(f"        {detail}")
    print()

    return all_passed


def test_seasonal_leaf_count() -> bool:
    """作業計画から日付ごとの葉枚数を引けるか。

    ここを年間固定（18枚/m²）にしていたのが、春の蒸散を過小評価する
    原因だった。実際は冬16枚〜春40枚/m² と倍以上変わる。
    """
    import datetime as dt

    from config import LEAF_COUNT_OUT_OF_SEASON_PER_M2
    from core.canopy import leaf_area_index_for_date, leaf_count_per_m2

    print("=" * 78)
    print("5. 葉枚数の季節変化（作業計画シートの取り込み）")
    print("=" * 78)

    results = {}

    print("  日付ごとの葉枚数と LAI:")
    checkpoints = [
        (dt.date(2025, 10, 15), 18.0, "秋（定植後）"),
        (dt.date(2025, 12, 15), 16.0, "厳冬期（減らす）"),
        (dt.date(2026, 2, 15), 28.0, "増やしはじめ"),
        (dt.date(2026, 4, 15), 40.0, "最大"),
        (dt.date(2026, 5, 15), 40.0, "最大"),
    ]
    matched = True
    for date, expected, label in checkpoints:
        actual = leaf_count_per_m2(date)
        lai = leaf_area_index_for_date(date)
        mark = "OK" if abs(actual - expected) < 0.1 else "NG"
        if mark == "NG":
            matched = False
        print(f"    [{mark}] {date}（{label}）"
              f" {actual:5.1f} 枚/m²  LAI {lai:.2f}")
    results["計画どおりの枚数が引ける"] = matched

    # 冬と春で倍以上ちがう
    winter = leaf_count_per_m2(dt.date(2025, 12, 15))
    spring = leaf_count_per_m2(dt.date(2026, 5, 15))
    ratio = spring / winter
    print(f"\n  春 ÷ 冬 = {spring:.0f} / {winter:.0f} = {ratio:.2f} 倍")
    results["春は冬の2倍以上の葉がある"] = ratio >= 2.0

    # 週の切れ目でなめらかにつながるか
    worst_jump = 0.0
    previous = None
    for offset in range(120):
        date = dt.date(2026, 1, 1) + dt.timedelta(days=offset)
        value = leaf_count_per_m2(date)
        if previous is not None:
            worst_jump = max(worst_jump, abs(value - previous))
        previous = value
    print(f"  隣り合う日の差の最大 {worst_jump:.2f} 枚/m²（1〜4月）")
    results["日ごとになめらかに変わる（1日の差が1枚未満）"] = worst_jump < 1.0

    # 作期外
    out_of_season = leaf_count_per_m2(dt.date(2026, 8, 15))
    print(f"  作期外（8月15日・週33）: {out_of_season:.1f} 枚/m²")
    results["作期外は既定値を返す"] = (
        abs(out_of_season - LEAF_COUNT_OUT_OF_SEASON_PER_M2) < 1e-9
    )

    print()
    for label, passed in results.items():
        print(f"  [{'OK' if passed else 'NG'}] {label}")
    print()
    return all(results.values())


def main() -> int:
    results = [
        test_layer_convergence(),
        test_beer_law(),
        test_lai_response(),
        test_realistic_magnitude(),
        test_seasonal_leaf_count(),
    ]

    print("=" * 78)
    if all(results):
        print("第3段階（群落スケーリング）は完了です。")
        print("=" * 78)
        return 0

    print("失敗した項目があります。先に進む前に修正が必要です。")
    print("=" * 78)
    return 1


if __name__ == "__main__":
    sys.exit(main())
