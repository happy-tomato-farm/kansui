"""
塩類濃度（EC）の検証。

【何を確かめるか】
1. ECの3つの基準の換算が往復するか
2. ★浸透ポテンシャルがマトリックポテンシャルより一桁大きいこと
3. 吸水の駆動力への効き方
4. 塩害の応答（収量と糖度が表と裏になっていること）
5. 肥料のN濃度から給液ECへ
6. 塩類収支が閉じるか
7. 定常状態の濃縮倍率
8. ★デルフィの指導の筋が通っているか（冬に高EC・春に低ECの意味）
9. 入力の検査が効くか
10. 土壌診断の実測値を通したときの桁
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (                                            # noqa: E402
    FEED_EC_BAND_BY_MONTH,
    FEED_EC_CONCENTRATION_FACTOR_ASSUMED,
    FIELD_CAPACITY_POTENTIAL_J_PER_KG,
    OPERATED_FEED_EC_BY_MONTH,
    RAW_WATER_EC_DS_PER_M,
    ROOT_ZONE_DEPTH_M,
    SALINITY_YIELD_THRESHOLD_ECE,
    SOIL_BULK_DENSITY_SENSITIVITY,
    SOIL_DIAGNOSIS_CHUO,
    SOIL_DIAGNOSIS_SAMPLE_DEPTH_M,
    SOIL_RETENTION,
)
from core.salinity import (                                      # noqa: E402
    DailySaltInput,
    clamp_daily_n_to_ec_band,
    daily_n_from_feed_ec,
    ec_1to5_to_saturated_extract,
    n_concentration_from_feed_ec,
    ec_1to5_to_soil_solution,
    feed_ec_from_n_concentration,
    gravimetric_water_content,
    mineral_n_to_g_per_m2,
    osmotic_potential_j_per_kg,
    salinity_response,
    simulate_salt,
    soil_solution_to_saturated_extract,
    steady_state_concentration_factor,
    total_soil_potential_j_per_kg,
    uptake_driving_force_ratio,
)
from core.soil import pf_from_potential, water_content_from_potential  # noqa: E402

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


theta_fc = water_content_from_potential(FIELD_CAPACITY_POTENTIAL_J_PER_KG)
matric_fc = FIELD_CAPACITY_POTENTIAL_J_PER_KG

# =============================================================================
print("=" * 80)
print("1. ECの3つの基準の換算")
print("=" * 80)
print(f"  圃場容水量の体積含水率 {theta_fc:.3f} → 重量含水率 "
      f"{gravimetric_water_content(theta_fc):.3f} g/g")
print(f"  飽和の体積含水率 {SOIL_RETENTION['THETA_S']:.3f} → 重量含水率 "
      f"{gravimetric_water_content(SOIL_RETENTION['THETA_S']):.3f} g/g")
print()
print(f"{'EC(1:5)':>9}{'土壌溶液EC':>12}{'ECe':>8}{'溶液→ECe':>11}")
worst = 0.0
for ec5 in (0.1, 0.18, 0.2, 0.47, 1.0):
    solution = ec_1to5_to_soil_solution(ec5, theta_fc)
    ece = ec_1to5_to_saturated_extract(ec5)
    # 溶液 → ECe の直接換算と、1:5 から直に出した ECe が合うこと
    direct = soil_solution_to_saturated_extract(solution, theta_fc)
    worst = max(worst, abs(direct - ece))
    print(f"{ec5:>9.2f}{solution:>12.2f}{ece:>8.2f}{direct:>11.2f}")
check(worst < 1e-9, "2つの経路で出した ECe が一致する", f"最大差 {worst:.2e}")

print()
print("  ★仮比重（実測していない）が換算倍率にそのまま効く")
print(f"{'仮比重':>8}{'EC(1:5) 0.2 → 土壌溶液EC':>26}")
for density in SOIL_BULK_DENSITY_SENSITIVITY:
    value = ec_1to5_to_soil_solution(0.2, theta_fc, density)
    print(f"{density:>8.1f}{value:>26.2f}")
span = (ec_1to5_to_soil_solution(0.2, theta_fc, min(SOIL_BULK_DENSITY_SENSITIVITY))
        / ec_1to5_to_soil_solution(0.2, theta_fc, max(SOIL_BULK_DENSITY_SENSITIVITY)))
check(span < 1.3, "仮比重の幅による振れは1.3倍未満", f"{span:.3f} 倍")

# =============================================================================
print()
print("=" * 80)
print("2. ★浸透ポテンシャルはマトリックポテンシャルより一桁大きい")
print("=" * 80)
print(f"  圃場容水量（pF 1.8）のマトリックポテンシャル {matric_fc} J/kg")
print()
print(f"{'土壌溶液EC':>11}{'浸透 J/kg':>11}{'マトリックの何倍':>18}"
      f"{'同じ効きのpF':>14}")
for ec in (0.5, 1.0, 2.0, 2.14, 4.0, 5.0):
    osmotic = osmotic_potential_j_per_kg(ec)
    print(f"{ec:>11.2f}{osmotic:>11.1f}{osmotic / matric_fc:>17.1f}倍"
          f"{pf_from_potential(osmotic):>14.2f}")

check(osmotic_potential_j_per_kg(2.0) / matric_fc > 10.0,
      "EC 2 dS/m の浸透はマトリックの10倍を超える",
      f"{osmotic_potential_j_per_kg(2.0) / matric_fc:.1f} 倍")

# 合計は足し算になること
total = total_soil_potential_j_per_kg(matric_fc, 2.0)
check(abs(total - (matric_fc + osmotic_potential_j_per_kg(2.0))) < 1e-12,
      "合計ポテンシャル = マトリック + 浸透", f"{total:.1f} J/kg")

print()
print("  ★つまり十分に湿った土では、根が感じる水ポテンシャルはECで決まる。")
print("    第9章の式9.20・9.22 は含水率だけで吸水を決めているので、")
print("    この項が落ちている。")

# =============================================================================
print()
print("=" * 80)
print("3. 吸水の駆動力への効き")
print("=" * 80)
print("  駆動力 = |ψ葉| − |ψ土|。葉内水ポテンシャル −1.0 MPa（1000 J/kg）のとき。")
print()
print(f"{'土壌溶液EC':>11}{'駆動力比':>10}{'減り方':>9}")
for ec in (0.0, 0.5, 1.0, 2.0, 4.0, 5.35):
    ratio = uptake_driving_force_ratio(matric_fc, ec)
    print(f"{ec:>11.2f}{ratio:>10.3f}{(1 - ratio) * 100:>8.1f}%")

check(abs(uptake_driving_force_ratio(matric_fc, 0.0) - 1.0) < 1e-12,
      "EC 0 なら駆動力は変わらない")
check(uptake_driving_force_ratio(matric_fc, 2.0)
      < uptake_driving_force_ratio(matric_fc, 1.0),
      "ECが高いほど駆動力は小さい")

# 葉内水ポテンシャルが低いほど、ECの相対的な影響は小さい
tight = uptake_driving_force_ratio(matric_fc, 2.0, leaf_potential_j_per_kg=1500.0)
loose = uptake_driving_force_ratio(matric_fc, 2.0, leaf_potential_j_per_kg=700.0)
print()
print(f"  葉内水ポテンシャル −1.5 MPa なら {tight:.3f}、"
      f"−0.7 MPa なら {loose:.3f}")
check(tight > loose,
      "葉をより低くできる（しおれに強い）ほどECの影響は小さい",
      f"{tight:.3f} > {loose:.3f}")
print("  ※果実の水ポテンシャルは葉より高い（−0.35〜−0.75 MPa）ので、")
print("    果実の肥大のほうがECに弱い。高EC栽培で糖度が上がる機序。")

# =============================================================================
print()
print("=" * 80)
print("4. 塩害の応答（収量と糖度は同じ現象の表と裏）")
print("=" * 80)
print(f"  収量のしきい値 ECe {SALINITY_YIELD_THRESHOLD_ECE} dS/m")
print()
print(f"{'ECe':>7}{'収量比':>9}{'糖度':>9}")
previous_yield = None
previous_brix = None
monotone = True
for ece in (1.0, 2.0, 2.5, 3.0, 3.5, 4.2, 6.0):
    response = salinity_response(ece)
    print(f"{ece:>7.2f}{response.yield_ratio * 100:>8.1f}%"
          f"{response.brix_change:>+9.2f}")
    if previous_yield is not None:
        if response.yield_ratio > previous_yield + 1e-12:
            monotone = False
        if response.brix_change < previous_brix - 1e-12:
            monotone = False
    previous_yield = response.yield_ratio
    previous_brix = response.brix_change

check(monotone, "ECeが上がると収量は下がり、糖度は上がる（表と裏）")
check(abs(salinity_response(2.0).yield_ratio - 1.0) < 1e-12,
      "しきい値より下では収量は落ちない")
check(salinity_response(100.0).yield_ratio == 0.0,
      "収量比は負にならない")

# =============================================================================
print()
print("=" * 80)
print("5. 肥料のN濃度から給液EC")
print("=" * 80)
print("  トケル養液配合1号（N比 0.10）のとき。実績の桁は12月 111・5月 18 mg-N/L。")
print()
print(f"{'N濃度':>9}{'肥料':>9}{'給液EC':>9}")
print(f"{'mg-N/L':>9}{'g/L':>9}{'dS/m':>9}")
for n_mg in (18.0, 50.0, 111.0, 200.0):
    ec = feed_ec_from_n_concentration(n_mg, 0.10)
    print(f"{n_mg:>9.0f}{n_mg / 1000.0 / 0.10:>9.2f}{ec:>9.2f}")

# 原水のECが切片になっていること
check(abs(feed_ec_from_n_concentration(0.0, 0.10) - 0.15) < 1e-12,
      "N濃度ゼロなら給液ECは原水のECになる",
      f"{feed_ec_from_n_concentration(0.0, 0.10):.3f} dS/m")
# N濃度に比例すること（切片を引けば）
a = feed_ec_from_n_concentration(100.0, 0.10) - 0.15
b = feed_ec_from_n_concentration(200.0, 0.10) - 0.15
check(abs(b - 2 * a) < 1e-12, "N濃度に比例する", f"{a:.3f} → {b:.3f}")

print()
print("  ★逆向き（目標ECから N量を決める）。段6の画面が使う経路。")
print()
print(f"{'給液EC':>9}{'N濃度':>10}{'潅水 1.12':>12}{'潅水 5.88':>12}")
print(f"{'dS/m':>9}{'mg-N/L':>10}{'kg-N/10a':>12}{'kg-N/10a':>12}")
worst = 0.0
for ec in (0.44, 0.47, 1.02, 1.93, 2.50):
    n_mg = n_concentration_from_feed_ec(ec, 0.10)
    # 往復するか（EC → N濃度 → EC）
    back = feed_ec_from_n_concentration(n_mg, 0.10)
    worst = max(worst, abs(back - ec))
    print(f"{ec:>9.2f}{n_mg:>10.1f}"
          f"{daily_n_from_feed_ec(ec, 1.12, 0.10):>12.3f}"
          f"{daily_n_from_feed_ec(ec, 5.88, 0.10):>12.3f}")
check(worst < 1e-9, "EC → N濃度 → EC が往復する", f"最大差 {worst:.2e}")

# 原水のECより低い目標は肥料ゼロ
check(n_concentration_from_feed_ec(0.10, 0.10) == 0.0,
      "原水のECより低い目標なら肥料は入れない（0 を返す）",
      "水だけでそのECなので、それ以上は下げようがない")

# ECを固定すると、潅水量に比例してN量が増える
a = daily_n_from_feed_ec(1.0, 2.0, 0.10)
b = daily_n_from_feed_ec(1.0, 4.0, 0.10)
check(abs(b - 2 * a) < 1e-12,
      "★ECを固定すると潅水量に比例してN量が増える",
      f"潅水 2 → {a:.3f} ／ 4 → {b:.3f} kg-N/10a。"
      f"レシピ側が「倍率を固定すればECが一定になり、潅水量を変えたぶんだけ"
      f"N量が増減する」と書いているのと同じ作り")

print()
print("  実績の運用値（12月 EC 1.93・潅水 1.12 ／ 5月 EC 0.47・潅水 5.88）:")
winter_n = daily_n_from_feed_ec(1.93, 1.12, 0.10)
spring_n = daily_n_from_feed_ec(0.47, 5.88, 0.10)
print(f"    12月 {winter_n:.3f} kg-N/10a（実績 0.164）")
print(f"     5月 {spring_n:.3f} kg-N/10a（実績 0.140）")
check(abs(winter_n - 0.164) < 0.05 and abs(spring_n - 0.140) < 0.05,
      "EC目標から出したN量が実績と近い",
      f"12月 {winter_n:.3f} vs 0.164 ／ 5月 {spring_n:.3f} vs 0.140")

# =============================================================================
print()
print("=" * 80)
print("6. 塩類収支が閉じるか")
print("=" * 80)

days = [
    DailySaltInput(
        date=f"{i}日目", irrigation_mm=6.0, drainage_mm=3.0,
        feed_ec_ds_per_m=1.0, n_uptake_g_per_m2=0.3,
        water_content=theta_fc)
    for i in range(60)
]
results = simulate_salt(days, 2.0, ROOT_ZONE_DEPTH_M)
worst = max(abs(r.balance_residual) for r in results)
print(f"  60日回した。EC {results[0].start_ec_solution:.2f} → "
      f"{results[-1].end_ec_solution:.2f} dS/m")
check(worst < 1e-9, "日ごとの塩類収支が閉じる", f"最大残差 {worst:.2e}")
check(all(r.end_salt >= 0.0 for r in results), "塩の量が負にならない")

# 流亡ゼロなら溜まり続ける
no_drain = [
    DailySaltInput(date=f"{i}日目", irrigation_mm=1.5, drainage_mm=0.0,
                   feed_ec_ds_per_m=1.9, n_uptake_g_per_m2=0.0,
                   water_content=theta_fc)
    for i in range(90)
]
piled = simulate_salt(no_drain, 2.0, ROOT_ZONE_DEPTH_M)
print(f"  流亡ゼロで90日（冬の条件）: EC 2.00 → "
      f"{piled[-1].end_ec_solution:.2f} dS/m")
check(piled[-1].end_ec_solution > 2.0,
      "流亡がなければ塩は溜まり続ける",
      f"{piled[-1].end_ec_solution:.2f} dS/m まで上がった")

# =============================================================================
print()
print("=" * 80)
print("7. 定常状態の濃縮倍率")
print("=" * 80)
print("  根圏EC ÷ 給液EC = 1 ÷ 流亡率")
print()
print(f"{'流亡率':>8}{'濃縮倍率':>10}{'給液1.0の落ち着き先':>22}")
for fraction in (0.5, 0.3, 0.2, 0.1):
    factor = steady_state_concentration_factor(fraction)
    print(f"{fraction * 100:>7.0f}%{factor:>10.1f}{factor:>22.1f}")
check(abs(steady_state_concentration_factor(0.5) - 2.0) < 1e-12,
      "流亡率50%なら2倍")
print()
print("  ★作物が塩を吸うぶんを数えていないので、これは上限。")

# =============================================================================
print()
print("=" * 80)
print("8. ★デルフィの指導の筋（冬に高EC・春に低EC）")
print("=" * 80)
print("  同じ駆動力の減り方でも、蒸散要求が大きい時期ほど痛い。")
print("  実績の月別（tools/時期別の妥当性を検証.py の出力）で確かめる。")
print()
print(f"{'時期':<10}{'根圏EC':>9}{'潜在蒸散':>10}{'駆動力の減り':>13}"
      f"{'失う蒸散':>10}")
print(f"{'':<10}{'dS/m':>9}{'mm/日':>10}{'%':>13}{'mm/日':>10}")
cases = [
    ("冬（12〜2月）", 4.74, 1.45),
    ("春（4〜6月）", 2.70, 4.88),
]
lost_by_season = {}
for label, ec, potential in cases:
    ratio = uptake_driving_force_ratio(matric_fc, ec)
    lost = potential * (1.0 - ratio)
    lost_by_season[label] = lost
    print(f"{label:<10}{ec:>9.2f}{potential:>10.2f}"
          f"{(1 - ratio) * 100:>12.1f}%{lost:>10.2f}")

winter_lost = lost_by_season["冬（12〜2月）"]
spring_lost = lost_by_season["春（4〜6月）"]
check(winter_lost < spring_lost,
      "★冬の高ECで失う蒸散は、春の低ECで失う量より少ない",
      f"冬 {winter_lost:.2f} < 春 {spring_lost:.2f} mm/日 "
      f"→ ECの負担を「水が要らない時期」に置く運用は筋が通っている")

# もし冬と春のECを入れ替えたらどうなるか
swapped_winter = 1.45 * (1.0 - uptake_driving_force_ratio(matric_fc, 2.70))
swapped_spring = 4.88 * (1.0 - uptake_driving_force_ratio(matric_fc, 4.74))
print()
print(f"  もし冬と春のECを入れ替えたら: 冬 {swapped_winter:.2f} ／ "
      f"春 {swapped_spring:.2f} mm/日")
print(f"  合計で失う蒸散  いまの運用 {winter_lost + spring_lost:.2f} ／ "
      f"入れ替え {swapped_winter + swapped_spring:.2f} mm/日")
check(winter_lost + spring_lost < swapped_winter + swapped_spring,
      "★いまの配り方のほうが、失う蒸散の合計が小さい",
      f"{winter_lost + spring_lost:.2f} < "
      f"{swapped_winter + swapped_spring:.2f} mm/日")

# =============================================================================
print()
print("=" * 80)
print("9. 入力の検査が効くか")
print("=" * 80)


def expect_error(label: str, call) -> None:
    try:
        call()
    except ValueError as error:
        check(True, label, str(error)[:52])
    else:
        check(False, label, "エラーにならなかった")


expect_error("ECが負なら止まる",
             lambda: osmotic_potential_j_per_kg(-1.0))
expect_error("含水率が範囲外なら止まる",
             lambda: ec_1to5_to_soil_solution(0.2, 1.5))
expect_error("仮比重が0なら止まる",
             lambda: ec_1to5_to_soil_solution(0.2, 0.42, 0.0))
expect_error("N比を % で渡したら止まる（0.10 であって 10 ではない）",
             lambda: feed_ec_from_n_concentration(100.0, 10.0))
expect_error("マトリックを負で渡したら止まる",
             lambda: total_soil_potential_j_per_kg(-1.0, 1.0))
expect_error("流亡率0では定常にならないので止まる",
             lambda: steady_state_concentration_factor(0.0))
expect_error("深さが0なら止まる",
             lambda: mineral_n_to_g_per_m2(5.0, 0.0))
expect_error("流亡が根圏の水量を超えたら止まる",
             lambda: simulate_salt(
                 [DailySaltInput(
                     date="x", irrigation_mm=500.0, drainage_mm=500.0,
                     feed_ec_ds_per_m=1.0, n_uptake_g_per_m2=0.0,
                     water_content=theta_fc)],
                 1.0, ROOT_ZONE_DEPTH_M))

# =============================================================================
print()
print("=" * 80)
print("10. 土壌診断の実測値を通したときの桁（中央ハウス）")
print("=" * 80)
print("  ★診断は6月採取なので『その作期の終わり』の値。")
print()
print(f"{'年':>6}{'EC(1:5)':>9}{'土壌溶液EC':>12}{'ECe':>7}{'浸透':>8}"
      f"{'駆動力比':>10}{'収量比':>8}{'無機態N':>9}")
print(f"{'':>6}{'dS/m':>9}{'dS/m':>12}{'dS/m':>7}{'J/kg':>8}{'':>10}"
      f"{'':>8}{'kg/10a':>9}")
for year, record in sorted(SOIL_DIAGNOSIS_CHUO.items()):
    ec5 = record["ec_1to5_ds_per_m"]
    solution = ec_1to5_to_soil_solution(ec5, theta_fc)
    ece = ec_1to5_to_saturated_extract(ec5)
    print(f"{year:>6}{ec5:>9.2f}{solution:>12.2f}{ece:>7.2f}"
          f"{osmotic_potential_j_per_kg(solution):>8.1f}"
          f"{uptake_driving_force_ratio(matric_fc, solution):>10.3f}"
          f"{salinity_response(ece).yield_ratio * 100:>7.1f}%"
          f"{mineral_n_to_g_per_m2(record['mineral_n_mg_per_100g'], SOIL_DIAGNOSIS_SAMPLE_DEPTH_M):>9.1f}")

# 2025年（EC 0.47）はしきい値を超えること
ece_2025 = ec_1to5_to_saturated_extract(
    SOIL_DIAGNOSIS_CHUO[2025]["ec_1to5_ds_per_m"])
check(ece_2025 > SALINITY_YIELD_THRESHOLD_ECE,
      "2025年の作終わりは塩害のしきい値を超えている",
      f"ECe {ece_2025:.2f} > {SALINITY_YIELD_THRESHOLD_ECE} dS/m")
ece_2026 = ec_1to5_to_saturated_extract(
    SOIL_DIAGNOSIS_CHUO[2026]["ec_1to5_ds_per_m"])
check(ece_2026 < SALINITY_YIELD_THRESHOLD_ECE,
      "2026年の作終わりはしきい値の下にある",
      f"ECe {ece_2026:.2f} < {SALINITY_YIELD_THRESHOLD_ECE} dS/m")

# =============================================================================
print()
print("=" * 80)
print("11. 給液ECの帯で挟む（量はNで決める）")
print("=" * 80)
print("  ★N需要から出したN量が、月別の帯から外れたら端まで戻す。")
print("    『ECを目標にしてN量を決める』のではない。向きが逆。")
print()

# --- 11-1. 帯の上限が ECe のしきい値と対応していること ---
# 給液EC × 濃縮倍率 → 根圏EC → ECe
for month, expected_ece in ((5, 2.5), (12, 4.0)):
    low, high = FEED_EC_BAND_BY_MONTH[month]
    ece_at_high = soil_solution_to_saturated_extract(
        high * FEED_EC_CONCENTRATION_FACTOR_ASSUMED, theta_fc)
    check(abs(ece_at_high - expected_ece) < 0.05,
          f"{month}月の帯の上限 {high:.2f} dS/m は ECe {expected_ece} に対応する",
          f"計算値 ECe {ece_at_high:.3f}（ねらい {expected_ece}）")

# 強光期の上限は低日射期より厳しいこと（尻腐果のため）
check(FEED_EC_BAND_BY_MONTH[5][1] < FEED_EC_BAND_BY_MONTH[12][1],
      "強光期（5月）の上限は低日射期（12月）より低い",
      f"5月 {FEED_EC_BAND_BY_MONTH[5][1]:.2f} < "
      f"12月 {FEED_EC_BAND_BY_MONTH[12][1]:.2f} dS/m")

# 9〜10月は盲点なので上げていないこと
for month in (9, 10):
    check(FEED_EC_BAND_BY_MONTH[month][1] <= FEED_EC_BAND_BY_MONTH[5][1],
          f"{month}月（盲点）の上限は強光期より上げていない",
          f"{month}月 {FEED_EC_BAND_BY_MONTH[month][1]:.2f} dS/m")

# すべての月で下限 ≦ 上限、かつ下限は原水ECより上（または栽培していない月）
for month, (low, high) in sorted(FEED_EC_BAND_BY_MONTH.items()):
    check(low <= high, f"{month}月の帯は下限 ≦ 上限",
          f"({low:.2f}, {high:.2f})")
    if high > 0.0:
        check(low > RAW_WATER_EC_DS_PER_M,
              f"{month}月の下限は原水EC {RAW_WATER_EC_DS_PER_M} より上",
              f"下限 {low:.2f} dS/m")

# --- 11-2. 帯の中なら素通りすること ---
n_fraction = 0.100
irrigation = 2.0   # L/m²
# 給液EC 0.8 dS/m になるN量を作って渡す
n_inside = daily_n_from_feed_ec(0.80, irrigation, n_fraction)
inside = clamp_daily_n_to_ec_band(n_inside, irrigation, 5, n_fraction)
check(not inside.is_clamped, "帯の中のN量は挟まれない",
      f"給液EC {inside.feed_ec_ds_per_m:.2f} dS/m / {inside.describe()}")
check(abs(inside.daily_n_kg_per_10a - n_inside) < 1e-12,
      "帯の中ならN量は変わらない",
      f"{n_inside:.4f} → {inside.daily_n_kg_per_10a:.4f} kg-N/10a")

# --- 11-3. 上限で切られること（5月に濃い液を渡す）---
n_too_much = daily_n_from_feed_ec(1.60, irrigation, n_fraction)
high_hit = clamp_daily_n_to_ec_band(n_too_much, irrigation, 5, n_fraction)
check(high_hit.clamped_by == "high", "5月に濃すぎるN量は上限で切られる",
      high_hit.describe())
check(high_hit.daily_n_kg_per_10a < n_too_much,
      "上限で切られたらN量は減る",
      f"{n_too_much:.4f} → {high_hit.daily_n_kg_per_10a:.4f} kg-N/10a")
check(abs(high_hit.feed_ec_ds_per_m - FEED_EC_BAND_BY_MONTH[5][1]) < 1e-9,
      "切ったあとの給液ECはちょうど上限",
      f"{high_hit.feed_ec_ds_per_m:.3f} dS/m")

# 同じN量を12月に渡すと、帯が広いので切られないこと
dec = clamp_daily_n_to_ec_band(n_too_much, irrigation, 12, n_fraction)
check(not dec.is_clamped, "同じN量を12月に渡すと帯の中に入る", dec.describe())

# --- 11-4. 下限で持ち上げられること ---
n_too_little = daily_n_from_feed_ec(0.20, irrigation, n_fraction)
low_hit = clamp_daily_n_to_ec_band(n_too_little, irrigation, 5, n_fraction)
check(low_hit.clamped_by == "low", "薄すぎるN量は下限で持ち上げられる",
      low_hit.describe())
check(low_hit.daily_n_kg_per_10a > n_too_little,
      "下限で持ち上げたらN量は増える",
      f"{n_too_little:.4f} → {low_hit.daily_n_kg_per_10a:.4f} kg-N/10a")

# --- 11-5. 潅水量が 0 の日は挟まずに通すこと（割り算ができない）---
zero_water = clamp_daily_n_to_ec_band(0.15, 0.0, 5, n_fraction)
check(not zero_water.is_clamped and
      abs(zero_water.daily_n_kg_per_10a - 0.15) < 1e-12,
      "潅水量 0 の日は挟まずにそのまま通す",
      f"N量 {zero_water.daily_n_kg_per_10a:.3f} kg-N/10a・"
      f"給液EC {zero_water.feed_ec_ds_per_m:.2f} dS/m")

# --- 11-6. 栽培していない月（帯が 0〜0）は液肥を出さないこと ---
off_season = clamp_daily_n_to_ec_band(0.15, irrigation, 7, n_fraction)
check(off_season.daily_n_kg_per_10a == 0.0,
      "7月（栽培していない）は液肥を出さない",
      f"N量 {off_season.daily_n_kg_per_10a:.3f} kg-N/10a")

# --- 11-7. 異常な入力はエラーになること（黙って既定値にしない）---
for label, args in (
    ("N量が負", (-0.1, irrigation, 5, n_fraction)),
    ("潅水量が負", (0.15, -1.0, 5, n_fraction)),
    ("帯に無い月", (0.15, irrigation, 13, n_fraction)),
):
    try:
        clamp_daily_n_to_ec_band(*args)
    except ValueError as err:
        check(True, f"{label} は ValueError になる", str(err).split("。")[0])
    else:
        check(False, f"{label} は ValueError になる", "例外が出なかった")

# --- 11-8. 12月の運用値は帯の上限を超えていること（いまの管理より下げる向き）---
check(OPERATED_FEED_EC_BY_MONTH[12] > FEED_EC_BAND_BY_MONTH[12][1],
      "12月の運用値は帯の上限を超えている（帯は下げる向きに効く）",
      f"運用 {OPERATED_FEED_EC_BY_MONTH[12]:.2f} > "
      f"上限 {FEED_EC_BAND_BY_MONTH[12][1]:.2f} dS/m。"
      f"★根圏ECの濃縮倍率は実測していないので暫定")

# 月別の帯と運用値を並べて見せる
print(f"{'月':>4}{'下限':>8}{'上限':>8}{'運用値':>9}{'上限でのECe':>12}{'判定':>14}")
for month in (9, 10, 11, 12, 1, 2, 3, 4, 5, 6):
    low, high = FEED_EC_BAND_BY_MONTH[month]
    operated = OPERATED_FEED_EC_BY_MONTH[month]
    ece_high = soil_solution_to_saturated_extract(
        high * FEED_EC_CONCENTRATION_FACTOR_ASSUMED, theta_fc)
    verdict = "運用が上限超え" if operated > high else "運用は帯の中"
    print(f"{month:>4}{low:>8.2f}{high:>8.2f}{operated:>9.2f}"
          f"{ece_high:>12.2f}{verdict:>14}")

# =============================================================================
print()
print("=" * 80)
if failures:
    print(f"NG が {len(failures)} 件ある:")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("すべて通過")
print("=" * 80)
