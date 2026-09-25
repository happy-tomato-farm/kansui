"""
土壌モジュールの検証。

【何を確かめるか】
1. 水分特性曲線が、土性の調査.xlsx の実測13点を再現できているか
2. 含水率 ⇔ 水ポテンシャルの往復変換が元に戻るか
3. 透水係数が単調で、圃場容水量の決め方が Ks の仮定に左右されないか
4. 式9.20・式9.22 が教科書どおりの値を返すか
5. 水収支が閉じているか（潅水 = 蒸散 + 蒸発 + 排水 + 貯留変化）
6. ユーザーの仮説「春先に圃場容水量以上を保てば蒸散が増える」の検証
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (
    FIELD_CAPACITY_ALTERNATIVE_J_PER_KG,
    FIELD_CAPACITY_POTENTIAL_J_PER_KG,
    MIN_AIR_FILLED_POROSITY,
    NEGLIGIBLE_DRAINAGE_MM_PER_DAY,
    ROOT_SYSTEM_MAX_UPTAKE_SENSITIVITY,
    ROOT_ZONE_DEPTH_SENSITIVITY_M,
    DRIP_WETTED_FRACTION_SENSITIVITY,
    SOIL_RETENTION,
)
from core import soil
from core.water_balance import (
    DailyWaterInput,
    WaterBalanceSettings,
    simulate_water_balance,
    summarize,
)

# 土性の調査.xlsx シルト質壌土（B〜D列）の実測値。
# B列の数値は表記が「J/kg」だが中身は cm H2O なので、換算してから使う。
# （根拠は config.py 第7節の注意書き）
MEASURED_HEAD_CM = [
    (1.0, 0.502), (5.0, 0.497), (10.0, 0.484), (20.0, 0.458), (40.0, 0.434),
    (80.0, 0.414), (160.0, 0.406), (345.0, 0.395), (690.0, 0.391),
    (2000.0, 0.354), (5000.0, 0.304), (10000.0, 0.268), (15000.0, 0.248),
]
MEASURED = [
    (head_cm * soil.J_PER_KG_PER_CM_WATER, theta)
    for head_cm, theta in MEASURED_HEAD_CM
]

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


# =============================================================================
print("=" * 78)
print("1. 水分特性曲線が実測値を再現するか")
print("=" * 78)
print("   h[cm]   ψ[J/kg]     pF     実測θ    計算θ    残差")

max_residual = 0.0
for (head_cm, theta_obs), (psi, _) in zip(MEASURED_HEAD_CM, MEASURED):
    theta_fit = soil.water_content_from_potential(psi)
    residual = theta_fit - theta_obs
    max_residual = max(max_residual, abs(residual))
    print(f"  {head_cm:7.0f} {psi:9.1f}  {soil.pf_from_potential(psi):5.2f}  "
          f"{theta_obs:.4f}  {theta_fit:.4f}  {residual:+.4f}")

rmse = (sum((soil.water_content_from_potential(p) - t) ** 2 for p, t in MEASURED)
        / len(MEASURED)) ** 0.5
print(f"\n  RMSE = {rmse:.5f}   最大残差 = {max_residual:.5f}")
check(rmse < 0.003, "RMSE が 0.003 未満", f"{rmse:.5f}")
check(max_residual < 0.005, "全点の残差が ±0.005 以内", f"最大 {max_residual:.5f}")

# =============================================================================
print()
print("=" * 78)
print("2. 含水率 ⇔ 水ポテンシャルの往復変換")
print("=" * 78)

worst = 0.0
for psi, _ in MEASURED:
    theta = soil.water_content_from_potential(psi)
    psi_back = soil.potential_from_water_content(theta)
    relative_error = abs(psi_back - psi) / psi
    worst = max(worst, relative_error)
print(f"  ψ → θ → ψ の最大相対誤差 = {worst:.2e}")
check(worst < 1e-6, "往復して元の値に戻る", f"{worst:.2e}")

# pF との往復
worst_pf = max(abs(soil.pf_from_potential(soil.potential_from_pf(pf)) - pf)
               for pf in [1.0, 1.8, 2.5, 3.0, 4.2])
check(worst_pf < 1e-9, "pF ⇔ ψ の往復が一致", f"{worst_pf:.2e}")

# 単調性（θ が増えれば ψ は減る）
thetas = [0.25 + 0.01 * i for i in range(25)]
potentials = [soil.potential_from_water_content(t) for t in thetas]
check(all(a > b for a, b in zip(potentials, potentials[1:])),
      "含水率が増えると水ポテンシャルは単調に減る")

# 範囲外はエラーにする（黙って既定値に置き換えない）
try:
    soil.potential_from_water_content(SOIL_RETENTION["THETA_R"] - 0.01)
    check(False, "残留含水率以下でエラーになる")
except ValueError as error:
    check("残留含水率" in str(error), "残留含水率以下でエラーになる", str(error)[:40] + "…")

# =============================================================================
print()
print("=" * 78)
print("3. 透水係数と圃場容水量")
print("=" * 78)
print("  pF    ψ[J/kg]     θ     空気率   K[mm/日]")

for pf in [1.0, 1.5, 1.8, 2.0, 2.2, 2.54, 2.7, 3.0, 3.5, 4.18]:
    psi = soil.potential_from_pf(pf)
    theta = soil.water_content_from_potential(psi)
    print(f"  {pf:4.2f} {psi:9.1f}  {theta:.4f}  {soil.air_filled_porosity(theta) * 100:5.1f}%  "
          f"{soil.hydraulic_conductivity_mm_per_day(theta):9.3f}")

ks_values = [soil.hydraulic_conductivity_mm_per_day(soil.water_content_from_potential(p))
             for p in [1.0, 10.0, 100.0, 1000.0]]
check(all(a > b for a, b in zip(ks_values, ks_values[1:])),
      "乾くほど透水係数が小さくなる")

fc = soil.FIELD_CAPACITY_WATER_CONTENT
pwp = soil.WILTING_POINT_WATER_CONTENT
print(f"\n  圃場容水量 θfc = {fc:.4f} "
      f"(pF {soil.pf_from_potential(FIELD_CAPACITY_POTENTIAL_J_PER_KG):.2f}, 空気率 "
      f"{soil.air_filled_porosity(fc) * 100:.1f}%)")
print(f"  永久しおれ点 θpwp = {pwp:.4f} (pF {soil.pf_from_potential(1500.0):.2f})")
print(f"  排水が {NEGLIGIBLE_DRAINAGE_MM_PER_DAY} mm/日 まで落ちる θ = "
      f"{soil.DRAINAGE_STOP_WATER_CONTENT:.4f} "
      f"(pF {soil.pf_from_potential(soil.potential_from_water_content(soil.DRAINAGE_STOP_WATER_CONTENT)):.2f})")
print(f"  有効水分量 θfc − θpwp = {fc - pwp:.4f}  → 根群域30cmで {(fc - pwp) * 300:.0f} mm")

# 圃場容水量は「排水が実質止まる含水率」。両者が近ければ定義と整合している。
check(abs(soil.DRAINAGE_STOP_WATER_CONTENT - fc) < 0.02,
      "採用した圃場容水量と「排水が止まる含水率」がほぼ一致する",
      f"差 {abs(soil.DRAINAGE_STOP_WATER_CONTENT - fc):.4f}")

# 文献のシルト質壌土に収まっているか（単位の読み違いを検出する網）
check(0.13 <= fc - pwp <= 0.20,
      "有効水分量がシルト質壌土の文献範囲 0.13〜0.20 に入る", f"{fc - pwp:.3f}")
check(0.22 <= pwp <= 0.28,
      "永久しおれ点が実測の最終点（θ=0.248）付近に来る", f"{pwp:.3f}")

# 国際慣用の 1/3 bar を使った場合との比較
theta_third_bar = soil.water_content_from_potential(FIELD_CAPACITY_ALTERNATIVE_J_PER_KG)
print(f"\n  参考: 国際慣用 1/3 bar（{FIELD_CAPACITY_ALTERNATIVE_J_PER_KG} J/kg, pF "
      f"{soil.pf_from_potential(FIELD_CAPACITY_ALTERNATIVE_J_PER_KG):.2f}）を使った場合")
print(f"    θ = {theta_third_bar:.4f}、空気率 "
      f"{soil.air_filled_porosity(theta_third_bar) * 100:.1f}%、"
      f"有効水分 {theta_third_bar - pwp:.4f}、"
      f"排水 {soil.hydraulic_conductivity_mm_per_day(theta_third_bar):.3f} mm/日")
check(soil.hydraulic_conductivity_mm_per_day(theta_third_bar)
      < NEGLIGIBLE_DRAINAGE_MM_PER_DAY,
      "1/3 bar ではこの土の排水はとうに止まっている（乾かしすぎ側）",
      f"{soil.hydraulic_conductivity_mm_per_day(theta_third_bar):.3f} mm/日")
check(abs((fc - pwp) - (theta_third_bar - pwp)) / (fc - pwp) < 0.20,
      "どちらの流儀を採っても有効水分量の差は2割以内",
      f"{(fc - pwp):.3f} vs {(theta_third_bar - pwp):.3f}")

# =============================================================================
print()
print("=" * 78)
print("4. 式9.20（有効水分度）と 式9.22（可能吸水速度）")
print("=" * 78)
print("   Aw     Up*    教科書の値")

# 教科書 式9.22: Up* = 1 − (1 + 1.37·Aw)^(−5)
# 期待値は式を電卓で直接計算したもの（実装とは別経路で求める）。
for aw in [1.0, 0.5, 0.3, 0.2, 0.1, 0.0]:
    want = 1.0 - 1.0 / (1.0 + 1.37 * aw) ** 5
    got = soil.potential_uptake_ratio(aw)
    print(f"  {aw:.1f}  {got:.4f}  {want:.4f}")
    check(abs(got - want) < 1e-9, f"Up*(Aw={aw}) が式どおり", f"{got:.4f} vs {want:.4f}")

check(abs(soil.available_water_fraction(soil.FIELD_CAPACITY_WATER_CONTENT) - 1.0) < 1e-9,
      "圃場容水量で Aw = 1")
check(abs(soil.available_water_fraction(soil.WILTING_POINT_WATER_CONTENT)) < 1e-9,
      "永久しおれ点で Aw = 0")
check(soil.available_water_fraction(SOIL_RETENTION["THETA_S"]) == 1.0,
      "圃場容水量を超えても Aw は 1 で頭打ち（余分な水は吸水能力を増やさない）")

# =============================================================================
print()
print("=" * 78)
print("5. 水収支が閉じているか")
print("=" * 78)

days = [DailyWaterInput(f"2026-04-{d:02d}", irrigation_mm=4.0,
                        potential_transpiration_mm=3.0) for d in range(1, 31)]
results = simulate_water_balance(days)
worst_residual = max(abs(r.balance_residual_mm) for r in results)
print(f"  30日ぶん計算。収支の残差の最大 = {worst_residual:.3e} mm")
check(worst_residual < 1e-9, "毎日 潅水 = 蒸散 + 蒸発 + 排水 + 貯留変化 が成り立つ",
      f"{worst_residual:.2e} mm")

total = summarize(results)
print(f"  潅水 {total.irrigation_mm:.1f} mm / 蒸散 {total.actual_transpiration_mm:.1f} mm "
      f"/ 排水 {total.drainage_mm:.1f} mm / 貯留変化 {total.storage_change_mm:+.1f} mm")
print(f"  流亡率 {total.drainage_fraction * 100:.1f}%")
check(total.drainage_fraction >= 0.0, "流亡率が負にならない")

# 入力の検査が効いているか
for bad, label in [
    ({"wetted_fraction": 0.0}, "濡れ面積率 0 を弾く"),
    ({"wetted_fraction": 1.5}, "濡れ面積率 1.5 を弾く"),
    ({"root_zone_depth_m": -0.1}, "根群域の深さが負なら弾く"),
    ({"daytime_window_hour": (18.0, 6.0)}, "時間帯が逆順なら弾く"),
]:
    try:
        WaterBalanceSettings(**bad)
        check(False, label)
    except ValueError:
        check(True, label)

# =============================================================================
print()
print("=" * 78)
print("6. 仮説の検証: 「土を常に圃場容水量以上に保てば蒸散が増えるか」")
print("=" * 78)
print("  春先を想定（可能蒸散 3.5 mm/日）。潅水量だけを変えて30日回す。")
print()
print("  潅水    蒸散    流亡    流亡率  酸欠時間  FC超過時間  平均pF")
print("  mm/日   mm/日   mm/日          時間/日   時間/日")

hypothesis_rows = []
for irrigation in [2.0, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]:
    days = [DailyWaterInput(f"day{d}", irrigation, 3.5) for d in range(30)]
    res = simulate_water_balance(days)
    agg = summarize(res)
    # 最後の10日（初期条件の影響が抜けたあと）で評価する
    tail = res[-10:]
    mean_pf = sum((r.start_pf + r.end_pf) / 2 for r in tail) / len(tail)
    transpiration_per_day = sum(r.actual_transpiration_mm for r in tail) / len(tail)
    drainage_per_day = sum(r.drainage_mm for r in tail) / len(tail)
    anoxic_per_day = sum(r.hours_anoxic for r in tail) / len(tail)
    fc_per_day = sum(r.hours_above_field_capacity for r in tail) / len(tail)
    fraction = drainage_per_day / irrigation
    hypothesis_rows.append((irrigation, transpiration_per_day, anoxic_per_day))
    print(f"  {irrigation:5.1f}  {transpiration_per_day:6.2f}  {drainage_per_day:6.2f}  "
          f"{fraction * 100:5.1f}%  {anoxic_per_day:7.1f}  {fc_per_day:9.1f}   {mean_pf:.2f}")

wet = [r for r in hypothesis_rows if r[0] >= 4.0]
transpiration_spread = max(r[1] for r in wet) - min(r[1] for r in wet)
print(f"\n  潅水 4〜10 mm/日 のあいだで、蒸散量の差は {transpiration_spread:.3f} mm/日")
check(transpiration_spread < 0.05,
      "十分な潅水を超えたら、さらに増やしても蒸散は増えない",
      f"差 {transpiration_spread:.3f} mm/日")
check(max(r[2] for r in wet) > 0.0,
      "過剰な潅水は酸欠時間を生む",
      f"最大 {max(r[2] for r in wet):.1f} 時間/日")

# =============================================================================
print()
print("=" * 78)
print("7. 感度解析: 実測できていない前提を振ると結論はどう動くか")
print("=" * 78)
print("  （潅水 6 mm/日・可能蒸散 3.5 mm/日 の30日、最後の10日で評価）")
print()
print("  根群域  濡れ面積  最大吸水   蒸散      流亡率   平均pF")
print("   m        率      mm/日     mm/日")


def run_case(depth: float, wetted: float, uptake: float,
             irrigation: float = 6.0, demand: float = 3.5):
    settings = WaterBalanceSettings(
        root_zone_depth_m=depth, wetted_fraction=wetted, max_uptake_mm_per_day=uptake
    )
    res = simulate_water_balance(
        [DailyWaterInput(f"day{d}", irrigation, demand) for d in range(30)], settings
    )
    tail = res[-10:]
    return (
        sum(r.actual_transpiration_mm for r in tail) / len(tail),
        sum(r.drainage_mm for r in tail) / len(tail),
        sum((r.start_pf + r.end_pf) / 2 for r in tail) / len(tail),
    )


all_transpiration = []
geometry_transpiration = []   # 根群域・濡れ面積だけを振ったもの（最大吸水は既定）
for depth in ROOT_ZONE_DEPTH_SENSITIVITY_M:
    for wetted in DRIP_WETTED_FRACTION_SENSITIVITY:
        for uptake in ROOT_SYSTEM_MAX_UPTAKE_SENSITIVITY:
            transpiration, drainage, mean_pf = run_case(depth, wetted, uptake)
            all_transpiration.append(transpiration)
            if uptake == 6.0:
                geometry_transpiration.append(transpiration)
            print(f"  {depth:5.2f}   {wetted:5.2f}    {uptake:5.1f}    "
                  f"{transpiration:6.2f}    {drainage / 6.0 * 100:5.1f}%   {mean_pf:.2f}")

geometry_spread = max(geometry_transpiration) - min(geometry_transpiration)
total_spread = max(all_transpiration) - min(all_transpiration)
print(f"\n  根群域の深さ・濡れ面積率だけを振った幅 = {geometry_spread:.4f} mm/日")
print(f"  最大吸水速度も含めた幅               = {total_spread:.4f} mm/日")
check(geometry_spread < 1e-6,
      "根群域の深さと濡れ面積率は、定常状態の水の配分を変えない",
      f"幅 {geometry_spread:.2e} mm/日")

# 水が足りているうちは最大吸水速度も効かないはず（Up* ≒ 1 で min() が可能蒸散側を選ぶ）
print()
print("  水が足りている場合に、最大吸水速度が効くかどうか:")
for uptake in ROOT_SYSTEM_MAX_UPTAKE_SENSITIVITY:
    transpiration, _, _ = run_case(0.35, 0.7, uptake, irrigation=6.0, demand=3.5)
    limited = "律速する" if transpiration < 3.5 - 1e-6 else "律速しない"
    print(f"    最大吸水 {uptake:5.1f} mm/日 → 蒸散 {transpiration:.2f} mm/日  （{limited}）")

check(abs(run_case(0.35, 0.7, 6.0)[0] - 3.5) < 1e-6,
      "潅水が足りていれば、蒸散は可能蒸散どおりになる（土壌は律速しない）")

# 最大吸水速度が可能蒸散を下回れば、当然そこで頭打ちになる。
# 既定の 6 mm/日 は可能蒸散 3.5 mm/日 より十分大きいので効かないが、
# 真夏に可能蒸散が 5〜6 mm/日 まで上がるとこの値が効きはじめる。
capped = run_case(0.35, 0.7, 3.0, irrigation=6.0, demand=3.5)[0]
print(f"\n    参考: 最大吸水 3.0 mm/日（可能蒸散 3.5 より小さい）→ 蒸散 {capped:.2f} mm/日")
check(capped < 3.5 - 1e-6,
      "最大吸水速度が可能蒸散を下回ると、そこで頭打ちになる",
      f"{capped:.2f} mm/日")

# 潅水が足りなければ土壌水分が律速する
dry = run_case(0.35, 0.7, 6.0, irrigation=2.0, demand=3.5)
print(f"    参考: 潅水 2.0 mm/日（可能蒸散 3.5 より少ない）→ 蒸散 {dry[0]:.2f} mm/日、"
      f"平均 pF {dry[2]:.2f}")
check(dry[0] < 3.5 - 1e-6, "潅水が不足すれば土壌水分が蒸散を抑える", f"{dry[0]:.2f} mm/日")

# --- この土に固有の、見逃せない事実 ---
print()
print("  この土の圃場容水量における空気率:")
air_at_fc = soil.air_filled_porosity(soil.FIELD_CAPACITY_WATER_CONTENT)
print(f"    θfc = {soil.FIELD_CAPACITY_WATER_CONTENT:.4f} のとき空気率 {air_at_fc * 100:.1f}%")
print(f"    根が必要とする目安は {MIN_AIR_FILLED_POROSITY * 100:.0f}% 以上")
theta_for_air = SOIL_RETENTION["THETA_S"] - MIN_AIR_FILLED_POROSITY
pf_for_air = soil.pf_from_potential(soil.potential_from_water_content(theta_for_air))
print(f"    空気率 {MIN_AIR_FILLED_POROSITY * 100:.0f}% を確保するには θ ≦ {theta_for_air:.4f}"
      f"（pF {pf_for_air:.2f}）まで乾かす必要がある")
check(air_at_fc < MIN_AIR_FILLED_POROSITY,
      "この土は圃場容水量でもすでに空気率の目安を下回る（＝水をやりすぎる余地が小さい）",
      f"{air_at_fc * 100:.1f}% < {MIN_AIR_FILLED_POROSITY * 100:.0f}%")

# =============================================================================
print()
print("=" * 78)
print("8. 水ポテンシャル勾配: 「土を湿らせれば蒸散が増える」は成り立つか")
print("=" * 78)
print("  日中の葉内水ポテンシャルを −1.0 MPa（1000 J/kg）としたときの駆動力")
print()
print("  状態                    pF    θ      ψ土[J/kg]  駆動力[J/kg]  飽和時比")

LEAF_POTENTIAL = 1000.0
reference = soil.soil_to_leaf_driving_force(SOIL_RETENTION["THETA_S"], LEAF_POTENTIAL)
driving_rows = []
for label, pf in [("ほぼ飽和 pF1.0", 1.0), ("圃場容水量 pF1.8", 1.8),
                  ("1/3bar pF2.54", 2.54),
                  ("pF3.0", 3.0), ("pF3.5", 3.5), ("pF4.0", 4.0)]:
    psi = soil.potential_from_pf(pf)
    theta = soil.water_content_from_potential(psi)
    force = soil.soil_to_leaf_driving_force(theta, LEAF_POTENTIAL)
    driving_rows.append((pf, force))
    print(f"  {label:<22} {pf:4.2f}  {theta:.4f}  {psi:9.1f}  {force:10.1f}  "
          f"{force / reference * 100:6.1f}%")

force_saturated = soil.soil_to_leaf_driving_force(
    soil.water_content_from_potential(soil.potential_from_pf(1.0)), LEAF_POTENTIAL)
force_fc = soil.soil_to_leaf_driving_force(
    soil.FIELD_CAPACITY_WATER_CONTENT, LEAF_POTENTIAL)
gain = (force_saturated - force_fc) / force_fc
print(f"\n  圃場容水量 → ほぼ飽和 まで湿らせて、駆動力の増加は {gain * 100:.1f}%")
check(gain < 0.05,
      "圃場容水量を超えて湿らせても、駆動力はほとんど増えない",
      f"+{gain * 100:.1f}%")

air_fc = soil.air_filled_porosity(soil.FIELD_CAPACITY_WATER_CONTENT)
air_sat = soil.air_filled_porosity(
    soil.water_content_from_potential(soil.potential_from_pf(1.0)))
print(f"  同じ範囲で空気率は {air_fc * 100:.1f}% → {air_sat * 100:.1f}% に下がる")
check(air_sat < air_fc / 2.0,
      "同じ範囲で空気率は半分以下になる（得るものより失うものが大きい）",
      f"{air_fc * 100:.1f}% → {air_sat * 100:.1f}%")

# =============================================================================
print()
print("=" * 78)
print("9. 同じ天気がつづいたときの落ち着き先（steady_state）")
print("=" * 78)
print("  1日だけの計算では、余った水は「土に溜まった」ことになって流亡に出ない。")
print("  毎日くり返すと土が満杯になり、入れた水がそのまま抜けるようになる。")

from core.soil import FIELD_CAPACITY_WATER_CONTENT
from core.water_balance import step_one_day, steady_state

TRANSPIRATION_MM = 2.91      # 5月・平年日射・葉18枚/m² のときの蒸散
settings = WaterBalanceSettings()

print()
print("   潅水[mm]  1日目の流亡率   定常の流亡率   落ち着くまで   落ち着いたか")
previous_fraction = -1.0
monotone = True
for irrigation in (3.1, 4.1, 5.2, 6.2, 8.3):
    first = step_one_day(
        FIELD_CAPACITY_WATER_CONTENT,
        DailyWaterInput(date="1", irrigation_mm=irrigation,
                        potential_transpiration_mm=TRANSPIRATION_MM),
        settings,
    )
    steady = steady_state(irrigation, TRANSPIRATION_MM, settings)
    print(f"   {irrigation:7.1f}  {first.drainage_fraction * 100:12.0f}%"
          f"  {steady.balance.drainage_fraction * 100:12.0f}%"
          f"  {steady.days:11d} 日  {'落ち着く' if steady.settled else '落ち着かない':>12s}")
    if steady.balance.drainage_fraction < previous_fraction:
        monotone = False
    previous_fraction = steady.balance.drainage_fraction

check(monotone, "潅水を増やすと定常の流亡率が単調に上がる")

# 定常では「入った水 ＝ 出ていく水」になる（土の貯留量が変わらない）
steady = steady_state(6.2, TRANSPIRATION_MM, settings)
balance = steady.balance
inflow = balance.irrigation_mm
outflow = (balance.actual_transpiration_mm + balance.drainage_mm
           + balance.soil_evaporation_mm)
check(steady.settled, "5月に 3 L/10MJ 相当（6.2mm）なら落ち着く",
      f"{steady.days} 日で落ち着く")
check(abs(inflow - outflow) < 0.02,
      "定常では 入った水 ＝ 蒸散 + 流亡（土に溜まる量がゼロになる）",
      f"入 {inflow:.3f} / 出 {outflow:.3f} mm")
check(abs(balance.storage_change_mm) < 0.02,
      "定常では土の貯留量の変化がほぼゼロ",
      f"{balance.storage_change_mm:+.4f} mm")

# 1日だけだと流亡が過小に出ることを、数字で押さえておく
first = step_one_day(
    FIELD_CAPACITY_WATER_CONTENT,
    DailyWaterInput(date="1", irrigation_mm=6.2,
                    potential_transpiration_mm=TRANSPIRATION_MM),
    settings,
)
check(first.drainage_fraction < balance.drainage_fraction / 5.0,
      "1日だけの計算は定常より流亡率がずっと小さく出る（5分の1未満）",
      f"1日目 {first.drainage_fraction * 100:.0f}% → "
      f"定常 {balance.drainage_fraction * 100:.0f}%")

# 実測から逆算した流亡率 55〜57% と合うか
check(0.45 < balance.drainage_fraction < 0.65,
      "3 L/10MJ の定常流亡率が、実測から逆算した 55〜57% の近くに来る",
      f"{balance.drainage_fraction * 100:.0f}%")

# 潅水が蒸散に大きく足りないと、しおれ点近くまで乾いてから落ち着く。
#
# 落ち着くのは「土が乾くと吸水が減る」（式9.22）ため。
# 水分が減るほど可能吸水速度が落ち、潅水量と釣り合うところで止まる。
# 落ち着いた＝健全、ではないので、pF を必ず見ること。
dry = steady_state(1.0, TRANSPIRATION_MM, settings, max_days=120)
print(f"\n  潅水 1.0 mm / 蒸散要求 {TRANSPIRATION_MM} mm のとき:")
print(f"    pF {dry.balance.end_pf:.2f} / 実際の蒸散 "
      f"{dry.balance.actual_transpiration_mm:.2f} mm / "
      f"{'落ち着く' if dry.settled else '落ち着かない'}（{dry.days} 日）")
check(dry.balance.end_pf > 3.5,
      "潅水が足りないと、しおれ点近く（pF 3.5超）まで乾く",
      f"pF {dry.balance.end_pf:.2f}（永久しおれ点は 4.18）")
check(dry.balance.actual_transpiration_mm < TRANSPIRATION_MM * 0.6,
      "そこでは蒸散が要求の6割未満に抑えられている（＝水ストレス）",
      f"{dry.balance.actual_transpiration_mm:.2f} / {TRANSPIRATION_MM} mm")

# わずかに足りない場合は、いつまでもゆっくり乾きつづける
slightly_dry = steady_state(1.02, 1.17, settings, max_days=120)
check(not slightly_dry.settled,
      "わずかに足りないだけなら、120日たっても落ち着かない（ゆっくり乾きつづける）",
      f"pF {slightly_dry.balance.end_pf:.2f} でまだ動いている")

for bad, label in [
    ({"irrigation_mm": -1.0}, "潅水量が負なら弾く"),
    ({"potential_transpiration_mm": -1.0}, "蒸散量が負なら弾く"),
    ({"max_days": 0}, "max_days 0 を弾く"),
]:
    kwargs = {"irrigation_mm": 5.0, "potential_transpiration_mm": 3.0}
    kwargs.update(bad)
    try:
        steady_state(**kwargs)
        check(False, label)
    except ValueError:
        check(True, label)

print()
print("=" * 78)
if failures:
    print(f"NG が {len(failures)} 件ある:")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("すべて通過")
print("=" * 78)
