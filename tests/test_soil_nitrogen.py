"""
根群域の無機態N収支の検証。

【何を確かめるか】
1. 収支が閉じるか（入れた分 − 出た分 = 残量の変化）
2. 濃度の計算が定義どおりか
3. 1日の刻み数が足りているか（順番の影響が消えているか）
4. 流亡が厳密解と合うか
5. 無機化の温度重みが正しいか
6. 入力の検査が効くか
7. ★冬に貯めたNは春まで残るか（根圏の洗い流しの半減期）
8. 需要が供給を超えたときの振る舞い
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import math                                                   # noqa: E402
from dataclasses import replace                               # noqa: E402

from core.soil_nitrogen import (                              # noqa: E402
    DailyNitrogenInput,
    SoilNitrogenSettings,
    concentration_mg_per_l,
    simulate,
    step_one_day,
    summarize,
    temperature_weights,
)

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


settings = SoilNitrogenSettings()
water_mm = settings.root_zone_water_mm()

# =============================================================================
print("=" * 78)
print("1. 収支が閉じるか")
print("=" * 78)
print(f"  根圏の水量（水槽の大きさ）{water_mm:.1f} L/m²"
      f"（深さ {settings.root_zone_depth_m} m × 含水率 "
      f"{settings.water_content:.3f}）")
print()

cases = [
    # (名前, 施肥, 需要, 潅水, 流亡, 無機化)
    ("12月（暗い・潅水少)", 0.164, 0.10, 1.5, 0.3, 0.05),
    ("5月（明るい・潅水多)", 0.140, 0.30, 9.0, 4.5, 0.15),
    ("無潅水の日", 0.0, 0.20, 0.0, 0.0, 0.10),
    ("基肥の日", 7.36, 0.0, 0.0, 0.0, 0.20),
]

worst = 0.0
for name, fertilizer, demand, irrigation, drainage, mineralization in cases:
    day = DailyNitrogenInput(
        date=name, fertilizer_n_g_per_m2=fertilizer,
        demand_n_g_per_m2=demand, irrigation_mm=irrigation,
        drainage_mm=drainage, mineralization_n_g_per_m2=mineralization)
    result = step_one_day(20.0, day, settings)
    worst = max(worst, abs(result.balance_residual_g_per_m2))
    print(f"  {name:<22} 残量 {result.start_n_g_per_m2:5.2f} → "
          f"{result.end_n_g_per_m2:5.2f}  吸収 {result.uptake_n_g_per_m2:.3f}  "
          f"流亡 {result.leaching_n_g_per_m2:.3f}  "
          f"濃度 {result.end_concentration_mg_per_l:5.1f} mg-N/L")
check(worst < 1e-9, "1日ぶんの収支が閉じる", f"最大残差 {worst:.2e} g/m²")

# =============================================================================
print()
print("=" * 78)
print("2. 濃度の計算")
print("=" * 78)
print("  濃度 [mg-N/L] = N量 [g/m²] ÷ 根圏の水量 [L/m²] × 1000")
print()
for nitrogen in (10.0, 20.0, 30.0):
    value = concentration_mg_per_l(nitrogen, water_mm)
    print(f"  N {nitrogen:5.1f} g/m² → {value:6.1f} mg-N/L")
expected = 20.0 / water_mm * 1000.0
check(abs(concentration_mg_per_l(20.0, water_mm) - expected) < 1e-12,
      "定義どおりの値になる", f"{expected:.3f} mg-N/L")
print()
print("  ※一般的な養液は 170〜220 mg-N/L。実績の液肥は12月 111・5月 18 mg-N/L。")

# =============================================================================
print()
print("=" * 78)
print("3. 1日の刻み数（順番の影響が消えているか）")
print("=" * 78)
print("  1日をまとめて計算すると「施肥・流亡・吸収のどれを先にするか」で")
print("  答えが変わる。刻めば差が縮むはず。")
print()

day = DailyNitrogenInput(
    date="5月", fertilizer_n_g_per_m2=0.14, demand_n_g_per_m2=0.30,
    irrigation_mm=9.0, drainage_mm=4.5, mineralization_n_g_per_m2=0.15)

reference = step_one_day(20.0, day, settings, substeps=2048)
print(f"{'刻み数':>8}{'残量':>12}{'2048刻みとの差':>18}")
for substeps in (1, 4, 24, 96, 2048):
    result = step_one_day(20.0, day, settings, substeps=substeps)
    gap = abs(result.end_n_g_per_m2 - reference.end_n_g_per_m2)
    relative = gap / reference.end_n_g_per_m2
    print(f"{substeps:>8}{result.end_n_g_per_m2:>12.5f}"
          f"{relative * 100:>16.4f}%")

result_24 = step_one_day(20.0, day, settings, substeps=24)
error_24 = abs(result_24.end_n_g_per_m2 - reference.end_n_g_per_m2) \
    / reference.end_n_g_per_m2
check(error_24 < 0.01, "24刻みで誤差1%未満",
      f"{error_24 * 100:.4f}%  → 既定の 24 で足りている")

# =============================================================================
print()
print("=" * 78)
print("4. 流亡が厳密解と合うか")
print("=" * 78)
print("  よく混ざった水槽を一定流量で流したときの残存率は exp(−D/W)。")
print("  吸収と入力をゼロにすれば、これだけが残る。")
print()

worst = 0.0
for drainage in (0.3, 1.0, 4.5, 9.0):
    day = DailyNitrogenInput(
        date=f"流亡{drainage}mm", fertilizer_n_g_per_m2=0.0,
        demand_n_g_per_m2=0.0, irrigation_mm=max(drainage, 0.0),
        drainage_mm=drainage, mineralization_n_g_per_m2=0.0)
    result = step_one_day(20.0, day, settings, substeps=24)
    exact = 20.0 * math.exp(-drainage / water_mm)
    worst = max(worst, abs(result.end_n_g_per_m2 - exact))
    print(f"  流亡 {drainage:4.1f} mm → 残量 {result.end_n_g_per_m2:8.5f} "
          f"（厳密解 {exact:8.5f}）  失った割合 "
          f"{(1 - result.end_n_g_per_m2 / 20.0) * 100:5.2f}%")
check(worst < 1e-9, "流亡が厳密解と一致する", f"最大差 {worst:.2e}")

# =============================================================================
print()
print("=" * 78)
print("5. 無機化の温度重み")
print("=" * 78)
print(f"  重み = Q10 ^ ((温度 − 基準) ÷ 10)、Q10 = "
      f"{settings.mineralization_q10}")
print()

# 作期の月別の日平均気温（ハウス内のおおよその水準）。
monthly_temp = {
    8: 27.0, 9: 24.0, 10: 19.0, 11: 15.0, 12: 13.0, 1: 13.0,
    2: 14.0, 3: 16.0, 4: 19.0, 5: 22.0, 6: 24.0,
}
months = list(monthly_temp)
weights = temperature_weights([monthly_temp[m] for m in months], settings)
total = sum(weights)
print(f"{'月':>4}{'気温':>8}{'重み':>9}{'無機化':>10}")
print(f"{'':>4}{'℃':>8}{'（相対）':>9}{'kg-N/10a':>10}")
for month, weight in zip(months, weights):
    share = weight / total
    print(f"{month:>4}{monthly_temp[month]:>8.1f}{weight / max(weights):>9.2f}"
          f"{share * settings.mineralization_total_g_per_m2:>10.2f}")
check(abs(total - 1.0) < 1e-12, "重みの合計が1", f"{total:.15f}")

warm = weights[months.index(8)]
cold = weights[months.index(12)]
check(warm > cold, "暖かい月のほうが重みが大きい",
      f"8月 {warm:.4f} > 12月 {cold:.4f}（{warm / cold:.2f} 倍）")

# 基準温度をずらしても結果が変わらないこと（正規化で消える）
shifted = temperature_weights(
    [monthly_temp[m] for m in months],
    replace(settings, mineralization_reference_temp_c=5.0))
worst = max(abs(a - b) for a, b in zip(weights, shifted))
check(worst < 1e-12, "基準温度を変えても重みは変わらない",
      f"最大差 {worst:.2e}")

print()
print("  ★無機化は夏と秋に多く、冬に少ない。運用（冬に濃く・春に薄く施肥）とは")
print("    逆向きなので、時期別の妥当性を見るうえでここが効く。")

# =============================================================================
print()
print("=" * 78)
print("6. 入力の検査が効くか")
print("=" * 78)


def expect_error(label: str, call) -> None:
    try:
        call()
    except ValueError as error:
        check(True, label, str(error)[:56])
    else:
        check(False, label, "エラーにならなかった")


expect_error("根群域の深さが0なら止まる",
             lambda: replace(settings, root_zone_depth_m=0.0))
expect_error("含水率が1以上なら止まる",
             lambda: replace(settings, water_content=1.2))
expect_error("Q10が1以下なら止まる",
             lambda: replace(settings, mineralization_q10=0.9))
expect_error("潅水の水のN濃度が負なら止まる",
             lambda: replace(settings, irrigation_water_n_mg_per_l=-1.0))
expect_error(
    "施肥Nが負なら止まる",
    lambda: DailyNitrogenInput(
        date="x", fertilizer_n_g_per_m2=-1.0, demand_n_g_per_m2=0.0,
        irrigation_mm=0.0, drainage_mm=0.0, mineralization_n_g_per_m2=0.0))
expect_error(
    "流亡が根圏の水量を超えたら止まる",
    lambda: step_one_day(
        1.0, DailyNitrogenInput(
            date="x", fertilizer_n_g_per_m2=0.0, demand_n_g_per_m2=0.0,
            irrigation_mm=500.0, drainage_mm=500.0,
            mineralization_n_g_per_m2=0.0), settings))

# ★流亡が潅水を超えるのは正常。前日に土へ溜まった水が、潅水しない日に
#   抜けていくからだ。水収支（core/water_balance.py）は貯留変化を持つので
#   潅水ゼロの日にも排水が出る。ここで弾くと実データが通らなくなる
#   （2026-10-01 に実際に起きた。2023-08-23 の潅水0・排水0.05 mm）。
ok_day = DailyNitrogenInput(
    date="潅水ゼロでも排水はある日", fertilizer_n_g_per_m2=0.0,
    demand_n_g_per_m2=0.1, irrigation_mm=0.0, drainage_mm=0.05,
    mineralization_n_g_per_m2=0.05)
result = step_one_day(20.0, ok_day, settings)
check(result.leaching_n_g_per_m2 > 0.0,
      "潅水ゼロの日に排水があっても通る（溜まった水が抜ける日）",
      f"流亡N {result.leaching_n_g_per_m2:.4f} g/m²")
expect_error(
    "前日の残量が負なら止まる",
    lambda: step_one_day(
        -1.0, DailyNitrogenInput(
            date="x", fertilizer_n_g_per_m2=0.0, demand_n_g_per_m2=0.0,
            irrigation_mm=0.0, drainage_mm=0.0,
            mineralization_n_g_per_m2=0.0), settings))
expect_error("集計する日がなければ止まる", lambda: summarize([]))

# =============================================================================
print()
print("=" * 78)
print("7. ★冬に貯めたNは春まで残るか（洗い流しの半減期）")
print("=" * 78)
print("  施肥の実績は「冬に濃く・春に薄く」なっている。これが")
print("  「冬に土へ貯めて春に引き出す」運用として成り立つのかを測る。")
print()
print("  入力も吸収もゼロにして、流亡だけで減る速さを見る。")
print("  半減期 = ln2 ÷ (流亡水量 ÷ 根圏の水量)")
print()

print(f"{'時期':<16}{'流亡 mm/日':>12}{'半減期 日':>11}{'100日後に残る':>14}")
half_lives = {}
for label, drainage in (("12月（暗い）", 0.3), ("2月", 1.0),
                        ("4月", 3.0), ("5月（明るい）", 4.5)):
    rate = drainage / water_mm
    half_life = math.log(2.0) / rate if rate > 0 else float("inf")
    remaining = math.exp(-rate * 100.0)
    half_lives[label] = half_life
    print(f"{label:<16}{drainage:>12.1f}{half_life:>11.0f}"
          f"{remaining * 100:>13.1f}%")

# 実際に回して確かめる
def flush_only(drainage_mm: float, days: int, start: float = 20.0) -> float:
    inputs = [
        DailyNitrogenInput(
            date=f"{i}日目", fertilizer_n_g_per_m2=0.0,
            demand_n_g_per_m2=0.0, irrigation_mm=drainage_mm,
            drainage_mm=drainage_mm, mineralization_n_g_per_m2=0.0)
        for i in range(days)
    ]
    return simulate(inputs, start, settings)[-1].end_n_g_per_m2


print()
spring_remaining = flush_only(4.5, 100) / 20.0
winter_remaining = flush_only(0.3, 100) / 20.0
print(f"  実際に100日回すと: 12月の条件で {winter_remaining * 100:.1f}% 残り、"
      f"5月の条件で {spring_remaining * 100:.1f}% 残る。")

check(winter_remaining > 0.7,
      "冬の条件ではNが土に残る",
      f"100日後に {winter_remaining * 100:.0f}%（半減期 "
      f"{half_lives['12月（暗い）']:.0f} 日）")
check(spring_remaining < 0.1,
      "★春の条件ではNが流れ去る",
      f"100日後に {spring_remaining * 100:.1f}%（半減期 "
      f"{half_lives['5月（明るい）']:.0f} 日）")

print()
print("  【読み方】")
print(f"  5月の半減期は {half_lives['5月（明るい）']:.0f} 日しかない。")
print("  つまり**冬に貯めたNは、春に入って1か月ほどで洗い流される**。")
print("  「冬に貯めて春に使う」という説明は成り立たない。")
print("  春の需要は、春に供給されるもの（施肥＋そのときの無機化）で")
print("  まかなわれているはずだ。無機化は暖かい春に増えるので（第5節）、")
print("  そちらが効いている可能性が高い。段5で日ごとのデータで確かめる。")

# =============================================================================
print()
print("=" * 78)
print("8. 需要が供給を超えたとき")
print("=" * 78)
print("  根圏にあるより多くは吸えない。足りない分は shortfall に出る。")
print()

inputs = [
    DailyNitrogenInput(
        date=f"{i}日目", fertilizer_n_g_per_m2=0.05,
        demand_n_g_per_m2=0.50, irrigation_mm=5.0,
        drainage_mm=2.5, mineralization_n_g_per_m2=0.05)
    for i in range(120)
]
results = simulate(inputs, 2.0, settings)
report = summarize(results)
print(report.describe())
check(abs(report.balance_residual) < 1e-9,
      "作期を通して収支が閉じる", f"残差 {report.balance_residual:.2e}")
check(report.short_days > 0, "N不足の日が検出される",
      f"{report.short_days} 日 / {report.days} 日")
check(all(r.end_n_g_per_m2 >= 0.0 for r in results),
      "残量が負にならない",
      f"最小 {min(r.end_n_g_per_m2 for r in results):.4f} g/m²")
check(report.uptake_ratio < 1.0, "吸えた割合が1未満になる",
      f"{report.uptake_ratio * 100:.1f}%")

# =============================================================================
print()
print("=" * 78)
if failures:
    print(f"NG が {len(failures)} 件ある:")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("すべて通過")
print("=" * 78)
