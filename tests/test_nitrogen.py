"""
窒素（N）需要の計算の検証。

【何を確かめるか】
1. 器官別の内訳の合計が、全身平均N濃度から出した値と一致するか
2. 単位の一致（1 g-N/m² = 1 kg-N/10a）が崩れていないか
3. 「果実1トンあたりN」が文献の範囲に入るか（独立した答え合わせ）
4. 2つの入口（糖から／収量から）が、較正係数を通したときに一致するか
5. 地力窒素の逆算が、施肥量の計算と往復して元に戻るか
6. 入力の検査が効くか
7. どの前提がいちばんN需要を動かすか（感度）
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import datetime as dt                                            # noqa: E402
from dataclasses import replace                                   # noqa: E402

from core.advisor import normal_radiation_mj                      # noqa: E402
from core.solar import clear_sky_tau, sensor_basis_radiation_mj   # noqa: E402

from config import (                                              # noqa: E402
    FERTILIZER_N_EFFICIENCY,
    N_PER_TONNE_FRUIT_LITERATURE_RANGE,
    ORGAN_N_CONTENT_SENSITIVITY,
    YIELD_CONVERSION,
    YIELD_CONVERSION_SENSITIVITY,
)
from core.nitrogen import (                                       # noqa: E402
    ORGANS,
    ORGAN_LABELS,
    NitrogenSettings,
    advise_fertilizer_n,
    advise_fertilizer_n_week,
    demand_from_dry_matter,
    demand_from_fresh_yield,
    demand_from_sugar,
    required_fertilizer_n,
    season_light_capture,
    solve_soil_n_supply,
)

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


settings = NitrogenSettings()

# =============================================================================
print("=" * 78)
print("1. 器官別の内訳が全身平均と合うか")
print("=" * 78)
print(settings.describe())
print()

# 分配率の合計は1でなければならない。
partition_total = sum(settings.dry_matter_partition.values())
check(abs(partition_total - 1.0) < 1e-12,
      "乾物の分配率の合計が1", f"{partition_total:.15f}")

worst = 0.0
for dry_matter in (500.0, 1000.0, 2850.0, 4000.0):
    demand = demand_from_dry_matter(dry_matter, settings)
    expected = dry_matter * settings.whole_plant_n_content
    worst = max(worst, abs(demand.uptake_g_per_m2 - expected))
    print(f"  全乾物 {dry_matter:7.0f} g/m² → 吸収N "
          f"{demand.uptake_g_per_m2:7.2f} g/m²（内訳の和）/ "
          f"{expected:7.2f}（平均N濃度から）")
check(worst < 1e-9, "内訳の和 = 全乾物 × 全身平均N濃度", f"最大差 {worst:.2e}")

# 果実N＋残渣N＝吸収N
demand = demand_from_dry_matter(2850.0, settings)
split = demand.removed_g_per_m2 + demand.residue_g_per_m2
check(abs(split - demand.uptake_g_per_m2) < 1e-9,
      "持ち出しN ＋ 残渣N ＝ 吸収N", f"差 {abs(split - demand.uptake_g_per_m2):.2e}")

# =============================================================================
print()
print("=" * 78)
print("2. 単位の一致（1 g-N/m² = 1 kg-N/10a）")
print("=" * 78)
print("  1 g/m² × 10000 m²/ha = 10000 g/ha = 10 kg/ha = 1 kg/10a")
check(abs(demand.uptake_g_per_m2 - demand.uptake_kg_per_10a) < 1e-12,
      "g/m² と kg/10a が同じ値",
      f"{demand.uptake_g_per_m2:.4f} / {demand.uptake_kg_per_10a:.4f}")

# =============================================================================
print()
print("=" * 78)
print("3. 果実1トンあたりN（独立した答え合わせ）")
print("=" * 78)

low, high = N_PER_TONNE_FRUIT_LITERATURE_RANGE
value = settings.n_per_tonne_fruit_kg
print(f"  いまの前提: {value:.2f} kg-N/t  （文献 {low}〜{high}）")

# 定義どおりか。収量から直接出した吸収Nを収量で割っても同じ値になるはず。
#
# ★もう1つ便利な単位の一致: 1 kg/m² = 1 t/10a
#   （1 kg/m² × 1000 m²/10a = 1000 kg/10a = 1 t/10a）
#   なので、収量 [kg/m²] の数字をそのまま [t/10a] として割ればよい。
for harvest in (20.0, 30.0, 31.4, 40.0):
    uptake = demand_from_fresh_yield(harvest, settings).uptake_kg_per_10a
    per_tonne = uptake / harvest
    print(f"  収量 {harvest:5.1f} kg/m²（= {harvest:.1f} t/10a）→ 吸収N "
          f"{uptake:6.2f} kg/10a → {per_tonne:.3f} kg-N/t")
    if abs(per_tonne - value) > 1e-9:
        check(False, "1トンあたりNが収量によらず一定", f"{per_tonne} ≠ {value}")
        break
else:
    check(True, "1トンあたりNが収量によらず一定", f"{value:.3f} kg-N/t")

check(low <= value <= high, "文献の範囲に入っている",
      f"{value:.2f} は {low}〜{high} の"
      + ("下端寄り" if value < low + (high - low) * 0.25 else "中ほど"))

# =============================================================================
print()
print("=" * 78)
print("4. 2つの入口（糖から／収量から）が一致するか")
print("=" * 78)
print("  糖から出した乾物と、収量から遡った乾物が同じなら、同じN需要になる。")
print("  ★較正係数は、糖から入る経路だけに掛かる。")
print()

harvest = 31.4
from_yield = demand_from_fresh_yield(harvest, settings)

# この収量を出すのに必要だった糖（第8節の逆算と同じ道すじ）
needed_sugar_g = (from_yield.total_dry_matter_g_per_m2
                  / settings.sugar_to_dry_matter)
from_sugar = demand_from_sugar(needed_sugar_g, settings, calibration_factor=1.0)
print(f"  収量 {harvest} kg/m² → 吸収N {from_yield.uptake_g_per_m2:.3f} g/m²")
print(f"  必要だった糖 {needed_sugar_g:.1f} g/m² → 吸収N "
      f"{from_sugar.uptake_g_per_m2:.3f} g/m²")
check(abs(from_yield.uptake_g_per_m2 - from_sugar.uptake_g_per_m2) < 1e-9,
      "2つの入口が一致する",
      f"差 {abs(from_yield.uptake_g_per_m2 - from_sugar.uptake_g_per_m2):.2e}")

# 較正係数は比例して効く
for factor in (1.0, 2.0, 3.0):
    scaled = demand_from_sugar(1000.0, settings, calibration_factor=factor)
    plain = demand_from_sugar(1000.0, settings, calibration_factor=1.0)
    ratio = scaled.uptake_g_per_m2 / plain.uptake_g_per_m2
    print(f"  較正係数 {factor:.1f} → N需要は素の {ratio:.3f} 倍")
    if abs(ratio - factor) > 1e-9:
        check(False, "較正係数が比例して効く", f"{ratio} ≠ {factor}")
        break
else:
    check(True, "較正係数が比例して効く")

print()
print("  ★素の光合成（較正係数1.0）は実収量の約1/3しか出ないので、")
print("    N需要も1/3で出る。絶対値には収量からの経路を使うこと。")

# =============================================================================
print()
print("=" * 78)
print("5. 地力窒素の逆算が往復するか")
print("=" * 78)

worst = 0.0
for uptake in (50.0, 62.0, 65.0, 70.0):
    for efficiency in (0.6, 0.7, 0.8, 1.0):
        for fertilizer in (30.0, 50.0, 70.0):
            soil = solve_soil_n_supply(uptake, fertilizer, efficiency)
            plan = required_fertilizer_n(uptake, soil, efficiency)
            # 土の供給が需要を超えていなければ、施肥量は元に戻るはず
            if soil < uptake:
                worst = max(worst, abs(plan.fertilizer_kg_per_10a - fertilizer))
print(f"  吸収N・利用率・施肥量を振って往復させた。最大差 {worst:.2e} kg/10a")
check(worst < 1e-9, "逆算 → 施肥量 が元に戻る", f"最大差 {worst:.2e}")

# 土の供給が需要を超えたら施肥量は0（負を返さない）
plan = required_fertilizer_n(50.0, 60.0, 0.7)
check(plan.fertilizer_kg_per_10a == 0.0,
      "土の供給が需要を超えたら施肥量は0",
      f"足りない分 {plan.shortfall_kg_per_10a:+.1f} → 施肥 "
      f"{plan.fertilizer_kg_per_10a:.1f} kg/10a")

# =============================================================================
print()
print("=" * 78)
print("6. 入力の検査が効くか")
print("=" * 78)


def expect_error(label: str, call) -> None:
    try:
        call()
    except (ValueError, KeyError) as error:
        check(True, label, str(error)[:58])
    else:
        check(False, label, "エラーにならなかった")


expect_error("N濃度を % で渡したら止まる（4 ではなく 0.04）",
             lambda: replace(settings, leaf_n=4.0))
expect_error("栄養器官の内訳の合計が1でないと止まる",
             lambda: replace(settings, leaf_share=0.9))
expect_error("果実分配率が1を超えたら止まる",
             lambda: replace(settings, fruit_allocation=1.5))
expect_error("全乾物が負なら止まる",
             lambda: demand_from_dry_matter(-1.0, settings))
expect_error("収量が負なら止まる",
             lambda: demand_from_fresh_yield(-1.0, settings))
expect_error("較正係数が0以下なら止まる",
             lambda: demand_from_sugar(1000.0, settings, calibration_factor=0.0))
expect_error("利用率が1を超えたら止まる",
             lambda: required_fertilizer_n(60.0, 20.0, 1.2))
expect_error("利用率が0なら止まる",
             lambda: solve_soil_n_supply(60.0, 50.0, 0.0))

# =============================================================================
print()
print("=" * 78)
print("7. どの前提がいちばんN需要を動かすか（感度）")
print("=" * 78)
print("  収量 31.4 kg/m²（2025作期・中央）での吸収N量。")
print("  ★幅が大きいものほど、実測する価値が高い。")
print()

base_uptake = demand_from_fresh_yield(31.4, settings).uptake_kg_per_10a
print(f"{'動かすもの':<16}{'値':>9}{'吸収N':>10}{'既定比':>9}")
print(f"{'（既定）':<16}{'—':>9}{base_uptake:>10.1f}{1.0:>9.3f}")

ranges: list[tuple[str, tuple, str]] = [
    ("果実乾物率", YIELD_CONVERSION_SENSITIVITY["FRUIT_DRY_MATTER_CONTENT"],
     "fruit_dry_matter_content"),
    ("果実分配率", YIELD_CONVERSION_SENSITIVITY["FRUIT_ALLOCATION"],
     "fruit_allocation"),
    ("果実のN濃度", ORGAN_N_CONTENT_SENSITIVITY["FRUIT"], "fruit_n"),
    ("葉のN濃度", ORGAN_N_CONTENT_SENSITIVITY["LEAF"], "leaf_n"),
    ("茎のN濃度", ORGAN_N_CONTENT_SENSITIVITY["STEM"], "stem_n"),
    ("根のN濃度", ORGAN_N_CONTENT_SENSITIVITY["ROOT"], "root_n"),
]

spread: dict[str, float] = {}
for label, values, field in ranges:
    print(f"{'':<16}{'':>9}{'':>10}{'':>9}")
    seen = []
    for value in values:
        custom = replace(settings, **{field: value})
        uptake = demand_from_fresh_yield(31.4, custom).uptake_kg_per_10a
        seen.append(uptake)
        shown = (f"{value * 100:.1f}%" if value < 0.1 else f"{value:.2f}")
        print(f"{label:<16}{shown:>9}{uptake:>10.1f}"
              f"{uptake / base_uptake:>9.3f}")
    spread[label] = max(seen) / min(seen)

print()
print("  効き方（幅の最大÷最小）の大きい順")
for label, ratio in sorted(spread.items(), key=lambda x: -x[1]):
    print(f"    {label:<14}{ratio:.3f} 倍")

biggest = max(spread, key=spread.get)
check(spread[biggest] > 1.1,
      "いちばん効く前提がはっきりしている", f"{biggest}（{spread[biggest]:.2f} 倍）")

# 器官別N濃度より、乾物側（乾物率・分配率）のほうが効くこと。
# これが成り立つなら、優先して測るのは乾物のほうだ。
dry_side = max(spread["果実乾物率"], spread["果実分配率"])
n_side = max(spread["果実のN濃度"], spread["葉のN濃度"],
             spread["茎のN濃度"], spread["根のN濃度"])
check(dry_side > n_side,
      "乾物側の前提のほうが、N濃度の前提より効く",
      f"乾物側 {dry_side:.3f} 倍 > N濃度側 {n_side:.3f} 倍 "
      f"→ 先に測るべきは乾物率と分配率")

# =============================================================================
print()
print("=" * 78)
print("8. 朝に使う『1日のN量』（段6）")
print("=" * 78)
print("  目標収量 → 作期の施肥N → その日の受光量の取り分 で割る。")
print("  ★配分の重みは受光量（日射 × 受光率）。糖は使わない（第4節の★）。")
print()

season = season_light_capture(2025)
print(f"  作期2025の平年受光量の合計 {season:.0f} MJ/m²")

example = advise_fertilizer_n(dt.date(2026, 1, 15), 8.0, 31.0)
print(f"  作期の吸収N需要   {example.season_uptake_kg_per_10a:6.1f} kg-N/10a")
print(f"  作期の施肥N       {example.season_fertilizer_kg_per_10a:6.1f} kg-N/10a"
      f"（実績 47〜51）")
check(47.0 <= example.season_fertilizer_kg_per_10a <= 51.0,
      "作期の施肥Nが実績の幅に入る",
      f"{example.season_fertilizer_kg_per_10a:.1f} kg-N/10a")

print()
print(f"{'日付':>12}{'日射':>7}{'取り分':>9}{'1日のN量':>11}")
print(f"{'':>12}{'MJ':>7}{'%':>9}{'kg-N/10a':>11}")
shares = 0.0
profile = {}
for month, day in ((10, 15), (12, 15), (2, 15), (4, 15), (5, 15)):
    year = 2025 if month >= 8 else 2026
    date = dt.date(year, month, day)
    radiation = normal_radiation_mj(date.timetuple().tm_yday)
    result = advise_fertilizer_n(date, radiation, 31.0)
    profile[month] = result.daily_n_kg_per_10a
    print(f"{str(date):>12}{radiation:>7.1f}{result.share * 100:>9.3f}"
          f"{result.daily_n_kg_per_10a:>11.3f}")

check(profile[5] > profile[12],
      "★N需要ベースでは春のほうが多くなる（運用と逆向き）",
      f"5月 {profile[5]:.3f} > 12月 {profile[12]:.3f} kg-N/10a。"
      f"運用は 12月 0.164 > 5月 0.140 なので逆。"
      f"★量はこちら（N需要）で決め、ECは帯で挟むだけにした"
      f"（config 第11-8節。以前は逆にしていたが撤回）")

# 取り分を作期ぶん足すと1になること（平年日射のとき）
date = dt.date(2025, 8, 1)
while date <= dt.date(2026, 7, 31):
    shares += advise_fertilizer_n(
        date, normal_radiation_mj(date.timetuple().tm_yday), 31.0).share
    date += dt.timedelta(days=1)
print()
check(abs(shares - 1.0) < 1e-9,
      "平年日射なら、作期の取り分の合計が1になる", f"{shares:.12f}")

# 日射が明るいほど取り分が増える
bright = advise_fertilizer_n(dt.date(2026, 5, 15), 20.0, 31.0)
dark = advise_fertilizer_n(dt.date(2026, 5, 15), 5.0, 31.0)
check(bright.daily_n_kg_per_10a > dark.daily_n_kg_per_10a,
      "明るい日は取り分が増える",
      f"20 MJ → {bright.daily_n_kg_per_10a:.3f} / "
      f"5 MJ → {dark.daily_n_kg_per_10a:.3f} kg-N/10a")

# 定植前はゼロ
before = advise_fertilizer_n(dt.date(2025, 8, 20), 15.0, 31.0)
check(before.before_planting and before.daily_n_kg_per_10a == 0.0,
      "定植前はゼロになる", f"{before.date} → {before.daily_n_kg_per_10a:.3f}")

# 目標収量に比例すること（作期の施肥Nは収量に対して線形ではない。
# 地力窒素を引いてから割るので切片がある）
low = advise_fertilizer_n(dt.date(2026, 1, 15), 8.0, 25.0)
high = advise_fertilizer_n(dt.date(2026, 1, 15), 8.0, 35.0)
print()
print(f"  目標収量 25 kg/m² → 作期の施肥N "
      f"{low.season_fertilizer_kg_per_10a:.1f} ／ "
      f"35 kg/m² → {high.season_fertilizer_kg_per_10a:.1f} kg-N/10a")
check(high.season_fertilizer_kg_per_10a > low.season_fertilizer_kg_per_10a,
      "目標収量を上げれば施肥Nも増える")
check(abs(high.season_fertilizer_kg_per_10a
          - low.season_fertilizer_kg_per_10a
          - (demand_from_fresh_yield(35.0).uptake_kg_per_10a
             - demand_from_fresh_yield(25.0).uptake_kg_per_10a)
          / FERTILIZER_N_EFFICIENCY) < 1e-9,
      "増え方は needs の差 ÷ 利用率（地力窒素は切片として効く）")

expect_error("日射が負なら止まる",
             lambda: advise_fertilizer_n(dt.date(2026, 1, 15), -1.0, 31.0))

# =============================================================================
print()
print("=" * 78)
print("9. 1週間で使うN量（レシピの「1週間で使いたいN量」）")
print("=" * 78)
print("  ★レシピ側は nday = nweek ÷ 7 と割るだけなので、")
print("    こちらも快晴7日ぶんの合計を渡す。向こうの割り算と合わせる。")
print()

# --- 9-1. 日別の合計と一致すること（式を二重に書いていないことの確認）---
start = dt.date(2026, 1, 15)
week_radiations = [8.0, 9.0, 7.5, 10.0, 6.0, 8.5, 9.5]
week = advise_fertilizer_n_week(start, week_radiations, 31.0)
by_day = sum(
    advise_fertilizer_n(start + dt.timedelta(days=i), mj, 31.0).daily_n_kg_per_10a
    for i, mj in enumerate(week_radiations)
)
check(abs(week.week_n_kg_per_10a - by_day) < 1e-12,
      "週の合計は日別の合計と一致する",
      f"週 {week.week_n_kg_per_10a:.6f} / 日別の足し上げ {by_day:.6f} kg-N/10a")
check(week.days == 7 and week.growing_days == 7,
      "1月の7日間はすべて定植後",
      f"days {week.days} / growing_days {week.growing_days}")
check(week.end_date == start + dt.timedelta(days=6),
      "最終日は初日＋6日", f"{week.start_date} 〜 {week.end_date}")

# --- 9-2. 7で割った値がレシピ側の1日N量になること ---
check(abs(week.daily_n_kg_per_10a - week.week_n_kg_per_10a / 7.0) < 1e-12,
      "1日N量は週N量 ÷ 7（レシピと同じ割り算）",
      f"{week.week_n_kg_per_10a:.4f} ÷ 7 = {week.daily_n_kg_per_10a:.5f}")

# --- 9-3. 日射が明るいほど週の目標は増えること ---
bright = advise_fertilizer_n_week(start, [12.0] * 7, 31.0)
dark = advise_fertilizer_n_week(start, [4.0] * 7, 31.0)
check(bright.week_n_kg_per_10a > dark.week_n_kg_per_10a,
      "明るい週のほうが目標は大きい",
      f"日射12 → {bright.week_n_kg_per_10a:.3f} / "
      f"日射4 → {dark.week_n_kg_per_10a:.3f} kg-N/10a")
check(abs(bright.week_n_kg_per_10a / dark.week_n_kg_per_10a - 3.0) < 1e-9,
      "日射が3倍なら目標も3倍（受光率は同じ週なので比例する）",
      f"比 {bright.week_n_kg_per_10a / dark.week_n_kg_per_10a:.6f}")

# --- 9-4. 定植前をまたぐ週は、定植前の日が入らないこと ---
# 作期は8月始まりだが定植は ISO 週38（9月中旬）
across = advise_fertilizer_n_week(dt.date(2025, 9, 12), [10.0] * 7, 31.0)
check(0 < across.growing_days < across.days,
      "定植をまたぐ週は一部の日だけ数える",
      f"{across.describe()}")
before_only = advise_fertilizer_n_week(dt.date(2025, 8, 10), [10.0] * 7, 31.0)
check(before_only.week_n_kg_per_10a == 0.0 and before_only.growing_days == 0,
      "全日が定植前なら週の目標はゼロ",
      before_only.describe())

# --- 9-5. 作期を通して足すと、作期の施肥Nに戻ること ---
# 平年日射で1年ぶんを週に切って足す。端数が出ないように日で回す。
season_total = 0.0
date_cursor = dt.date(2025, 8, 1)
while date_cursor <= dt.date(2026, 7, 31):
    season_total += advise_fertilizer_n(
        date_cursor, normal_radiation_mj(date_cursor.timetuple().tm_yday), 31.0
    ).daily_n_kg_per_10a
    date_cursor += dt.timedelta(days=1)
check(abs(season_total - week.season_fertilizer_kg_per_10a) < 0.01,
      "平年日射で作期を通して足すと作期の施肥Nに戻る",
      f"足し上げ {season_total:.3f} / 作期の施肥N "
      f"{week.season_fertilizer_kg_per_10a:.3f} kg-N/10a")

# --- 9-6. 快晴基準の週の目標が、平年基準より大きいこと ---
# 快晴が7日つづくのは平年より明るい想定なので、目標は上に出る
normal_week = advise_fertilizer_n_week(
    start,
    [normal_radiation_mj((start + dt.timedelta(days=i)).timetuple().tm_yday)
     for i in range(7)],
    31.0)
clear_week = advise_fertilizer_n_week(
    start,
    [sensor_basis_radiation_mj(
        start + dt.timedelta(days=i),
        tau=clear_sky_tau((start + dt.timedelta(days=i)).timetuple().tm_yday))
     for i in range(7)],
    31.0)
check(clear_week.week_n_kg_per_10a > normal_week.week_n_kg_per_10a,
      "快晴7日の目標は平年7日より大きい",
      f"快晴 {clear_week.week_n_kg_per_10a:.3f} / "
      f"平年 {normal_week.week_n_kg_per_10a:.3f} kg-N/10a")

# 月ごとに快晴7日の目標を並べる（管理の目安として見る表）
print(f"{'初日':>12}{'快晴日射の平均':>16}{'週の目標':>11}{'1日ぶん':>10}")
print(f"{'':>12}{'MJ/m²':>16}{'kg-N/10a':>11}{'kg-N/10a':>10}")
for month, year in ((10, 2025), (11, 2025), (12, 2025),
                    (1, 2026), (2, 2026), (3, 2026),
                    (4, 2026), (5, 2026), (6, 2026)):
    first_day = dt.date(year, month, 1)
    mjs = [sensor_basis_radiation_mj(
        first_day + dt.timedelta(days=i),
        tau=clear_sky_tau((first_day + dt.timedelta(days=i)).timetuple().tm_yday))
        for i in range(7)]
    advice_week = advise_fertilizer_n_week(first_day, mjs, 31.0)
    print(f"{first_day.isoformat():>12}{sum(mjs) / 7:>16.1f}"
          f"{advice_week.week_n_kg_per_10a:>11.3f}"
          f"{advice_week.daily_n_kg_per_10a:>10.4f}")

# --- 9-7. 異常な入力は止まること ---
expect_error("日射の並びが空なら止まる",
             lambda: advise_fertilizer_n_week(start, [], 31.0))
expect_error("並びの中に負の日射があれば止まる",
             lambda: advise_fertilizer_n_week(start, [8.0, -1.0, 8.0], 31.0))

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
