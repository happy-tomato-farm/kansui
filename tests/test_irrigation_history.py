"""
過去の潅水実績の読み込みを検証する。

【これは何のためのデータか】
モデルが出した潅水量を、自分の過去の実績と並べて確かめるためのもの。
★モデルの計算にはいっさい使っていない。答え合わせ用の参考値。

【いちばん大事な性質】
同じ時期でも、明るい日と暗い日では潅水量がまるで違う。
だから **明るさをそろえて比べられること** がこのモジュールの肝。
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import datetime as dt

from core.advisor import normal_radiation_mj
from core.irrigation_history import HISTORY_PATH, describe_source, lookup
from core.solar import clear_sky_sensor_mj, cloudy_sensor_mj

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


# =============================================================================
print("=" * 78)
print("1. データが読めるか")
print("=" * 78)

check(HISTORY_PATH.exists(), "潅水実績のファイルがある", str(HISTORY_PATH.name))
print(f"  出どころ: {describe_source()}")

# =============================================================================
print()
print("=" * 78)
print("2. 明るさをそろえて引けるか（このモジュールの肝）")
print("=" * 78)
print("  同じ日付でも、快晴の日と曇天の日では潅水量が違うはず。")
print()
print("   日付     空     渡した日射  日数  実績の日射  ふだん   最大")

brightness_samples = {}
for month, day in [(10, 15), (12, 15), (2, 15), (4, 15), (5, 15)]:
    year = 2025 if month >= 9 else 2026
    date = dt.date(year, month, day)
    day_of_year = date.timetuple().tm_yday
    per_sky = {}
    for label, radiation in [
        ("快晴", clear_sky_sensor_mj(date)),
        ("平年", normal_radiation_mj(day_of_year)),
        ("曇天", cloudy_sensor_mj(date)),
    ]:
        entry = lookup(date, "中央", sensor_radiation_mj=radiation)
        if entry is None:
            print(f"   {month:2d}/{day}  {label}  {radiation:8.1f}   （記録が足りない）")
            continue
        per_sky[label] = entry
        print(f"   {month:2d}/{day}  {label}  {radiation:8.1f} {entry.days:5d}"
              f" {entry.median_radiation_mj:10.1f} {entry.median_l_per_m2:7.2f}"
              f" {entry.max_l_per_m2:6.2f}")
    brightness_samples[month] = per_sky
    print()

# 快晴の日のほうが、平年なみの日より多く出しているはず。
ordered = 0
compared = 0
for month, per_sky in brightness_samples.items():
    if "快晴" in per_sky and "平年" in per_sky:
        compared += 1
        if per_sky["快晴"].median_l_per_m2 > per_sky["平年"].median_l_per_m2:
            ordered += 1
check(compared > 0 and ordered == compared,
      "どの時期でも、快晴の日のほうが平年なみの日より多く潅水している",
      f"{ordered} / {compared} 時期")

# 選ばれた日の日射が、渡した日射に近いこと（±35%以内）。
worst_label, worst_gap = None, 0.0
for month, per_sky in brightness_samples.items():
    for label, entry in per_sky.items():
        if not entry.matched_brightness:
            continue
        date = dt.date(2025 if month >= 9 else 2026, month, 15)
        asked = {"快晴": clear_sky_sensor_mj(date),
                 "平年": normal_radiation_mj(date.timetuple().tm_yday),
                 "曇天": cloudy_sensor_mj(date)}[label]
        gap = abs(entry.median_radiation_mj - asked) / asked
        if gap > worst_gap:
            worst_label, worst_gap = f"{month}月 {label}", gap
check(worst_gap < 0.35, "選ばれた日の日射が、渡した日射に近い（±35%以内）",
      f"いちばん離れたのは {worst_label} の {worst_gap * 100:.0f}%")

# =============================================================================
print("=" * 78)
print("3. 季節の形が実態に合っているか")
print("=" * 78)
print("  快晴の日どうしで比べる。春に多く、冬に少ないはず。")
print()
print("   日付        ふだん   最大    平均   日数   作期")

season_samples = {}
for month, day in [(10, 15), (12, 15), (2, 15), (4, 15), (5, 15), (6, 15)]:
    year = 2025 if month >= 9 else 2026
    date = dt.date(year, month, day)
    entry = lookup(date, "中央",
                   sensor_radiation_mj=clear_sky_sensor_mj(date))
    if entry is None:
        print(f"   {date}  （記録が足りない）")
        continue
    season_samples[month] = entry
    print(f"   {date}  {entry.median_l_per_m2:6.2f}  {entry.max_l_per_m2:6.2f}"
          f"  {entry.mean_l_per_m2:6.2f}  {entry.days:5d}   {entry.season_range}")

check(len(season_samples) >= 5, "主な時期の記録がそろっている",
      f"{len(season_samples)} 時期")

if 12 in season_samples and 5 in season_samples:
    winter = season_samples[12].median_l_per_m2
    spring = season_samples[5].median_l_per_m2
    check(spring > winter * 3,
          "5月の潅水量は12月の3倍を超える（季節の形が出ている）",
          f"12月 {winter:.2f} → 5月 {spring:.2f} L/m²")

for month, entry in season_samples.items():
    if entry.max_l_per_m2 < entry.median_l_per_m2:
        check(False, f"{month}月: 最大が中央値を下回っている")
        break
else:
    check(True, "どの時期でも 最大 ≧ 中央値 になっている")

# =============================================================================
print()
print("=" * 78)
print("4. 値がありえない大きさになっていないか")
print("=" * 78)
print("  ★以前、古い年の「潅水時間(分)」を L/m² として読んでいたことがある。")
print("    5月の最大が 30 L/m² と出ていたのがそれ（実際は30分）。")

worst_month, worst_max = None, 0.0
for month in range(1, 13):
    year = 2025 if month >= 9 else 2026
    date = dt.date(year, month, 15)
    entry = lookup(date, "中央")
    if entry is None:
        continue
    if entry.max_l_per_m2 > worst_max:
        worst_month, worst_max = month, entry.max_l_per_m2
check(worst_max < 20.0, "どの時期でも潅水量が 20 L/m² を超えない",
      f"最大は {worst_month}月 の {worst_max:.2f} L/m²")

# =============================================================================
print()
print("=" * 78)
print("5. 東ハウス・記録のない時期・誤った指定")
print("=" * 78)

east = lookup(dt.date(2026, 5, 15), "東",
              sensor_radiation_mj=clear_sky_sensor_mj(dt.date(2026, 5, 15)))
central = lookup(dt.date(2026, 5, 15), "中央",
                 sensor_radiation_mj=clear_sky_sensor_mj(dt.date(2026, 5, 15)))
check(east is not None, "東ハウスの記録がある")
if east and central:
    print(f"  5月15日（快晴）  中央 {central.median_l_per_m2:.2f} / "
          f"東 {east.median_l_per_m2:.2f} L/m²")
    check(east.median_l_per_m2 > 0, "東の値が正")

summer = lookup(dt.date(2026, 8, 15), "中央")
check(summer is None, "作期外（8月）は None を返す")

try:
    lookup(dt.date(2026, 5, 15), "西")
    check(False, "記録にないハウス名はエラーになる")
except KeyError as error:
    check("記録がない" in str(error), "記録にないハウス名はエラーになる")

try:
    lookup(dt.date(2026, 5, 15), "中央", sensor_radiation_mj=-1.0)
    check(False, "日射が負ならエラーになる")
except ValueError as error:
    check("日射" in str(error), "日射が負ならエラーになる")

# 日射を渡さなければ、天気を問わない値が返る
plain = lookup(dt.date(2026, 5, 15), "中央")
check(plain is not None and not plain.matched_brightness,
      "日射を渡さなければ天気を問わない値になる",
      f"{plain.days} 日・{plain.basis}" if plain else "")

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
