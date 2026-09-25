"""
「今日の潅水量のめやす」の検証。

【何を確かめるか】
1. 日射の日変化パターンが、日積算を正しく保っているか
2. 日射予測を上げると蒸散も潅水も増えるか（単調性）
3. 日積算が同じなら、日変化の形が多少違っても蒸散量は大きく変わらないか
4. 実データの月別平均と突き合わせて、ずれが許せる範囲か
5. 入力の検査が効くか
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dataclasses import replace

from config import FILM_DEGRADATION_FACTOR, HOUSE_SPECS, RECIPE_RADIATION_COEF
from core.advisor import (
    STEPS_PER_DAY,
    DayForecast,
    _diurnal_shape,
    advise,
    forecast_from_date,
    forecast_from_day_of_year,
    forecast_from_month,
    normal_radiation_mj,
)

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


# =============================================================================
print("=" * 78)
print("1. 日射の日変化パターン")
print("=" * 78)

forecast = forecast_from_month(12.0, 4)
shape = _diurnal_shape(forecast)
print(f"  ステップ数 {len(shape)}（1日288ステップ＝5分刻み）")
print(f"  合計 {sum(shape):.10f}")
check(len(shape) == STEPS_PER_DAY, "288ステップある")
check(abs(sum(shape) - 1.0) < 1e-12, "配分の合計がちょうど1", f"{sum(shape):.2e}")
check(all(w >= 0.0 for w in shape), "負の配分がない")

sunrise, sunset = forecast.daylight_hours
night = [w for i, w in enumerate(shape)
         if not sunrise <= (i + 0.5) * 24.0 / STEPS_PER_DAY < sunset]
check(all(w == 0.0 for w in night), "夜間の配分はゼロ", f"{len(night)}ステップ")

peak_index = shape.index(max(shape))
peak_hour = (peak_index + 0.5) * 24.0 / STEPS_PER_DAY
print(f"  山の位置 {peak_hour:.1f}時（日中 {sunrise:.0f}〜{sunset:.0f}時 の中央）")
check(abs(peak_hour - (sunrise + sunset) / 2.0) < 0.2, "山が日中の真ん中に来る")

# =============================================================================
print()
print("=" * 78)
print("2. 日射予測に対する単調性")
print("=" * 78)
print("  日射    蒸散     10MJあたり  潅水めやす  流亡率")

previous_transpiration = -1.0
previous_irrigation = -1.0
monotonic = True
for mj in (2.0, 5.0, 8.0, 11.0, 14.0, 18.0):
    advice = advise(forecast_from_month(mj, 4))
    print(f"  {mj:5.1f}  {advice.transpiration_l_per_m2:6.2f}  "
          f"{advice.water_per_10mj_l_per_m2:9.2f}  "
          f"{advice.recommended_irrigation_l_per_m2:9.2f}  "
          f"{advice.drainage_fraction * 100:6.1f}%")
    if advice.transpiration_l_per_m2 <= previous_transpiration:
        monotonic = False
    if advice.recommended_irrigation_l_per_m2 <= previous_irrigation:
        monotonic = False
    previous_transpiration = advice.transpiration_l_per_m2
    previous_irrigation = advice.recommended_irrigation_l_per_m2

check(monotonic, "日射が増えると蒸散も潅水も必ず増える")

zero = advise(forecast_from_month(0.0, 4))
check(zero.transpiration_l_per_m2 >= 0.0, "日射ゼロでも蒸散が負にならない",
      f"{zero.transpiration_l_per_m2:.4f} L/m²")

# =============================================================================
print()
print("=" * 78)
print("3. 日射の読みかえが既存アプリと同じか")
print("=" * 78)

advice = advise(forecast_from_month(10.0, 4))
expected_outside = 10.0 / RECIPE_RADIATION_COEF
expected_effective = expected_outside * FILM_DEGRADATION_FACTOR
print(f"  センサー 10.0 → {advice.outside_radiation_mj:.2f} "
      f"→ 実質 {advice.effective_radiation_mj:.2f} MJ/m²")
print(f"  既存アプリの式: 10.0 ÷ {RECIPE_RADIATION_COEF} × {FILM_DEGRADATION_FACTOR} "
      f"= {expected_effective:.2f}")
check(abs(advice.outside_radiation_mj - expected_outside) < 1e-9,
      "分母の元が「センサー値 ÷ レシピ換算係数」になっている")
check(abs(advice.effective_radiation_mj - expected_effective) < 1e-9,
      "実質日射が「それ × フィルム劣化」になっている")

# 10MJあたり潅水量の定義どおりか
implied = advice.recommended_irrigation_l_per_m2 / advice.effective_radiation_mj * 10.0
check(abs(advice.water_per_10mj_l_per_m2 - implied) < 1e-9,
      "10MJあたり潅水量が「潅水量 ÷ 実質日射 × 10」になっている")

# =============================================================================
print()
print("=" * 78)
print("4. 日変化の形をずらしても日積算が同じなら結果が近いか")
print("=" * 78)
print("  （曇天日は山が低くなだらかになる。その影響の大きさを見る）")
print("  日の出〜日没   蒸散 [L/m²/日]   基準との差")

base = advise(forecast_from_month(12.0, 4))
print(f"  6〜17時（基準）  {base.transpiration_l_per_m2:10.3f}        —")
差 = []
for window in [(5.0, 18.0), (7.0, 16.0), (6.5, 17.5)]:
    forecast = replace(forecast_from_month(12.0, 4), daylight_hours=window)
    advice = advise(forecast)
    difference = (advice.transpiration_l_per_m2 - base.transpiration_l_per_m2) \
        / base.transpiration_l_per_m2
    差.append(abs(difference))
    print(f"  {window[0]:.1f}〜{window[1]:.1f}時    {advice.transpiration_l_per_m2:10.3f}  "
          f"{difference * 100:+9.1f}%")

check(max(差) < 0.12,
      "日射の時間帯を±1.5時間ずらしても蒸散量の差は12%以内",
      f"最大 {max(差) * 100:.1f}%")

# =============================================================================
print()
print("=" * 78)
print("5. 実データの月別平均との突き合わせ")
print("=" * 78)
print("  （左が実データを5分値で流した結果、右が日積算だけから見積もった結果）")
print()
print("  月  日射    実測ベース  予測ベース   比")

# 実データを流したときの月平均蒸散量（scratchpad/run_water_balance.py の出力）
MEASURED = {
    11: (6.49, 1.13), 12: (5.83, 0.93), 1: (7.46, 1.05), 2: (8.27, 1.36),
    3: (10.48, 1.96), 4: (11.96, 2.73), 5: (14.60, 3.39), 6: (11.78, 2.66),
}
ratios = []
for month in (11, 12, 1, 2, 3, 4, 5, 6):
    radiation, measured = MEASURED[month]
    predicted = advise(forecast_from_month(radiation, month)).transpiration_l_per_m2
    ratio = predicted / measured
    ratios.append(ratio)
    print(f"  {month:2d}  {radiation:5.2f}  {measured:9.2f}  {predicted:10.2f}  {ratio:6.2f}")

mean_ratio = sum(ratios) / len(ratios)
print(f"\n  比の平均 {mean_ratio:.2f}、幅 {min(ratios):.2f}〜{max(ratios):.2f}")
check(0.7 < mean_ratio < 1.3,
      "日積算だけの見積りが、5分値を流した結果の±30%に収まる",
      f"平均 {mean_ratio:.2f}倍")

# =============================================================================
print()
print("=" * 78)
print("5-2. 日付ベースの推定が日ごとに滑らかに変わるか")
print("=" * 78)
print("  （月区切りだと上旬と下旬が同じ値になってしまう）")
print()
print("  日付    平年日射  推定気温  推定湿度  推定CO2")

import datetime as _dt

samples = []
for month, day in [(4, 5), (4, 15), (4, 25), (5, 5), (5, 15), (5, 25),
                   (11, 5), (11, 25), (1, 15)]:
    date = _dt.date(2026, month, day)
    doy = date.timetuple().tm_yday
    normal = normal_radiation_mj(doy)
    estimate = forecast_from_date(normal, date)
    samples.append((month, day, estimate))
    print(f"  {month:2d}/{day:02d}  {normal:8.2f}  {estimate.mean_temp_c:8.1f}  "
          f"{estimate.mean_relative_humidity:8.3f}  {estimate.mean_co2_ppm:7.0f}")

# 同じ月の上旬と下旬で、推定値がちゃんと違うこと
april_early = next(e for m, d, e in samples if (m, d) == (4, 5))
april_late = next(e for m, d, e in samples if (m, d) == (4, 25))
temp_gap = abs(april_late.mean_temp_c - april_early.mean_temp_c)
rh_gap = abs(april_late.mean_relative_humidity - april_early.mean_relative_humidity)
print()
print(f"  4月上旬と下旬の差: 気温 {temp_gap:.2f} ℃ / 湿度 {rh_gap:.3f}")
check(temp_gap > 0.3, "同じ月でも上旬と下旬で推定気温が違う", f"{temp_gap:.2f} ℃")
check(rh_gap > 0.01, "同じ月でも上旬と下旬で推定湿度が違う", f"{rh_gap:.3f}")

# 日をまたいで滑らかか（隣り合う日で飛ばない）
jumps = []
for doy in range(1, 366):
    a = forecast_from_day_of_year(normal_radiation_mj(doy), doy)
    b = forecast_from_day_of_year(normal_radiation_mj(doy % 365 + 1), doy % 365 + 1)
    jumps.append(abs(b.mean_temp_c - a.mean_temp_c))
print(f"  隣り合う日の推定気温の差: 最大 {max(jumps):.3f} ℃")
check(max(jumps) < 0.2, "日をまたいでも推定値が飛ばない（滑らかにつながる）",
      f"最大 {max(jumps):.3f} ℃")

# 年をまたいでつながるか（12/31 と 1/1）
end_of_year = forecast_from_day_of_year(normal_radiation_mj(365), 365)
start_of_year = forecast_from_day_of_year(normal_radiation_mj(1), 1)
wrap = abs(end_of_year.mean_temp_c - start_of_year.mean_temp_c)
print(f"  12/31 と 1/1 の差: {wrap:.3f} ℃")
check(wrap < 0.3, "年をまたいでもつながる", f"{wrap:.3f} ℃")

# forecast_from_month は「その月の15日」と同じはず
for month in (1, 4, 7, 11):
    by_month = forecast_from_month(10.0, month)
    by_date = forecast_from_date(10.0, _dt.date(2026, month, 15))
    same = abs(by_month.mean_temp_c - by_date.mean_temp_c) < 1e-9
    check(same, f"forecast_from_month({month}) が15日と一致する")

# =============================================================================
print()
print("=" * 78)
print("5-3. 日変化の組み立て（2026-09-25 に直したところ）")
print("=" * 78)
print("  湿り具合は水蒸気圧で持ち、夜は暖房の効いた水準を使う。")
print("  CO2 は日射に対して指数で減衰する（直線ではない）。")
print()

# --- 夜の値が入っていること ---
sample = forecast_from_date(10.0, _dt.date(2026, 1, 15))
check(sample.mean_vapor_pressure_kpa is not None, "日中の水蒸気圧が入る",
      f"{sample.mean_vapor_pressure_kpa:.3f} kPa")
check(sample.night_temp_c is not None, "夜の気温が入る",
      f"{sample.night_temp_c:.1f} ℃")

# 夜は日中より涼しく、湿っている（＝飽差が小さい）
from core.psychrometry import saturation_vapor_pressure_kpa as _es

day_vpd = _es(sample.mean_temp_c) - sample.mean_vapor_pressure_kpa
night_vpd = _es(sample.night_temp_c) - sample.night_vapor_pressure_kpa
print(f"  1月15日: 日中の飽差 {day_vpd:.3f} kPa ／ 夜の飽差 {night_vpd:.3f} kPa")
check(sample.night_temp_c < sample.mean_temp_c, "夜は日中より涼しい",
      f"{sample.night_temp_c:.1f} < {sample.mean_temp_c:.1f} ℃")
check(night_vpd < day_vpd * 0.6, "夜の飽差は日中よりずっと小さい",
      f"{night_vpd:.3f} vs {day_vpd:.3f} kPa")
check(10.0 < sample.night_temp_c < 20.0,
      "夜の気温が暖房の設定温度あたりに収まる", f"{sample.night_temp_c:.1f} ℃")

# --- 気温と湿度を手で渡したら、水蒸気圧もその2つから作られること ---
manual = forecast_from_date(10.0, _dt.date(2026, 1, 15),
                            temp_c=22.0, relative_humidity=0.70)
expected_ea = _es(22.0) * 0.70
check(abs(manual.mean_vapor_pressure_kpa - expected_ea) < 1e-9,
      "気温と湿度を手で渡すと、水蒸気圧もその2つから作られる",
      f"{manual.mean_vapor_pressure_kpa:.4f} kPa")

# --- CO2 が日射に対して飽和すること ---
print()
print("  日射に対する CO2 の動き（1月15日）")
print("    日射[MJ]   CO2[ppm]   前の段からの下がり方")
previous_co2 = None
drops = []
for mj in [0.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0]:
    co2 = forecast_from_date(mj, _dt.date(2026, 1, 15)).mean_co2_ppm
    drop = "" if previous_co2 is None else f"{co2 - previous_co2:+8.0f}"
    if previous_co2 is not None:
        drops.append(previous_co2 - co2)
    print(f"    {mj:7.1f} {co2:10.0f}   {drop}")
    previous_co2 = co2

check(all(a > b for a, b in zip(drops, drops[1:])),
      "★明るくなるほど CO2 の下がり方がゆるむ（飽和する）",
      f"最初 {drops[0]:.0f} ppm → 最後 {drops[-1]:.0f} ppm")

dark_co2 = forecast_from_date(1.0, _dt.date(2026, 1, 15)).mean_co2_ppm
bright_co2 = forecast_from_date(20.0, _dt.date(2026, 1, 15)).mean_co2_ppm
check(dark_co2 - bright_co2 > 200.0,
      "暗い日と晴天日で CO2 が大きく違う", f"{dark_co2:.0f} → {bright_co2:.0f} ppm")
check(bright_co2 > 380.0, "晴天日でも外気濃度を下回らない", f"{bright_co2:.0f} ppm")

# --- 夜を別扱いにしたことが蒸散に効いているか ---
# 夜の気温を日中と同じにすると、夜間蒸散が増えて全体も増えるはず
from dataclasses import replace as _replace

base = forecast_from_date(12.0, _dt.date(2026, 5, 15))
with_warm_night = _replace(base, night_temp_c=base.mean_temp_c,
                           night_vapor_pressure_kpa=base.mean_vapor_pressure_kpa)
normal_night = advise(base, leaching_fraction=0.0).transpiration_l_per_m2
warm_night = advise(with_warm_night, leaching_fraction=0.0).transpiration_l_per_m2
print()
print(f"  5月15日 12MJ: 夜を実測水準にすると {normal_night:.3f} L/m²、"
      f"日中と同じにすると {warm_night:.3f} L/m²")
check(warm_night > normal_night * 1.05,
      "★夜を日中と同じ空気にすると蒸散が目に見えて増える（夜の扱いが効いている）",
      f"{(warm_night / normal_night - 1) * 100:+.0f}%")

# =============================================================================
print()
print("=" * 78)
print("6. 入力の検査")
print("=" * 78)

for call, label in [
    (lambda: advise(forecast_from_month(-1.0, 4)), "日射予測が負なら弾く"),
    (lambda: advise(forecast_from_month(10.0, 4), house="西"), "知らないハウス名を弾く"),
    (lambda: advise(forecast_from_month(10.0, 4), leaching_fraction=-0.1),
     "塩を流す上乗せが負なら弾く"),
    (lambda: forecast_from_month(10.0, 13), "月が13なら弾く"),
    (lambda: forecast_from_day_of_year(10.0, 0), "通日が0なら弾く"),
    (lambda: forecast_from_day_of_year(10.0, 400), "通日が400なら弾く"),
    (lambda: advise(replace(forecast_from_month(10.0, 4), daylight_hours=(18.0, 6.0))),
     "日の出と日没が逆なら弾く"),
]:
    try:
        call()
        check(False, label)
    except ValueError:
        check(True, label)

# =============================================================================
print()
print("=" * 78)
print("7. 気をつけることが言葉で出るか")
print("=" * 78)

heavy = advise(forecast_from_month(12.0, 4), leaching_fraction=1.5)
print(f"  塩を流す上乗せ150%（＝蒸散の2.5倍潅水）のとき:")
for warning in heavy.warnings:
    print(f"    ・{warning}")
check(len(heavy.warnings) > 0, "やりすぎたときに警告が出る",
      f"{len(heavy.warnings)}件")
check(heavy.is_over_irrigating or heavy.drainage_fraction > 0.40,
      "やりすぎが数字でも判定される",
      f"空気率 {heavy.min_air_filled_porosity * 100:.1f}% / "
      f"流亡率 {heavy.drainage_fraction * 100:.0f}%")

gentle = advise(forecast_from_month(12.0, 4), leaching_fraction=0.10)
print(f"\n  塩を流す上乗せ10% のとき: 警告 {len(gentle.warnings)}件")
check(len(gentle.warnings) <= len(heavy.warnings),
      "控えめな潅水のほうが警告が少ない")

# ハウス別の換算
central = advise(forecast_from_month(12.0, 4), house="中央")
east = advise(forecast_from_month(12.0, 4), house="東")
print(f"\n  中央 {central.recommended_irrigation_l_per_house:8.0f} L/棟")
print(f"  東   {east.recommended_irrigation_l_per_house:8.0f} L/棟")
check(abs(central.recommended_irrigation_l_per_m2
          - east.recommended_irrigation_l_per_m2) < 1e-9,
      "m²あたりは同じ（床面積が同じなので棟あたりも同じ）")

# =============================================================================
print()
print("=" * 78)
print("8. レシピとの受け渡し係数（換算係数・フィルム劣化）")
print("=" * 78)
print("  ★どちらも物理量ではなく、液肥混入機レシピと数字を合わせるための取り決め。")
print("  実測の透過率（core/solar.py・季節で 0.62〜0.72）とは別物。")

base_forecast = forecast_from_month(12.0, 4)
base_advice = advise(base_forecast)

# 換算係数: センサー値を「実質日射」に直す割り算。小さいほど分母が大きくなる。
print()
print("   換算係数   実質日射の元[MJ]   実質日射[MJ]   10MJあたり潅水量[L/m²]")
previous_outside = 0.0
coef_ok = True
for coef in (0.55, 0.65, 0.78, 0.90):
    trial = advise(base_forecast, recipe_radiation_coef=coef)
    print(f"   {coef:8.2f}   {trial.outside_radiation_mj:16.1f}"
          f"   {trial.effective_radiation_mj:12.1f}"
          f"   {trial.water_per_10mj_l_per_m2:22.3f}")
    if previous_outside and trial.outside_radiation_mj >= previous_outside:
        coef_ok = False
    previous_outside = trial.outside_radiation_mj

check(coef_ok, "換算係数を上げると分母は小さくなる（割り算なので）")

half = advise(base_forecast, recipe_radiation_coef=0.65 / 2.0)
check(abs(half.outside_radiation_mj - base_advice.outside_radiation_mj * 2.0) < 1e-9,
      "換算係数を半分にすると分母はちょうど2倍になる",
      f"{base_advice.outside_radiation_mj:.2f} → {half.outside_radiation_mj:.2f} MJ/m²")

# 蒸散はセンサー値だけで決まるので、換算係数を変えても動かない
transpirations = [
    advise(base_forecast, recipe_radiation_coef=c).transpiration_l_per_m2
    for c in (0.55, 0.65, 0.78, 0.90)
]
check(max(transpirations) - min(transpirations) < 1e-9,
      "換算係数を変えても蒸散量は変わらない（作物が浴びる光はセンサー値そのもの）",
      f"{transpirations[0]:.4f} L/m²")

# ★最も大事な性質: 換算係数がずれても、実際に出る水は変わらない。
#   レシピは 潅水量 = need × 面積 × 実質日射 ÷ 10 で水を出す。
#   こちらが need を同じ係数で作っていれば、分子と分母で打ち消し合う。
print()
print("   換算係数を変えても、レシピが出す水は変わらないか")
print("   換算係数   10MJあたり潅水量   レシピが出す水[L/m²]")
delivered = []
for coef in (0.55, 0.65, 0.78, 0.90):
    trial = advise(base_forecast, recipe_radiation_coef=coef)
    # レシピ側の計算をそのまま真似る
    recipe_effective = base_forecast.sensor_radiation_mj / coef * 1.00
    water = trial.water_per_10mj_l_per_m2 * recipe_effective / 10.0
    delivered.append(water)
    print(f"   {coef:8.2f}   {trial.water_per_10mj_l_per_m2:16.3f}"
          f"   {water:18.4f}")
check(max(delivered) - min(delivered) < 1e-9,
      "★換算係数が何であれ、レシピが出す水は同じ（分子と分母で打ち消し合う）",
      f"{delivered[0]:.4f} L/m²")

# フィルム劣化: 実質日射（＝分母）を決める掛け算。
print()
print("   フィルム係数   実質日射[MJ]   10MJあたり潅水量[L/m²]")
for film in (0.85, 0.95, 1.00):
    trial = advise(base_forecast, film_degradation_factor=film)
    print(f"   {film:12.2f}   {trial.effective_radiation_mj:12.1f}"
          f"   {trial.water_per_10mj_l_per_m2:22.3f}")

no_film = advise(base_forecast, film_degradation_factor=1.00)
old_film = advise(base_forecast, film_degradation_factor=0.95)
check(abs(no_film.effective_radiation_mj - no_film.outside_radiation_mj) < 1e-9,
      "フィルム係数 1.00 なら実質日射 = ハウス外の日射")
check(abs(old_film.water_per_10mj_l_per_m2
          / no_film.water_per_10mj_l_per_m2 - 1.0 / 0.95) < 1e-9,
      "フィルム係数 0.95 → 1.00 で「10MJあたり潅水量」は 1/0.95 倍だけ小さくなる",
      f"{old_film.water_per_10mj_l_per_m2:.3f} → "
      f"{no_film.water_per_10mj_l_per_m2:.3f} L/m²"
      f"（{(no_film.water_per_10mj_l_per_m2 / old_film.water_per_10mj_l_per_m2 - 1) * 100:+.1f}%）")

# 実際に出る水の量は、係数を変えても変わらない（分母と分子が同じだけ動く）
check(abs(no_film.recommended_irrigation_l_per_m2
          - old_film.recommended_irrigation_l_per_m2) < 1e-9,
      "係数を変えても、潅水量そのもの（L/m²）は変わらない",
      f"{no_film.recommended_irrigation_l_per_m2:.3f} L/m²")

# 入力の検査
for bad, label in [
    ({"recipe_radiation_coef": 0.0}, "換算係数 0 を弾く"),
    ({"recipe_radiation_coef": 1.5}, "換算係数 1.5 を弾く"),
    ({"film_degradation_factor": 0.0}, "フィルム係数 0 を弾く"),
    ({"film_degradation_factor": -0.5}, "フィルム係数が負なら弾く"),
]:
    try:
        advise(base_forecast, **bad)
        check(False, label)
    except ValueError:
        check(True, label)

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
