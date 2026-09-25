"""
天文計算による日射量の検証。

【何を確かめるか】
1. 天文値が正しいか（赤緯・可照時間・南中時刻）
2. 日照予測.html と同じ数字が出るか
3. 物理として筋が通っているか（大気圏外を超えない・τ の効き方）
4. センサー基準に直したとき、実測の上限を包むか
5. 入力の検査が効くか
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import datetime as dt

from config import CLOUDY_TAU, COVER_TRANSMITTANCE, SOLAR_SITE
from core.solar import (
    clear_sky_sensor_mj,
    clear_sky_tau,
    cloudy_sensor_mj,
    cover_transmittance_for_date,
    daily_radiation_mj,
    extraterrestrial_radiation_w,
    instantaneous_radiation,
    pressure_from_altitude_kpa,
    sensor_basis_radiation_mj,
    solar_day,
    solar_day_from_date,
)

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


# =============================================================================
print("=" * 78)
print("1. 天文値")
print("=" * 78)

solstice_summer = solar_day_from_date(dt.date(2026, 6, 21))
solstice_winter = solar_day_from_date(dt.date(2025, 12, 22))
equinox_spring = solar_day_from_date(dt.date(2026, 3, 20))

print(f"  夏至  赤緯 {solstice_summer.declination_deg:+.2f}°  "
      f"可照 {solstice_summer.daylength_h:.2f} h  "
      f"南中 {solstice_summer.solar_noon_h:.2f} 時")
print(f"  冬至  赤緯 {solstice_winter.declination_deg:+.2f}°  "
      f"可照 {solstice_winter.daylength_h:.2f} h  "
      f"南中 {solstice_winter.solar_noon_h:.2f} 時")
print(f"  春分  赤緯 {equinox_spring.declination_deg:+.2f}°  "
      f"可照 {equinox_spring.daylength_h:.2f} h  "
      f"南中 {equinox_spring.solar_noon_h:.2f} 時")

# 地軸の傾きは 23.44°。夏至・冬至でここまで振れる。
check(abs(solstice_summer.declination_deg - 23.44) < 0.1,
      "夏至の赤緯が +23.44°（地軸の傾き）",
      f"{solstice_summer.declination_deg:+.2f}°")
check(abs(solstice_winter.declination_deg + 23.44) < 0.1,
      "冬至の赤緯が −23.44°",
      f"{solstice_winter.declination_deg:+.2f}°")
check(abs(equinox_spring.declination_deg) < 1.0,
      "春分の赤緯がほぼ 0°",
      f"{equinox_spring.declination_deg:+.2f}°")

# 宇都宮あたりの可照時間。大気による屈折を見ていないので、
# 実際の日の出日の入り（屈折で前後に約4分ずつ延びる）より少し短く出る。
check(14.3 < solstice_summer.daylength_h < 14.7,
      "夏至の可照時間が 14.5 時間前後",
      f"{solstice_summer.daylength_h:.2f} h")
check(9.3 < solstice_winter.daylength_h < 9.7,
      "冬至の可照時間が 9.5 時間前後",
      f"{solstice_winter.daylength_h:.2f} h")
check(abs(equinox_spring.daylength_h - 12.0) < 0.2,
      "春分の可照時間が約12時間",
      f"{equinox_spring.daylength_h:.2f} h")

# 東経139.76°は明石（135°）より4.76°東。1°で4分早くなるので約19分早い。
# 均時差（最大±16分）が乗るので 11.4〜12.0 時のあいだに収まる。
noons = [solar_day(j).solar_noon_h for j in range(1, 366)]
check(all(11.3 < n < 12.1 for n in noons),
      "南中時刻が年間を通して 11.3〜12.1 時に収まる",
      f"{min(noons):.2f}〜{max(noons):.2f} 時")

# 南中は日の出と日の入りのちょうど真ん中
check(abs((solstice_summer.sunrise_h + solstice_summer.sunset_h) / 2.0
          - solstice_summer.solar_noon_h) < 1e-9,
      "南中は日の出と日の入りの真ん中にある")

# =============================================================================
print()
print("=" * 78)
print("2. 日照予測.html と同じ数字が出るか")
print("=" * 78)
print("  （Dropbox/作業記録/アプリ/日照予測.html の式11.1〜11.13 の移植）")

check(abs(pressure_from_altitude_kpa(0.0) - 101.3) < 1e-6,
      "標高 0 m の大気圧が 101.3 kPa",
      f"{pressure_from_altitude_kpa(0.0):.4f} kPa")
check(abs(pressure_from_altitude_kpa(150.0) - 99.511) < 0.01,
      "標高 150 m の大気圧が約 99.51 kPa",
      f"{pressure_from_altitude_kpa(150.0):.3f} kPa")

# 快晴 τ の季節変化。
# 宇都宮気象台の実測341日から当てはめた値で、山は1月中旬・谷は7月中旬。
# 冬に高いのは大気が乾いて澄むため。
tau_january = clear_sky_tau(dt.date(2026, 1, 15).timetuple().tm_yday)
tau_july = clear_sky_tau(dt.date(2026, 7, 15).timetuple().tm_yday)
print(f"  1月15日の快晴τ {tau_january:.3f} / 7月15日 {tau_july:.3f}")
check(abs(tau_january - 0.867) < 0.01, "1月中旬の快晴τが 0.87 前後",
      f"{tau_january:.3f}")
check(abs(tau_july - 0.570) < 0.01, "7月中旬の快晴τが 0.57 前後",
      f"{tau_july:.3f}")
check(all(0.55 <= clear_sky_tau(j) <= 0.88 for j in range(1, 366)),
      "快晴τは年間を通して物理的にありうる範囲に収まる")

# 山が元日ではなく1月中旬にある（sin 項が入っている証拠）。
tau_by_day = [(j, clear_sky_tau(j)) for j in range(1, 366)]
peak_day = max(tau_by_day, key=lambda pair: pair[1])[0]
check(5 <= peak_day <= 30, "快晴τの山が1月上〜中旬にある",
      f"第 {peak_day} 日")

# 日積算（ハウス外・快晴）
print()
print("   月   快晴(外) MJ   曇天(外) MJ")
for month in range(1, 13):
    date = dt.date(2026, month, 15)
    print(f"   {month:2d}   {daily_radiation_mj(date):10.1f}   "
          f"{daily_radiation_mj(date, tau=CLOUDY_TAU):10.1f}")

summer_mj = daily_radiation_mj(dt.date(2026, 6, 21))
winter_mj = daily_radiation_mj(dt.date(2025, 12, 22))
check(26.0 < summer_mj < 29.0, "夏至の快晴日射が 27 MJ/m² 前後",
      f"{summer_mj:.1f} MJ/m²")
check(9.5 < winter_mj < 11.5, "冬至の快晴日射が 10 MJ/m² 前後",
      f"{winter_mj:.1f} MJ/m²")

# =============================================================================
print()
print("=" * 78)
print("3. 物理として筋が通っているか")
print("=" * 78)

# 全天日射は、大気圏外日射（物理の天井）を絶対に超えない
day = solar_day_from_date(dt.date(2026, 6, 21))
worst_ratio = 0.0
for minute in range(0, 1441, 5):
    hour = minute / 60.0
    ground = instantaneous_radiation(hour, day, tau=0.95).global_w
    space = extraterrestrial_radiation_w(hour, day)
    if space > 1.0:
        worst_ratio = max(worst_ratio, ground / space)
check(worst_ratio <= 1.0,
      "地表の全天日射は、τ を 0.95 まで上げても大気圏外を超えない",
      f"最大比 {worst_ratio:.3f}")

# τ を上げると日射は増える（単調）
taus = [0.30, 0.45, 0.60, 0.70, 0.80, 0.90]
values = [daily_radiation_mj(dt.date(2026, 5, 15), tau=t) for t in taus]
print("   τ ごとの5月15日の日積算（ハウス外）:")
for tau_value, mj in zip(taus, values):
    print(f"     τ={tau_value:.2f}  {mj:5.1f} MJ/m²")
check(all(a < b for a, b in zip(values, values[1:])),
      "τ を上げると日積算が単調に増える")

# 標高が高いほど空気が薄く、日射は強い
low = daily_radiation_mj(dt.date(2026, 5, 15), altitude_m=0.0)
high = daily_radiation_mj(dt.date(2026, 5, 15), altitude_m=2000.0)
check(high > low, "標高が高いほど日射が強い",
      f"0 m {low:.1f} → 2000 m {high:.1f} MJ/m²")

# 南中時が1日でいちばん強い
noon_w = instantaneous_radiation(day.solar_noon_h, day, tau=0.70).global_w
other = [instantaneous_radiation(h / 4.0, day, tau=0.70).global_w
         for h in range(0, 96)]
check(noon_w >= max(other) - 1e-6, "南中時が1日でいちばん日射が強い",
      f"南中 {noon_w:.0f} W/m²")

# 夜は日射ゼロ
midnight = instantaneous_radiation(0.0, day, tau=0.70)
check(midnight.global_w == 0.0, "真夜中の日射はゼロ")
check(midnight.cos_zenith <= 0.0, "真夜中は太陽が地平線の下にある")

# =============================================================================
print()
print("=" * 78)
print("4. センサー基準に直したとき、実測の上限を包むか")
print("=" * 78)
print(f"  センサー基準 = 快晴計算 × その日の透過率（水準 {COVER_TRANSMITTANCE} × 季節）")
print()

# 透過率の季節変化。実測（宇都宮気象台との突き合わせ）では
# 5月中旬が最大・11月中旬が最小で、幅は 1.13 倍。
print("   月   透過率   実測(快晴日・直近3年)")
OBSERVED_TRANSMITTANCE = {
    1: 0.684, 2: 0.654, 3: 0.736, 4: 0.715, 5: 0.720, 6: 0.725,
    9: 0.712, 10: 0.616, 11: 0.647, 12: 0.654,
}
transmittances = {}
for month in range(1, 13):
    doy = dt.date(2026, month, 15).timetuple().tm_yday
    value = cover_transmittance_for_date(doy)
    transmittances[month] = value
    observed = OBSERVED_TRANSMITTANCE.get(month)
    tail = f"{observed:.3f}" if observed else "  —  "
    print(f"   {month:2d}   {value:.3f}    {tail}")

check(transmittances[5] > transmittances[11],
      "透過率は5月のほうが11月より高い（太陽高度が高いぶん反射が少ない）",
      f"5月 {transmittances[5]:.3f} > 11月 {transmittances[11]:.3f}")
spread = max(transmittances.values()) / min(transmittances.values())
check(1.10 < spread < 1.17, "透過率の年間の幅が実測どおり 1.13 倍前後",
      f"{spread:.3f} 倍")

# 9月・10月は当てはめが実測から外れることが分かっている（config 第4-1節）。
# その2か月を除けば、実測との差は 0.04 以内に収まるはず。
worst_month, worst_gap = None, 0.0
for month, observed in OBSERVED_TRANSMITTANCE.items():
    if month in (9, 10):
        continue
    gap = abs(transmittances[month] - observed)
    if gap > worst_gap:
        worst_month, worst_gap = month, gap
check(worst_gap < 0.04, "9〜10月を除けば実測との差は 0.04 以内",
      f"最も外れたのは {worst_month}月 の {worst_gap:.3f}")
print()

#: 月ごとに最も明るかった日 [年, 月, 日, 実測の日積算 MJ/m²]（ハウス内センサー）。
#: 出典: 作業日誌「ハウス日射量」列、2024年8月〜2026年6月の632日。
#:
#: ★必ず同じ日どうしで比べること。
#:   以前は「月の最大」対「その月の15日の快晴計算」で比べていたため、
#:   月初の長い日が最大になる月（11月4日など）で4割も超えて見えていた。
#:   実際に同じ日で比べれば 1.27 倍で、それは下の説明どおりの想定内。
OBSERVED_BRIGHTEST = [
    (2025,  1, 30, 10.34), (2025,  2, 27, 12.16), (2025,  3, 21, 16.51),
    (2025,  4, 30, 20.51), (2026,  5, 30, 22.82), (2025,  6,  5, 22.12),
    (2025,  8, 21, 16.43), (2024,  9,  5, 16.85), (2024, 10,  2, 11.90),
    (2025, 11,  4, 11.47), (2024, 12, 24,  8.08),
]

# 【「快晴」が何を指すかが変わった — 読むときの注意】
# 較正係数 1.20 があったころ、この曲線は「9年でいちばん明るかった日」に
# 近い上限だった（374日のうち超えたのは7日だけ）。
# いまは τ を実測の快晴日341日の中央値に当てはめたので、
# **ふつうに晴れた日**を指す。だから飛び抜けて明るい日は曲線を超える。
# 毎朝の既定値としては、上限よりふつうの快晴日のほうが使いやすい。
# 例外的に澄んだ日は、画面で τ を上げて対応する。
print("   月    日付      快晴(計算)   実測の最大   実測/快晴")
excesses = []
for year, month, day, observed in OBSERVED_BRIGHTEST:
    date = dt.date(year, month, day)
    clear = clear_sky_sensor_mj(date)
    ratio = observed / clear
    excesses.append((month, ratio))
    print(f"   {month:2d}   {date}   {clear:8.1f}   {observed:10.1f}   {ratio:9.3f}")

worst_month, worst_ratio = max(excesses, key=lambda pair: pair[1])
check(worst_ratio < 1.30,
      "いちばん明るかった日でも快晴の計算値を3割は超えない",
      f"最大は {worst_month}月 の {worst_ratio:.3f} 倍")

# 逆に、曲線が高すぎてもいけない。どの月かで必ず超える日があるはず。
over_count = sum(1 for _, ratio in excesses if ratio > 1.0)
check(over_count >= 6,
      "半分以上の月で、最も明るい日が快晴の計算値に届く（曲線が高すぎない）",
      f"{over_count} / {len(excesses)} か月")

# 曇天は快晴より必ず暗い
for month in (1, 5, 9):
    date = dt.date(2026, month, 15)
    check(cloudy_sensor_mj(date) < clear_sky_sensor_mj(date),
          f"{month}月: 曇天は快晴より暗い",
          f"{cloudy_sensor_mj(date):.1f} < {clear_sky_sensor_mj(date):.1f} MJ/m²")

# 透過率の水準を上げるとセンサー基準の値も比例して増える
base = sensor_basis_radiation_mj(dt.date(2026, 5, 15), cover_transmittance=0.65)
higher = sensor_basis_radiation_mj(dt.date(2026, 5, 15), cover_transmittance=0.78)
check(abs(higher / base - 0.78 / 0.65) < 1e-9,
      "透過率の水準に比例してセンサー基準の日射が変わる",
      f"{base:.2f} → {higher:.2f} MJ/m²")

# 画面の換算が往復して戻るか（気象予報 → センサー基準 → 気象予報）。
#
# アプリでは「気象予報の全天日射量」を
#     センサー基準 = 予報 × その日の透過率
# で読みかえ、表示するときは逆をたどる。ここが食い違うと、
# 基準を切り替えただけで結果が変わってしまう。
print()
print("  画面の換算の往復（気象予報 → センサー基準 → 気象予報）:")
worst_roundtrip = 0.0
for month in (1, 5, 9, 12):
    date = dt.date(2026, month, 15)
    physics_outside = daily_radiation_mj(date)
    to_sensor = cover_transmittance_for_date(date.timetuple().tm_yday)
    sensor = physics_outside * to_sensor
    back = sensor / to_sensor
    worst_roundtrip = max(worst_roundtrip, abs(back - physics_outside))
    print(f"    {month:2d}月  予報 {physics_outside:5.1f} → センサー基準 {sensor:5.1f}"
          f" → 予報 {back:5.1f} MJ/m²（透過率 {to_sensor:.3f}）")
check(worst_roundtrip < 1e-9, "気象予報 ⇄ センサー基準 の換算が往復して元に戻る",
      f"最大差 {worst_roundtrip:.2e}")

# 快晴のセンサー基準値は、物理計算 × その日の透過率 とちょうど一致する
date = dt.date(2026, 5, 15)
expected = daily_radiation_mj(date) * cover_transmittance_for_date(
    date.timetuple().tm_yday)
check(abs(clear_sky_sensor_mj(date) - expected) < 1e-9,
      "快晴のセンサー基準値 = 物理計算 × その日の透過率")

# 較正係数は廃止した。引数として渡されたら黙って無視せずエラーにする。
try:
    sensor_basis_radiation_mj(date, calibration_factor=1.2)
    check(False, "廃止した較正係数を渡すとエラーになる")
except TypeError:
    check(True, "廃止した較正係数を渡すとエラーになる",
          "黙って無視されない")

# =============================================================================
print()
print("=" * 78)
print("5. 入力の検査")
print("=" * 78)

for bad_day, label in [(0, "通日 0 を弾く"), (367, "通日 367 を弾く")]:
    try:
        solar_day(bad_day)
        check(False, label)
    except ValueError as error:
        check("通日" in str(error), label)

for bad_tau, label in [(0.0, "τ = 0 を弾く"), (1.0, "τ = 1 を弾く"),
                       (-0.5, "τ が負なら弾く")]:
    try:
        daily_radiation_mj(dt.date(2026, 5, 15), tau=bad_tau)
        check(False, label)
    except ValueError as error:
        check("透過率" in str(error), label)

try:
    solar_day(180, latitude_deg=95.0)
    check(False, "緯度 95° を弾く")
except ValueError as error:
    check("緯度" in str(error), "緯度 95° を弾く")

try:
    sensor_basis_radiation_mj(dt.date(2026, 5, 15), cover_transmittance=0.0)
    check(False, "透過率 0 を弾く")
except ValueError as error:
    check("透過率" in str(error), "透過率 0 を弾く")

try:
    cover_transmittance_for_date(136, base_transmittance=0.0)
    check(False, "透過率の水準 0 を弾く")
except ValueError as error:
    check("透過率" in str(error), "透過率の水準 0 を弾く")

for bad_day, label in [(0, "透過率: 通日 0 を弾く"),
                       (400, "透過率: 通日 400 を弾く")]:
    try:
        cover_transmittance_for_date(bad_day)
        check(False, label)
    except ValueError as error:
        check("通日" in str(error), label)

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
