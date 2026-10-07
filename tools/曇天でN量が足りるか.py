"""曇天でN量が足りなくなるかを、実測の日射と潅水量で確かめる。

【問い】（2026-10-07 ユーザー）
　いまの設計は、その日の受光量に比例してN量を配っている。
　曇天日は潅水量を減らし、ECをそろえるためN量も減る。
　この「減り」が効きすぎて、作期を通した必要N量に届かなくなることはないか。
　もしそうなら、快晴日のN量を増やす調整が要るのではないか。

【なぜ起こりうるのか】
　収量は毎年 30〜32 kg/m² で安定している。天気が悪くても落ちていない。
　つまり**N需要は毎年ほぼ一定**のはず。
　ところが配分は日射に比例しているので、暗い年は配る量が減る。
　この2つがずれていれば、暗い年ほど足りなくなる。

【漏れの道は3つある】
　(1) 暗い年   … 受光量が平年より少ないと、取り分の合計が1に届かない
　(2) 潅水ゼロの日 … 水を出さない日はN量もゼロになる
　(3) ECの帯   … 上限で切られた日はその分だけ減る（10月は35%カット）

【使うデータ】
　data/潅水実績.json … 作期2023〜2025の日別 [通日, 年, 日射MJ, 潅水L/m²]
　これは**実測**なので、仮定を置かずに数えられる。

使い方:
    py tools\曇天でN量が足りるか.py
"""
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import (                                            # noqa: E402
    FEED_EC_BAND_BY_MONTH,
    FERTILIZER_N_EFFICIENCY,
    FERTILIZER_N_FRACTION_DEFAULT,
    SOIL_N_SUPPLY_KG_PER_10A,
    TARGET_FRESH_YIELD_KG_PER_M2,
)
from core.advisor import normal_radiation_mj                     # noqa: E402
from core.nitrogen import (                                      # noqa: E402
    _before_planting,
    advise_fertilizer_n,
    daily_light_capture,
    season_light_capture,
)
from core.salinity import clamp_daily_n_to_ec_band               # noqa: E402

HISTORY_PATH = Path(__file__).resolve().parents[1] / "data" / "潅水実績.json"
HOUSE = "中央"


def date_of(day_of_year: int, year: int) -> dt.date:
    """通日と年から日付を作る。"""
    return dt.date(year, 1, 1) + dt.timedelta(days=day_of_year - 1)


def season_of(date: dt.date) -> int:
    """作期（8月始まり）。2025年8月〜2026年7月 なら 2025。"""
    return date.year if date.month >= 8 else date.year - 1


def load_days(house: str) -> list[tuple[dt.date, float, float]]:
    """[(日付, 日射MJ, 潅水L/m²), ...] を日付順で返す。"""
    if not HISTORY_PATH.exists():
        raise FileNotFoundError(
            f"潅水実績が無い: {HISTORY_PATH}\n"
            f"`py tools/潅水実績を作り直す.py` で作ること。"
        )
    data = json.loads(HISTORY_PATH.read_text(encoding="utf-8"))
    houses = data.get("ハウス", {})
    if house not in houses:
        raise KeyError(
            f"潅水実績に '{house}' が無い。あるのは {list(houses)}。"
        )
    rows = []
    for record in houses[house]:
        day_of_year, year, radiation, water = (
            int(record[0]), int(record[1]), float(record[2]), float(record[3]))
        rows.append((date_of(day_of_year, year), radiation, water))
    return sorted(rows)


# =============================================================================
print("=" * 86)
print("曇天でN量が足りなくなるか（実測の日射と潅水量で確かめる）")
print("=" * 86)
print()
print(f"  ハウス: {HOUSE}　目標収量: {TARGET_FRESH_YIELD_KG_PER_M2} kg/m²")
print(f"  地力窒素 {SOIL_N_SUPPLY_KG_PER_10A} ／ 利用率 {FERTILIZER_N_EFFICIENCY}"
      f" ／ 肥料のN比 {FERTILIZER_N_FRACTION_DEFAULT}")
print()

days = load_days(HOUSE)
print(f"  実測データ: {len(days)} 日（{days[0][0]} 〜 {days[-1][0]}）")
print()

# 作期ごとにまとめる
seasons: dict[int, list[tuple[dt.date, float, float]]] = {}
for date, radiation, water in days:
    if _before_planting(date):
        continue
    seasons.setdefault(season_of(date), []).append((date, radiation, water))

# =============================================================================
print("=" * 86)
print("1. その作期は平年より明るかったか暗かったか")
print("=" * 86)
print("  ★記録のある日だけで比べる。記録の無い日は両方から外す（公平に比べるため）。")
print()
print(f"{'作期':>6}{'記録日数':>9}{'実測の受光量':>14}{'平年の受光量':>14}"
      f"{'明るさ':>9}{'作期の平年合計':>15}{'記録の割合':>11}")
print(f"{'':>6}{'日':>9}{'MJ/m²':>14}{'MJ/m²':>14}{'':>9}{'MJ/m²':>15}{'':>11}")

brightness: dict[int, float] = {}
coverage: dict[int, float] = {}
for season_year in sorted(seasons):
    rows = seasons[season_year]
    actual = sum(daily_light_capture(d, r) for d, r, _ in rows)
    normal = sum(
        daily_light_capture(d, normal_radiation_mj(d.timetuple().tm_yday))
        for d, _, _ in rows)
    whole = season_light_capture(season_year)
    brightness[season_year] = actual / normal if normal > 0 else float("nan")
    coverage[season_year] = normal / whole if whole > 0 else float("nan")
    print(f"{season_year:>6}{len(rows):>9}{actual:>14.0f}{normal:>14.0f}"
          f"{brightness[season_year] * 100:>8.1f}%{whole:>15.0f}"
          f"{coverage[season_year] * 100:>10.1f}%")
print()
print("  → 「明るさ」が100%より小さければ、その作期は平年より暗かったということ。")
print("    「記録の割合」は、作期のうち何割の日が記録に残っているか。")
print("    これが小さいと、下の合計は作期の一部しか見ていないことになる。")
print()

# =============================================================================
print("=" * 86)
print("2. N量はどれだけ配られるか（計画 → 帯で切る → 実際に入る）")
print("=" * 86)
print("  ★記録のある日だけの合計。作期全体ではないので、絶対値ではなく")
print("    『計画に対して何%落ちたか』を見ること。")
print()
print(f"{'作期':>6}{'計画':>10}{'帯で切った後':>14}{'帯の損失':>11}"
      f"{'潅水ゼロ':>10}{'ゼロでの損失':>14}{'残る割合':>10}")
print(f"{'':>6}{'kg-N/10a':>10}{'kg-N/10a':>14}{'':>11}{'日':>10}"
      f"{'kg-N/10a':>14}{'':>10}")

results = {}
for season_year in sorted(seasons):
    rows = seasons[season_year]
    planned = 0.0        # 受光量に比例して配る量（帯を通す前）
    delivered = 0.0      # 帯で切ったあと、実際に入る量
    zero_water_loss = 0.0
    zero_water_days = 0
    clamped_days = 0
    for date, radiation, water in rows:
        nday = advise_fertilizer_n(
            date, radiation, TARGET_FRESH_YIELD_KG_PER_M2).daily_n_kg_per_10a
        planned += nday
        if water <= 0.0:
            # 水を出さない日はNも入らない
            zero_water_days += 1
            zero_water_loss += nday
            continue
        banded = clamp_daily_n_to_ec_band(
            nday, water, date.month, FERTILIZER_N_FRACTION_DEFAULT)
        if banded.is_clamped:
            clamped_days += 1
        delivered += banded.daily_n_kg_per_10a
    band_loss = planned - delivered - zero_water_loss
    results[season_year] = {
        "planned": planned, "delivered": delivered,
        "band_loss": band_loss, "zero_loss": zero_water_loss,
        "zero_days": zero_water_days, "clamped_days": clamped_days,
        "days": len(rows),
    }
    print(f"{season_year:>6}{planned:>10.1f}{delivered:>14.1f}"
          f"{band_loss:>10.1f}{zero_water_days:>10}"
          f"{zero_water_loss:>14.1f}{delivered / planned * 100:>9.1f}%")
print()

# =============================================================================
print("=" * 86)
print("3. ★問いへの答え: 作期を通して足りるのか")
print("=" * 86)
print()
print("  ★★このデータには偏りがある★★")
print(f"    245日前後の記録に『潅水ゼロの日』が1日も入っていない。")
print("    つまり記録は**潅水した日しか残っていない**。")
print("    記録に無い15〜18%の日は、雨天・曇天で水を出さなかった日が")
print("    多いと考えるのが自然。**まさに今回問題にしている日が抜けている。**")
print()
print("  だから1つの数字では答えられない。両端で挟む。")
print()
print("   上端 … 記録に無い日も、記録のある日と同じようにNが入ったとみなす")
print("   下端 … 記録に無い日は、水を出さずNもゼロだったとみなす")
print()
print("  本当の値はこのあいだにある。")
print()

target = advise_fertilizer_n(
    dt.date(2026, 1, 15), 10.0, TARGET_FRESH_YIELD_KG_PER_M2
).season_fertilizer_kg_per_10a

print(f"{'作期':>6}{'作期の目標':>12}{'下端の実入り':>15}{'下端':>9}"
      f"{'上端の実入り':>15}{'上端':>9}")
print(f"{'':>6}{'kg-N/10a':>12}{'kg-N/10a':>15}{'':>9}{'kg-N/10a':>15}{'':>9}")
low_ratios = []
high_ratios = []
for season_year in sorted(seasons):
    r = results[season_year]
    share = coverage[season_year]
    # 下端: 記録のある日に入った分だけ。記録に無い日はゼロ扱い
    low = r["delivered"]
    # 上端: 記録に無い日も同じ調子だったとみなして引き伸ばす
    high = r["delivered"] / share if share > 0 else float("nan")
    low_ratios.append(low / target)
    high_ratios.append(high / target)
    print(f"{season_year:>6}{target:>12.1f}{low:>15.1f}{low / target * 100:>8.1f}%"
          f"{high:>15.1f}{high / target * 100:>8.1f}%")
print()

worst = min(low_ratios)
best = max(high_ratios)
print(f"  → 作期の目標に対して **{worst * 100:.0f}〜{best * 100:.0f}%** のあいだ。")
print()
if worst >= 0.95:
    print("  ◎ 下端でも足りている。いまのままでよい。")
elif worst >= 0.85:
    print("  △ 下端では1割ほど足りない。ただし下端は")
    print("    『記録に無い日はすべて潅水ゼロ』という最も厳しい見方。")
    print("    実際にはその中にも潅水した日があるはずなので、")
    print("    **本当の値はもっと上**。走りながら直す仕組み（第6節(b)）があれば確実。")
else:
    print("  ✗ 下端で明らかに足りない。快晴日のN量を増やす調整が要る。")
print()
print("  ★3作期とも平年より明るかった（113〜122%）のも、この偏りの表れ。")
print("    暗い日が記録から抜けているので、残った日は明るい側に寄る。")
print("    **『暗い年に足りなくなるか』は、このデータでは直接は確かめられない。**")
print()

# =============================================================================
print("=" * 86)
print("4. 足りないとしたら、どこで漏れているか")
print("=" * 86)
print()
print(f"{'作期':>6}{'暗さによる減り':>16}{'ECの帯':>11}{'潅水ゼロ':>11}{'合計の減り':>13}")
for season_year in sorted(seasons):
    r = results[season_year]
    share = coverage[season_year]
    # 暗さによる減り = 平年どおりなら配られたはずの量 − 実際に計画された量
    rows = seasons[season_year]
    normal_plan = sum(
        advise_fertilizer_n(
            d, normal_radiation_mj(d.timetuple().tm_yday),
            TARGET_FRESH_YIELD_KG_PER_M2).daily_n_kg_per_10a
        for d, _, _ in rows)
    dark_loss = normal_plan - r["planned"]
    total_loss = dark_loss + r["band_loss"] + r["zero_loss"]
    print(f"{season_year:>6}{dark_loss / share:>16.1f}"
          f"{r['band_loss'] / share:>11.1f}{r['zero_loss'] / share:>11.1f}"
          f"{total_loss / share:>13.1f}")
print()
print("  （作期ぶんに引き伸ばした値。単位はすべて kg-N/10a）")
print()

# =============================================================================
print("=" * 86)
print("5. 実績の施肥Nと比べる")
print("=" * 86)
print("  ★実際に入れていた量（作業記録より）と並べる。")
print("    モデルの『実入り』がこれに近ければ、設計は現実と同じ水準にある。")
print()
actual_records = {2023: 49.0, 2024: 48.9, 2025: 50.6}   # README「16作期の実績」より
print(f"{'作期':>6}{'実績の施肥N':>14}{'モデルの実入り':>18}{'比':>9}")
for season_year in sorted(seasons):
    if season_year not in actual_records:
        continue
    r = results[season_year]
    scaled = r["delivered"] / coverage[season_year]
    print(f"{season_year:>6}{actual_records[season_year]:>14.1f}"
          f"{scaled:>18.1f}{scaled / actual_records[season_year]:>9.2f}")
print()

# =============================================================================
print("=" * 86)
print("6. 快晴日のN量を増やすとしたら、どれだけか")
print("=" * 86)
print()
need = 1.0 / worst if worst > 0 else float("nan")
print(f"  下端（いちばん厳しい見方）を目標ちょうどに乗せるには "
      f"**{need:.2f} 倍**（{(need - 1) * 100:+.0f}%）。")
print(f"  上端で見るとすでに {best * 100:.0f}% あるので、"
      f"一律に上げると明るい年は出しすぎになる。")
print()
print("  ただし、増やし方には2つある。")
print()
print("   (a) 一律に増やす … 快晴日も曇天日も同じ比率で増やす。")
print("       明るい年は目標を超える。ECも上がるので帯に当たりやすくなる。")
print()
print("   (b) 走りながら直す … 作期の途中で『ここまでに入った量』と")
print("       『ここまでに入るはずだった量』を比べ、不足分を残りの日に配り直す。")
print("       暗い年だけ自動で増え、明るい年は増えない。")
print("       レシピ側に直近7日の実績を合計する仕組みが**すでにある**ので、")
print("       そこに作期の通算を足せば作れる。")
print()
print("=" * 86)
