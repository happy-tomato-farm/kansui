"""作業日誌から潅水量の実績を読み、アプリが読むJSONにする。

毎回 Excel を開くと起動が遅くなるので、先に書き出しておく。

【★列の読み方に注意（2026-09-22 に直した）】
作業日誌は1年ぶん34列のブロックが横に並んでいる。
ブロックによって潅水の列の作りが違う。

    古い年（2018-08〜2020-07）  14列目=東の潅水時間(分)   15列目=東のL/m²
                               16列目=中央の潅水時間(分) 17列目=中央のL/m²
    新しい年（2020-08以降）     14列目=15列目=東のL/m²
                               16列目=17列目=中央のL/m²

つまり **15列目と17列目を読めば全年で L/m² になる**。
以前は14・16列目を読んでいたので、2018〜2020年だけ「分」が
L/m² として混ざっていた（5月の最大が 30 L/m² と出ていたのはこれ）。

日ごとの日射（ハウス内センサー）も一緒に書き出す。
アプリ側で「同じくらいの明るさだった日」だけを選んで比べるために使う。
"""
import datetime as dt
import json
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

XLSX = r"c:\Users\kimij\Dropbox\作業記録\過去\作業記録\2026作業記録.xlsx"
OUT = Path(__file__).resolve().parent.parent / "data" / "潅水実績.json"
BLOCK = 34

#: ブロックの先頭からの位置。★東=15・中央=17（L/m²）。上のコメントを読むこと。
COL = {"date": 1, "東": 15, "中央": 17, "sensor": 21}

#: 潅水量として受けつける上限 [L/m²/日]。これを超えたら列の読み違いを疑う。
MAX_PLAUSIBLE_L_PER_M2 = 25.0

#: ここより前の作期は使わない（作期は8月始まり。2023 なら 2023年8月〜2024年7月）。
#:
#: 【なぜ古い年を切るか】（2026-09-22 ユーザー確認）
#: 記録の取り方も潅水のやり方も当時とは違うので、並べても参考にならない。
#: ついでに、列の作りが違って読み違えやすい 2018〜2020年のブロックも外れる。
#:
#: 【使える作期】
#:   2023-08〜2024-07  248日
#:   2024-08〜2025-07  263日
#:   2025-08〜2026-07  248日
#:   ※2021-08〜2023-07 はファイルに無い（欠測）ので、実質3作期
FIRST_SEASON = 2023

workbook = openpyxl.load_workbook(XLSX, data_only=True, read_only=True)
rows = list(workbook["作業日誌"].iter_rows(values_only=True))
width = max(len(row) for row in rows)

records = []
for base in range(0, width, BLOCK):
    for row in rows:
        def cell(key):
            index = base + COL[key]
            return row[index] if index < len(row) else None

        date = cell("date")
        if isinstance(date, dt.datetime):
            records.append({
                "date": date.date(),
                "東": cell("東"),
                "中央": cell("中央"),
                "sensor": cell("sensor"),
            })

frame = pd.DataFrame(records).drop_duplicates("date").sort_values("date")
for column in ("東", "中央", "sensor"):
    frame[column] = pd.to_numeric(frame[column], errors="coerce")

# ハウス日射量の単位が年度で変わっている（j/cm² の年と MJ/m² の年がある）。
# 0.15〜2380 の幅があるので、100 を超える行は j/cm² とみなして 0.01 を掛ける。
frame["sensor_mj"] = np.where(
    frame["sensor"] > 100, frame["sensor"] * 0.01, frame["sensor"])

frame["doy"] = pd.to_datetime(frame["date"]).dt.dayofyear
frame["year"] = pd.to_datetime(frame["date"]).dt.year

# 作期は8月始まり。8月以降はその年、7月以前は前年の作期に属する。
months = pd.to_datetime(frame["date"]).dt.month
frame["season"] = np.where(months >= 8, frame["year"], frame["year"] - 1)

print(f"読み取り {len(frame)} 日  {frame['date'].min()} 〜 {frame['date'].max()}")
before = len(frame)
frame = frame[frame["season"] >= FIRST_SEASON].copy()
print(f"{FIRST_SEASON}年8月以降にしぼる: {before} → {len(frame)} 日")
if len(frame) == 0:
    raise ValueError(
        f"FIRST_SEASON = {FIRST_SEASON} にすると使える日が1日も残らない。"
        f"作業日誌の日付を確かめること。"
    )

payload_houses = {}
for house in ("中央", "東"):
    watered = frame[frame[house].notna() & (frame[house] > 0)]

    # 読み違いの検査。分を L/m² と取り違えると、ここで引っかかる。
    too_large = watered[watered[house] > MAX_PLAUSIBLE_L_PER_M2]
    if len(too_large) > 0:
        raise ValueError(
            f"{house}: 潅水量が {MAX_PLAUSIBLE_L_PER_M2} L/m² を超える日が "
            f"{len(too_large)} 日ある（最大 {watered[house].max():.1f}）。"
            f"列の読み違いの可能性が高い。"
            f"最初の数日: {list(too_large['date'].head(3))}"
        )

    with_sun = watered[watered["sensor_mj"].notna() & (watered["sensor_mj"] > 0)]
    print(f"  {house}: 潅水の記録 {len(watered)} 日"
          f"（うち日射もそろう {len(with_sun)} 日）"
          f"  {watered['date'].min()} 〜 {watered['date'].max()}")
    print(f"        中央値 {watered[house].median():.2f} L/m²"
          f"  最大 {watered[house].max():.2f}")

    # 日射がそろっている日だけを書き出す。
    # 明るさで絞って比べるのが目的なので、日射のない日は使えない。
    payload_houses[house] = [
        [int(row["doy"]), int(row["year"]),
         round(float(row["sensor_mj"]), 2), round(float(row[house]), 2)]
        for _, row in with_sun.iterrows()
    ]

payload = {
    "説明": (
        "作業日誌の潅水量の実績。1件が1日ぶんで "
        "[通日, 年, ハウス内日射 MJ/m², 潅水量 L/m²]。"
        "アプリ側で「同じ時期・同じくらいの明るさ」の日だけを選んで比べる。"
    ),
    "出典": "Dropbox/作業記録/過去/作業記録/2026作業記録.xlsx 作業日誌シート",
    "列の読み方": (
        "東=15列目・中央=17列目（L/m²）。古い年は14・16列目が潅水時間(分)なので"
        "読んではいけない。"
    ),
    "作成日": dt.date.today().isoformat(),
    "形式": ["通日", "年", "日射MJ", "潅水L"],
    "期間": f"{frame['date'].min()} 〜 {frame['date'].max()}",
    "最初の作期": FIRST_SEASON,
    "作期": sorted(int(s) for s in frame["season"].unique()),
    "ハウス": payload_houses,
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
print(f"\n書き出した: {OUT}  ({OUT.stat().st_size / 1024:.0f} KB)")

# --- 確認: 月ごと・明るさごとの代表値 ---
print("\n月ごとの潅水量（中央・日射の明るさで3つに分けて）")
print("  ※ 明るさは「その月のハウス内日射の上位3分の1／中／下位3分の1」")
print("   月   日数    明るい日   ふつう    暗い日")
central = frame[frame["中央"].notna() & (frame["中央"] > 0)
                & frame["sensor_mj"].notna() & (frame["sensor_mj"] > 0)].copy()
central["month"] = pd.to_datetime(central["date"]).dt.month
for month in range(1, 13):
    group = central[central["month"] == month]
    if len(group) < 10:
        continue
    low, high = group["sensor_mj"].quantile([1 / 3, 2 / 3])
    bright = group[group["sensor_mj"] >= high]["中央"].median()
    middle = group[(group["sensor_mj"] >= low)
                   & (group["sensor_mj"] < high)]["中央"].median()
    dark = group[group["sensor_mj"] < low]["中央"].median()
    print(f"   {month:2d}  {len(group):5d}   {bright:8.2f}  {middle:8.2f}  {dark:8.2f}")
