"""過去の作業記録から施肥の実績を読み、アプリが読むJSONにする。

潅水と同じ作法。毎回 Excel を8冊開くと遅いので、先に書き出しておく。

【シートの作り】（2019〜2026年度の8冊で同一であることを確認済み）
「追肥管理 (改良)」シートは、日付が**横**に並ぶ形をしている。

    3行目   E列から日付が1日1列ずつ横に伸びる（D列は作期の合計）
    8〜17行  中央ハウスの投入資材。資材名と量が2行ずつ5組
    33〜42行 東ハウスの投入資材。同じ作り
    61〜63行 中央の N・P・K [kg/10a]（シート側の計算結果）
    67〜69行 東の N・P・K [kg/10a]
    78行〜   資材マスタ。C列=資材名 D列=N比 E列=P比 F列=K比

**2018作業記録.xlsx は対象外。**「追肥管理」という別の作りのシートしか
持っておらず、行の位置が違う。8作期あれば十分なので取り込まない。

【どちらの数字を正とするか】
このシートには同じ量が2通りで入っている。

    (a) 資材名と投入量（8〜17行・33〜42行）… 人が手で書いた一次記録
    (b) N・P・K [kg/10a]（61〜63行・67〜69行）… シート側の計算式の結果

**このツールは (a) から計算し直した値を正とする。**理由は、(b) に
計算式の抜けがあるため。2019-06-01 と 2020-04-07 の中央ハウスで、
資材（くみあい液肥1号など）が記録されているのに N(/10a) 行が 0 になっている。
計算式がその列に入っていない。

抜けの大きさは16作期ぶんで2日・各 0.3 kg/10a で、作期合計の 0.4% にすぎない。
つまり **(b) はほぼ正しいが完全ではない**。一次記録から1つの規則で
計算し直せば、この手の穴は生じない。

【それでも (b) と突き合わせる】
突き合わせは残す。行をひとつずらして読んでいても数字は「もっともらしく」
出てしまい、人間の目では気づけないからだ。潅水で「潅水時間(分)を L/m² と
して読んでいた」のと同じ型の誤りを防ぐ仕掛けである。
差が既知の抜け（1日 1 kg/10a 未満・作期合計の2%未満）を超えたら止める。

【★堆肥は施肥Nに数えない】（2026-10-01 にこう決めた）
マスタの堆肥は N比 1.0・P比 21・K比 25 で、比率としてありえない。
投入量 2.00〜4.00 も kg か t か袋数か分からない（ユーザーも
「堆肥の記録はかなり怪しい。参考程度に」と確認済み）。

そこで堆肥は **施肥Nではなく地力窒素の一部**として扱い、別枠で書き出す。
肥効の出かたからしてもそのほうが正しい。堆肥のNは有機態で、
作付けの年に無機化するのは20〜30%にすぎず、残りは何年もかけて出てくる。
「その作期に施した肥料」ではなく「土が持っている供給力」の側だ。
詳しくは config 第10-3-2節。

**この扱いが結論を左右する。**東ハウスの作期2021〜2023は基肥が堆肥だけで、
固形のN肥料を入れていない。地力窒素の下限を決める位置にいる作期なので、
扱いを間違えると答えが動く。

【★ハウス面積は 17.4 a で固定する】
2019・2020年度のファイルには 17.0 と書かれているが**記載間違い**
（2026-10-01 ユーザー確認）。ファイルの値は読んで照合するだけにし、
換算には config.HOUSE_AREA_A を使う。
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (                                      # noqa: E402
    HOUSE_AREA_A,
    NITROGEN_AS_SOIL_SUPPLY_MATERIALS,
)

# =============================================================================
# 1. 設定
# =============================================================================

XLSX_DIR = Path(r"c:\Users\kimij\Dropbox\作業記録\過去\作業記録")
OUT = Path(__file__).resolve().parent.parent / "data" / "施肥実績.json"

SHEET = "追肥管理 (改良)"

#: 取り込む年度のファイル名。2018 は別の作りなので入れない（上の説明）。
YEARS = (2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026)

#: 行の位置（1から数える）。上の「シートの作り」と対応する。
ROW_DATE = 3
ROW_MATERIAL = {"中央": 8, "東": 33}      # 資材名の1組目。以降2行ずつ5組
MATERIAL_SLOTS = 5
ROW_NPK_PER_10A = {"中央": 61, "東": 67}  # N・P・K の3行がここから
ROW_MASTER_START = 78

#: 日付が並び始める列（1から数える）。D列は合計なので使わない。
COL_FIRST_DATE = 5

#: 資材マスタの列。
COL_MASTER_NAME = 3
COL_MASTER_N = 4
COL_MASTER_P = 5
COL_MASTER_K = 6

#: 突き合わせの許容差。これを超えたら行の読み違いとみなして止める。
#:
#: 既知の抜け（計算式が入っていない日）は1日あたり最大 0.29 kg/10a なので、
#: 1日の許容を 1.0 に置けば既知の穴は通し、行ずれは捕まえられる。
#: 行を1つずらすと資材名と量の対応が崩れて桁が変わるので、必ず引っかかる。
CHECK_TOLERANCE_DAY_KG_PER_10A = 1.0
CHECK_TOLERANCE_SEASON_RATIO = 0.02

#: 1日の施肥Nとしてありうる上限 [kg/10a]。
#: 基肥は1日に10 kg/10a ほど入るので、それを通して余裕を見た値。
#: これを超えたら行の読み違いを疑う。
MAX_PLAUSIBLE_N_KG_PER_10A = 30.0

#: 単位が他の資材とそろっていない疑いのある資材。警告を出すために持つ。
SUSPECT_MATERIALS = NITROGEN_AS_SOIL_SUPPLY_MATERIALS


# =============================================================================
# 2. 読む（副作用のない関数として組む）
# =============================================================================

def cell(rows: list[tuple], row_1based: int, col_1based: int):
    """1から数える行・列で値を取る。範囲外なら None。"""
    row = rows[row_1based - 1] if row_1based - 1 < len(rows) else ()
    index = col_1based - 1
    return row[index] if index < len(row) else None


def read_dates(rows: list[tuple]) -> dict[int, dt.date]:
    """3行目から「列番号 → 日付」を作る。"""
    header = rows[ROW_DATE - 1]
    dates = {}
    for index in range(COL_FIRST_DATE - 1, len(header)):
        value = header[index]
        if isinstance(value, dt.datetime):
            dates[index + 1] = value.date()
        elif isinstance(value, dt.date):
            dates[index + 1] = value
    return dates


def read_master(rows: list[tuple]) -> dict[str, dict[str, float]]:
    """資材マスタ（名前 → N比・P比・K比）を読む。

    空欄は 0 とみなす。その成分を含まない資材という意味だ。
    """
    master = {}
    for row_number in range(ROW_MASTER_START, len(rows) + 1):
        name = cell(rows, row_number, COL_MASTER_NAME)
        if not isinstance(name, str) or not name.strip():
            continue

        def ratio(col: int) -> float:
            value = cell(rows, row_number, col)
            return float(value) if isinstance(value, (int, float)) else 0.0

        master[name.strip()] = {
            "N": ratio(COL_MASTER_N),
            "P": ratio(COL_MASTER_P),
            "K": ratio(COL_MASTER_K),
        }
    return master


def read_area_a(rows: list[tuple]) -> dict[str, float]:
    """1行目のハウス面積 [a] を読む。中央=E列・東=H列。"""
    out = {}
    for house, col in (("中央", 5), ("東", 8)):
        value = cell(rows, 1, col)
        if not isinstance(value, (int, float)) or value <= 0:
            raise ValueError(
                f"{house}ハウスの面積が1行目から読めなかった（値: {value!r}）。"
                f"シートの作りが変わった可能性がある。"
            )
        out[house] = float(value)
    return out


def read_materials(
    rows: list[tuple], house: str, dates: dict[int, dt.date],
) -> list[tuple[dt.date, str, float]]:
    """その日に入れた資材を (日付, 資材名, 量) の並びで返す。"""
    base = ROW_MATERIAL[house]
    out = []
    for col, date in dates.items():
        for slot in range(MATERIAL_SLOTS):
            name = cell(rows, base + slot * 2, col)
            quantity = cell(rows, base + slot * 2 + 1, col)
            if (isinstance(name, str) and name.strip()
                    and isinstance(quantity, (int, float)) and quantity):
                out.append((date, name.strip(), float(quantity)))
    return out


def read_npk(
    rows: list[tuple], house: str, dates: dict[int, dt.date],
) -> dict[dt.date, tuple[float, float, float]]:
    """シート側が計算した N・P・K [kg/10a] を日付ごとに返す。"""
    base = ROW_NPK_PER_10A[house]
    out = {}
    for col, date in dates.items():
        values = []
        for offset in range(3):          # N・P・K の3行
            value = cell(rows, base + offset, col)
            values.append(float(value) if isinstance(value, (int, float)) else 0.0)
        if any(values):
            out[date] = tuple(values)
    return out


# =============================================================================
# 3. 検算
# =============================================================================

def recompute_npk(
    materials: list[tuple[dt.date, str, float]],
    master: dict[str, dict[str, float]],
    area_a: float = HOUSE_AREA_A,
    include_soil_supply_materials: bool = False,
) -> dict[dt.date, tuple[float, float, float]]:
    """資材の量とマスタの成分比から、N・P・K [kg/10a] を計算する。

    これが施肥実績の正の値。1つの規則で通すので、シート側の
    計算式の抜けに影響されない。

    Args:
        materials: (日付, 資材名, 量) の並び
        master: 資材マスタ
        area_a: ハウス面積 [a]。既定は config.HOUSE_AREA_A（17.4）。
            ★ファイルに書かれた面積は使わない（2019・2020年度は記載間違い）。
        include_soil_supply_materials:
            False（既定）… 堆肥などを除く。これが「施肥N」
            True          … 含める。堆肥由来を別に出すときに使う
    """
    out: dict[dt.date, list[float]] = {}
    for date, name, quantity in materials:
        if name not in master:
            raise KeyError(
                f"資材マスタに無い資材が使われている: '{name}'（{date}）。"
                f"マスタ（78行目以降）に行が足りないか、名前の書き方が違う。"
            )
        is_soil_supply = name in NITROGEN_AS_SOIL_SUPPLY_MATERIALS
        if is_soil_supply != include_soil_supply_materials:
            continue
        bucket = out.setdefault(date, [0.0, 0.0, 0.0])
        for index, element in enumerate(("N", "P", "K")):
            bucket[index] += quantity * master[name][element] / area_a * 10.0
    return {date: tuple(values) for date, values in out.items()}


def compare_with_sheet(
    sheet_npk: dict[dt.date, tuple[float, float, float]],
    own_npk: dict[dt.date, tuple[float, float, float]],
    label: str,
) -> list[tuple[dt.date, float, float]]:
    """シート側のN量と突き合わせ、食い違う日を返す。大きければ止める。

    Returns:
        [(日付, シートの値, 自分の値)] の並び。食い違いが無ければ空。
    """
    gaps = []
    for date in sorted(set(sheet_npk) | set(own_npk)):
        sheet_value = sheet_npk.get(date, (0.0, 0.0, 0.0))[0]
        own_value = own_npk.get(date, (0.0, 0.0, 0.0))[0]
        if abs(sheet_value - own_value) > 0.02:
            gaps.append((date, sheet_value, own_value))

    for date, sheet_value, own_value in gaps:
        if abs(sheet_value - own_value) > CHECK_TOLERANCE_DAY_KG_PER_10A:
            raise ValueError(
                f"{label} {date}: シートのN量 {sheet_value:.3f} と、"
                f"資材から計算した {own_value:.3f} kg/10a が "
                f"{CHECK_TOLERANCE_DAY_KG_PER_10A} 以上離れている。"
                f"行の位置を読み違えているか、マスタの成分比が変わった。"
            )

    sheet_total = sum(values[0] for values in sheet_npk.values())
    own_total = sum(values[0] for values in own_npk.values())
    if sheet_total > 0:
        ratio = abs(own_total - sheet_total) / sheet_total
        if ratio > CHECK_TOLERANCE_SEASON_RATIO:
            raise ValueError(
                f"{label}: 作期合計のN量が合わない。"
                f"シート {sheet_total:.2f} / 自分の計算 {own_total:.2f} kg/10a "
                f"（差 {ratio * 100:.1f}%、許容 "
                f"{CHECK_TOLERANCE_SEASON_RATIO * 100:.0f}%）。"
            )

    if gaps:
        print(f"      シート側に計算式の抜け {len(gaps)} 日 "
              f"（計 {own_total - sheet_total:+.2f} kg/10a）: "
              f"{', '.join(str(d) for d, _, _ in gaps)}")
    else:
        print(f"      シート側と一致（計 {own_total:.2f} kg/10a）")
    return gaps


def season_of(date: dt.date) -> int:
    """作期の年を返す。作期は8月始まり（8月以降はその年、7月以前は前年）。"""
    return date.year if date.month >= 8 else date.year - 1


# =============================================================================
# 4. 書き出す（ここだけが副作用）
# =============================================================================

def main() -> None:
    master: dict[str, dict[str, float]] = {}
    area_by_year: dict[int, dict[str, float]] = {}
    daily: dict[str, list[list]] = {"中央": [], "東": []}
    breakdown: dict[str, list[list]] = {"中央": [], "東": []}
    compost: dict[str, list[list]] = {"中央": [], "東": []}
    suspect_seasons: set[int] = set()
    sheet_gaps: list[str] = []

    for year in YEARS:
        path = XLSX_DIR / f"{year}作業記録.xlsx"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} が無い。作業記録は年度ごとに別ファイルなので、"
                f"置き場所が変わっていないか確かめること。"
            )

        workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
        if SHEET not in workbook.sheetnames:
            workbook.close()
            raise ValueError(
                f"{path.name} に『{SHEET}』シートが無い。"
                f"あるのは {workbook.sheetnames}。"
            )
        rows = list(workbook[SHEET].iter_rows(values_only=True))
        workbook.close()

        dates = read_dates(rows)
        if not dates:
            raise ValueError(
                f"{path.name}: 3行目から日付が1つも読めなかった。"
                f"シートの作りが変わった可能性がある。"
            )

        # ★ファイルに書かれた面積は照合だけに使う。
        #   2019・2020年度は 17.0 と書かれているが記載間違い
        #   （2026-10-01 ユーザー確認。正しくは両ハウスとも 17.4 a）。
        #   換算には config.HOUSE_AREA_A を使う。
        year_master = read_master(rows)
        year_area = read_area_a(rows)
        area_by_year[year] = year_area
        for house, value in year_area.items():
            if abs(value - HOUSE_AREA_A) > 1e-9:
                print(f"  ※ ファイルの{house}の面積 {value} a は "
                      f"{HOUSE_AREA_A} a と違う（記載間違い）。"
                      f"換算には {HOUSE_AREA_A} を使う。")

        # 資材マスタは年度で成分比が書き換わることがある。新しい年を採る。
        for name, ratios in year_master.items():
            if name in master and master[name] != ratios:
                print(f"  ※ 資材『{name}』の成分比が前年と違う "
                      f"{master[name]} → {ratios}。新しい年の値を採る。")
            master[name] = ratios

        print(f"{path.name}  日付 {len(dates)} 列  "
              f"{min(dates.values())} 〜 {max(dates.values())}")

        for house in ("中央", "東"):
            materials = read_materials(rows, house, dates)
            sheet_npk = read_npk(rows, house, dates)

            # (a) 行の読み方の検査。シートと同じ条件（ファイルの面積・堆肥を含む）で
            #     計算し直して突き合わせる。ここが合えば行の位置は正しい。
            as_sheet = {}
            for include in (False, True):
                part = recompute_npk(materials, year_master,
                                     year_area[house], include)
                for date, values in part.items():
                    merged = as_sheet.setdefault(date, [0.0, 0.0, 0.0])
                    for index in range(3):
                        merged[index] += values[index]
            as_sheet = {date: tuple(v) for date, v in as_sheet.items()}
            gaps = compare_with_sheet(sheet_npk, as_sheet, f"{year} {house}")
            sheet_gaps.extend(
                [f"{year} {house} {date}" for date, _, _ in gaps])

            # (b) 書き出す値。面積は 17.4 固定、堆肥は除く（＝施肥N）。
            own_npk = recompute_npk(materials, year_master)
            # (c) 堆肥由来は別枠。地力窒素の側に数える。
            compost_npk = recompute_npk(
                materials, year_master, include_soil_supply_materials=True)
            for date, (n, _p, _k) in sorted(compost_npk.items()):
                compost[house].append([
                    date.timetuple().tm_yday, date.year, round(n, 4)])

            for date, (n, p, k) in sorted(own_npk.items()):
                if n > MAX_PLAUSIBLE_N_KG_PER_10A:
                    raise ValueError(
                        f"{year} {house} {date}: 1日の施肥Nが "
                        f"{n:.1f} kg/10a で、上限 "
                        f"{MAX_PLAUSIBLE_N_KG_PER_10A} を超えた。"
                        f"行の読み違いを疑うこと。"
                    )
                daily[house].append([
                    date.timetuple().tm_yday, date.year,
                    round(n, 4), round(p, 4), round(k, 4),
                ])

            for date, name, quantity in materials:
                breakdown[house].append([
                    date.timetuple().tm_yday, date.year, name, quantity,
                ])
                if name in SUSPECT_MATERIALS:
                    suspect_seasons.add(season_of(date))

    # --- 作期ごとの合計。基肥（8月）と追肥（9月以降）を分ける ---
    totals: dict[str, dict[str, dict[str, float]]] = {}
    for house in ("中央", "東"):
        per_season: dict[str, dict[str, float]] = {}
        for doy, year, n, p, k in daily[house]:
            date = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
            key = str(season_of(date))
            bucket = per_season.setdefault(
                key, {"基肥N": 0.0, "追肥N": 0.0, "N": 0.0, "P": 0.0, "K": 0.0})
            # 実績では8月の投入は定植前の1日に集中している（基肥）。
            bucket["基肥N" if date.month == 8 else "追肥N"] += n
            bucket["N"] += n
            bucket["P"] += p
            bucket["K"] += k
        # 堆肥由来Nを別枠で足す（施肥Nには入れない）
        for doy, year, n in compost[house]:
            date = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
            bucket = per_season.setdefault(
                str(season_of(date)),
                {"基肥N": 0.0, "追肥N": 0.0, "N": 0.0, "P": 0.0, "K": 0.0})
            bucket["堆肥N"] = bucket.get("堆肥N", 0.0) + n

        totals[house] = {
            season: {name: round(value, 3) for name, value in bucket.items()}
            for season, bucket in sorted(per_season.items())
        }

    # --- ★記録漏れの検査 ---
    #
    # 同じ作期で、片方のハウスには基肥Nがあり、もう片方が0なら記録漏れを疑う。
    # 両ハウスは同じ作型・同じ作業で動かしているので、片方だけ基肥を
    # 入れないことは実際には起きない（2026-10-01 ユーザー確認）。
    #
    # 実例: 東ハウスの作期2021〜2023は基肥が堆肥だけで、固形のN肥料の
    # 記録が無い。シート側のハウス別集計も「東計 0」なので読み落としではなく、
    # 記録そのものが抜けている。
    #
    # ★この作期を地力窒素の逆算に使ってはいけない。施肥Nが実際より
    #   少なく出るので、土が出したN量を過大に見積もることになる。
    incomplete: list[dict] = []
    for season in sorted({s for h in totals.values() for s in h}):
        base = {
            house: totals[house].get(season, {}).get("基肥N", 0.0)
            for house in ("中央", "東")
        }
        for house, value in base.items():
            other = "東" if house == "中央" else "中央"
            if value <= 0.0 < base[other]:
                incomplete.append({
                    "作期": int(season),
                    "ハウス": house,
                    "理由": (
                        f"基肥Nの記録が無い（{other}は "
                        f"{base[other]:.2f} kg/10a）。記録漏れと判断する"
                    ),
                })

    payload = {
        "説明": (
            "作業記録『追肥管理 (改良)』シートの施肥実績。1件が1日ぶんで "
            "[通日, 年, N kg/10a, P kg/10a, K kg/10a]。"
            "N収支の答え合わせに使う。"
        ),
        "出典": (
            "Dropbox/作業記録/過去/作業記録/{2019〜2026}作業記録.xlsx "
            "追肥管理 (改良) シート"
        ),
        "読み方の注意": (
            "N・P・K はシートの計算式の結果ではなく、資材名と投入量から"
            "計算し直した値（シート側には計算式の抜けがある）。"
            "堆肥は施肥Nに含めず『堆肥N』として別枠にしてある"
            "（記録の単位が怪しく、肥効も地力窒素の側だから。config 第10-3-2節）。"
            "2018作業記録.xlsx は『追肥管理』という別の作りのシートなので"
            "取り込んでいない。"
        ),
        "シート側に計算式の抜けがあった日": sheet_gaps,
        "作成日": dt.date.today().isoformat(),
        "形式": ["通日", "年", "N kg/10a", "P kg/10a", "K kg/10a"],
        "基肥の定義": "8月中の投入を基肥とみなす。実績では8/16〜8/23の1日に集中している。",
        "ハウス面積_a": HOUSE_AREA_A,
        "ファイルに書かれた面積_a": {
            str(year): area for year, area in area_by_year.items()},
        "面積の注意": (
            f"換算にはすべて {HOUSE_AREA_A} a を使った。"
            "2019・2020年度のファイルには 17.0 と書かれているが記載間違い"
            "（2026-10-01 ユーザー確認）。"
        ),
        "作期": sorted({int(s) for h in totals.values() for s in h}),
        "資材マスタ": master,
        "堆肥を含む作期": sorted(suspect_seasons),
        "★記録が不完全な作期": incomplete,
        "記録が不完全な作期の扱い": (
            "基肥Nの記録が片方のハウスにしか無い作期。両ハウスは同じ作型で"
            "動かしているので、片方だけ基肥を入れないことは実際には起きない。"
            "施肥Nが実際より少なく出るので、地力窒素の逆算には使わないこと。"
        ),
        "ハウス": daily,
        "堆肥由来N": compost,
        "資材の内訳": breakdown,
        "作期の合計": totals,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    print(f"\n書き出した: {OUT}  ({OUT.stat().st_size / 1024:.0f} KB)")

    # --- 確認用の表示 ---
    print("\n作期ごとの施肥N [kg/10a]（基肥＋追肥）")
    print("  ※ 作期は8月始まり。作期2025 = 2025年8月〜2026年7月 "
          "= ファイル名の 2026作業記録.xlsx")
    print("  ※ 堆肥は施肥Nに含めず、別枠（堆肥N）にしてある")
    flagged = {(row["作期"], row["ハウス"]) for row in incomplete}
    print(f"{'作期':>6}{'ファイル':>9}"
          + "".join(f"{h + '基肥':>9}{h + '追肥':>9}{h + '計':>8}{h + '堆肥':>8}"
                    for h in ("中央", "東")))
    for season in payload["作期"]:
        line = f"{season:>6}{season + 1:>9}"
        for house in ("中央", "東"):
            bucket = totals[house].get(season if isinstance(season, str)
                                       else str(season))
            if bucket is None:
                line += f"{'—':>9}{'—':>9}{'—':>8}{'—':>8}"
            else:
                mark = "★" if (season, house) in flagged else ""
                line += (f"{bucket['基肥N']:>8.2f}{mark:<1}"
                         f"{bucket['追肥N']:>9.2f}{bucket['N']:>8.2f}"
                         f"{bucket.get('堆肥N', 0.0):>8.2f}")
        print(line)

    if incomplete:
        print("\n★ 記録が不完全な作期（地力窒素の逆算に使ってはいけない）")
        for row in incomplete:
            print(f"    作期{row['作期']} {row['ハウス']}ハウス: {row['理由']}")

    if suspect_seasons:
        print(f"\n※ 堆肥を含む作期: {sorted(suspect_seasons)}。"
              f"堆肥Nは別枠にしてあるので施肥Nの数字は汚れていない。"
              f"地力窒素の側で参考値として扱うこと。")


if __name__ == "__main__":
    main()
