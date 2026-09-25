"""
作業日誌Excelから潅水実績と気象台データを読み込む。

【元データ】
Dropbox/作業記録/過去/作業記録/YYYY作業記録.xlsx の「作業日誌」シート。

このシートは年度ごとのブロックが横方向に並んでいる。
2026作業記録.xlsx の場合、B〜AH列が2026年度（2025年8月〜2026年7月）、
AJ〜BP列が2025年度、というように左から新しい順に並ぶ。

【潅水実績の扱いについて】
潅水量は L/m²/日 で記録されている。ただしこれは
「実際に必要だった量」ではなく「経験にもとづいて流した量」であり、
適正量である保証はない（利用者からの申し送り）。

さらに土耕栽培のため、次の関係になる。

    潅水量 = 蒸散量 + 土壌蒸発 + 深部浸透 ± 土壌水分の増減

土壌が水を蓄えるので、日単位では潅水量と蒸散量は一致しない。
したがって潅水実績は「正解データ」ではなく、
オーダーと季節傾向を確かめるための参考値として扱う。
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# 作業日誌シートの列位置（0始まり）。
# 年度ブロックの先頭（2026作業記録.xlsx の B列＝index 1）を基準とした相対位置。
COLUMN_OFFSETS = {
    "date": 0,               # B列 日付
    "outside_temp_mean": 2,  # D列 平均気温（宇都宮）
    "outside_temp_max": 3,   # E列 最高気温
    "outside_temp_min": 4,   # F列 最低気温
    "precipitation": 5,      # G列 降水量
    "outside_radiation": 6,  # H列 日射量（宇都宮）[MJ/m²]
    "sunrise": 7,            # I列 日の出
    "sunset": 8,             # J列 日の入
    "irrigation_east": 14,   # P列 灌水 東 [L/m²]
    "irrigation_central": 16,  # R列 灌水 中央 [L/m²]
    "house_radiation": 20,   # V列 ハウス日射量
    "effective_radiation": 21,  # W列 実質日射（外部推定）
    "inside_radiation_east": 26,   # AB列 ハウス内実質日射量予測値 東
    "inside_radiation_central": 27,  # AC列 同 中央
}

# データが始まる行（1始まり。ヘッダーが4行あり、5行目から）
FIRST_DATA_ROW = 5
LAST_DATA_ROW = 369

# 年度ブロックの先頭列（1始まり）。B列＝2。
CURRENT_YEAR_START_COLUMN = 2


def load_work_diary(
    excel_path: Path | str,
    sheet_name: str = "作業日誌",
) -> pd.DataFrame:
    """
    作業日誌シートから、最新年度のブロックを読み込む。

    Args:
        excel_path: 作業記録Excelのパス
        sheet_name: シート名

    Returns:
        date, outside_radiation, irrigation_central, irrigation_east,
        inside_radiation_central などの列を持つ DataFrame
    """
    excel_path = Path(excel_path)
    if not excel_path.exists():
        raise FileNotFoundError(f"作業記録Excelが見つかりません: {excel_path}")

    # openpyxl で読む。数式の結果（値）が欲しいので data_only=True。
    raw = pd.read_excel(
        excel_path,
        sheet_name=sheet_name,
        header=None,
        engine="openpyxl",
    )

    # 年度ブロックを切り出す（0始まりに直す）
    start = CURRENT_YEAR_START_COLUMN - 1

    rows = []
    for row_index in range(FIRST_DATA_ROW - 1, min(LAST_DATA_ROW, len(raw))):
        record = {}
        for name, offset in COLUMN_OFFSETS.items():
            column = start + offset
            if column < raw.shape[1]:
                record[name] = raw.iat[row_index, column]
            else:
                record[name] = None
        rows.append(record)

    result = pd.DataFrame(rows)

    # 日付が入っていない行（余白）は落とす
    result = result[pd.notna(result["date"])].copy()
    result["date"] = pd.to_datetime(result["date"], errors="coerce").dt.date
    result = result[pd.notna(result["date"])].reset_index(drop=True)

    # 数値であるべき列を数値に直す。
    # 文字列やメモが混ざっている場合は NaN になる（errors="coerce"）。
    numeric_columns = [
        "outside_temp_mean", "outside_temp_max", "outside_temp_min",
        "precipitation", "outside_radiation",
        "irrigation_east", "irrigation_central",
        "house_radiation", "effective_radiation",
        "inside_radiation_east", "inside_radiation_central",
    ]
    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")

    return result


def get_irrigation(
    diary: pd.DataFrame,
    house: str,
) -> pd.DataFrame:
    """
    指定したハウスの潅水実績を取り出す。

    潅水量が 0 または未記入の日は「記録なし」として除外する。
    （実際に潅水しなかったのか、記録し忘れたのかを区別できないため）

    Args:
        diary: load_work_diary() の戻り値
        house: "中央" または "東"

    Returns:
        date, irrigation_l_per_m2, outside_radiation_mj を持つ DataFrame
    """
    column_map = {"中央": "irrigation_central", "東": "irrigation_east"}
    if house not in column_map:
        raise ValueError(
            f"ハウス名は '中央' か '東' を指定してください: '{house}'"
        )

    column = column_map[house]
    inside_column = (
        "inside_radiation_central" if house == "中央" else "inside_radiation_east"
    )

    result = diary[["date", column, "outside_radiation", inside_column]].copy()
    result.columns = [
        "date", "irrigation_l_per_m2", "outside_radiation_mj",
        "inside_radiation_mj",
    ]

    # 0 と欠測を除く
    result = result[
        pd.notna(result["irrigation_l_per_m2"])
        & (result["irrigation_l_per_m2"] > 0)
    ].reset_index(drop=True)

    return result


def compare_with_model(
    daily_model: pd.DataFrame,
    irrigation: pd.DataFrame,
) -> pd.DataFrame:
    """
    モデルが出した蒸散量と、潅水実績を日付で突き合わせる。

    Args:
        daily_model: daily.aggregate_daily() の戻り値
        irrigation: get_irrigation() の戻り値

    Returns:
        両者を結合し、比（モデル ÷ 実績）を加えた DataFrame
    """
    merged = daily_model.merge(irrigation, on="date", how="inner")

    # モデルの蒸散量が潅水実績の何割にあたるか。
    # 土耕なので日単位では一致しないが、傾向を見るために計算する。
    merged["model_to_irrigation_ratio"] = (
        merged["transpiration_l_per_m2"] / merged["irrigation_l_per_m2"]
    )

    return merged


def summarize_comparison(comparison: pd.DataFrame) -> pd.DataFrame:
    """比較結果を月別にまとめる。"""
    work = comparison.copy()
    work["month"] = pd.to_datetime(work["date"]).dt.strftime("%Y-%m")

    complete = work[work["is_complete"]] if "is_complete" in work else work

    return complete.groupby("month").agg(
        日数=("date", "count"),
        モデル蒸散_L_m2=("transpiration_l_per_m2", "mean"),
        潅水実績_L_m2=("irrigation_l_per_m2", "mean"),
        比=("model_to_irrigation_ratio", "mean"),
        ハウス内日射_MJ=("daily_radiation_mj", "mean"),
        外部日射_MJ=("outside_radiation_mj", "mean"),
    ).reset_index()
