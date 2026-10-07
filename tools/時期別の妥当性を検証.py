"""施肥の時期配分（冬に濃く・春に薄く）が妥当かを、日ごとに回して確かめる。

【答える問い】（2026-10-01 ユーザー聞き取り）
デルフィのコーチから「冬は窒素を多めに流して EC を上げる、春先は EC を
下げる」という指導を受けて、いまの運用になっている。これは妥当か。

**狙っているのは EC** で、Nはその手段。だから2つを別々に見る。

    (A) N は足りているか   … 根圏の無機態Nが需要に届くか
    (B) EC は適正か        … 浸透ポテンシャルで吸水を削っていないか、
                             塩害のしきい値を超えていないか

【中央ハウスだけを見る】
東ハウスの過去の記録は当てにならない（2026-10-01 ユーザー確認）。
基肥Nの記録が抜けている作期もある（data/施肥実績.json の
「★記録が不完全な作期」）。

【データの流れ】

    作業日誌        → 日別の潅水量・日射（無潅水の日も含む）
    施肥実績.json   → 日別の施肥N
    data.js         → 作期の収量 → 作期のN需要
    advisor         → 日別の蒸散量・糖（★重いのでキャッシュする）
    water_balance   → 日別の流亡水量・含水率
    nitrogen        → 作期のN需要を、糖の日別配分で割る
    soil_nitrogen   → 根圏の無機態N・濃度
    salinity        → 給液EC → 根圏のEC → 浸透ポテンシャル → 吸水・収量・糖度

【使い方】
    py tools/時期別の妥当性を検証.py

初回は advisor を750日ぶん回すので2〜3分かかる。結果は
data/日別の蒸散と糖.json に残すので、2回目からはすぐ終わる。
作り直したいときはそのファイルを消す。
"""
from __future__ import annotations

import datetime as dt
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import openpyxl                                               # noqa: E402

from config import (                                          # noqa: E402
    EXTINCTION_COEFFICIENT_K,
    FIELD_CAPACITY_POTENTIAL_J_PER_KG,
    LEAF_AREA_PER_LEAF_M2,
    ROOT_ZONE_DEPTH_M,
    SALINITY_YIELD_THRESHOLD_ECE,
    SOIL_DIAGNOSIS_CHUO,
    SOIL_DIAGNOSIS_EC_TARGET_1TO5,
    SOIL_DIAGNOSIS_SAMPLE_DEPTH_M,
    SOIL_N_SUPPLY_KG_PER_10A,
)
from core.advisor import advise, forecast_from_date            # noqa: E402
from core.canopy import leaf_count_per_m2                     # noqa: E402
from core.nitrogen import demand_from_fresh_yield              # noqa: E402
from core.salinity import (                                    # noqa: E402
    DailySaltInput,
    ec_1to5_to_soil_solution,
    feed_ec_from_n_concentration,
    mineral_n_to_g_per_m2,
    osmotic_potential_j_per_kg,
    salinity_response,
    simulate_salt,
    soil_solution_to_saturated_extract,
    steady_state_concentration_factor,
    uptake_driving_force_ratio,
)
from core.soil import water_content_from_potential             # noqa: E402
from core.soil_nitrogen import (                              # noqa: E402
    DailyNitrogenInput,
    SoilNitrogenSettings,
    simulate,
    summarize,
    temperature_weights,
)
from core.water_balance import (                               # noqa: E402
    DailyWaterInput,
    simulate_water_balance,
)

ROOT = Path(__file__).resolve().parent.parent
SHIHI_JSON = ROOT / "data" / "施肥実績.json"
CACHE_JSON = ROOT / "data" / "日別の蒸散と糖.json"
XLSX = Path(r"c:\Users\kimij\Dropbox\作業記録\過去\作業記録\2026作業記録.xlsx")
YIELD_JS = Path(r"c:\Users\kimij\Dropbox\作業記録\data.js")

HOUSE = "中央"

#: 作業日誌の列の位置（ブロックの先頭からの位置、1から数える）。
#: tools/潅水実績を作り直す.py と同じ。★中央の潅水は17列目（15は東）。
BLOCK = 34
COL_DATE = 1
COL_IRRIGATION = 17
COL_RADIATION = 21

#: 使う作期。潅水の記録が信用できるのは2023年8月以降の3作期
#: （tools/潅水実績を作り直す.py の FIRST_SEASON と同じ理由）。
SEASONS = (2023, 2024, 2025)

#: 作期の月の並び（8月始まり）。7月は栽培していない。
SEASON_MONTHS = (8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6)

#: 使っている肥料のN比。トケル養液配合1号（液肥混入機レシピの既定）。
FERTILIZER_N_FRACTION = 0.10

#: 定植の ISO 週。config.LEAF_COUNT_PLAN_PER_M2 が値を持ち始める週。
#: これより前は作物がいないので、N需要をゼロとして扱う。
PLANTING_ISO_WEEK = 38


# =============================================================================
# 1. 日ごとのデータを集める
# =============================================================================

def season_of(date: dt.date) -> int:
    """作期の年（8月始まり）。"""
    return date.year if date.month >= 8 else date.year - 1


def before_planting(date: dt.date) -> bool:
    """定植前かどうか。

    作期は8月始まりだが、定植は9月中旬。config の LEAF_COUNT_PLAN_PER_M2 は
    ISO週38から値を持っており、週27〜37（7月上旬〜9月中旬）は栽培していない。
    その週に入るまでは作物がいないので、N需要はゼロとして扱う。

    ★基肥は8月に入るので、収支そのものは8月から回す。需要だけを止める。
    """
    return date.month in (7, 8) or (
        date.month == 9 and date.isocalendar().week < PLANTING_ISO_WEEK)


def read_diary() -> dict[dt.date, dict[str, float]]:
    """作業日誌から、日別の潅水量 [L/m²] と日射 [MJ/m²] を読む。

    ★無潅水の日も残す。潅水実績.json は潅水>0 の日だけなので、
      これで比を取ると分母から「水をやらなかった日」が抜ける
      （README「★無潅水の日を分母に入れること」）。
    """
    if not XLSX.exists():
        raise FileNotFoundError(
            f"{XLSX} が無い。作業記録は年度ごとに別ファイルなので、"
            f"置き場所が変わっていないか確かめること。"
        )
    workbook = openpyxl.load_workbook(XLSX, data_only=True, read_only=True)
    rows = list(workbook["作業日誌"].iter_rows(values_only=True))
    workbook.close()
    width = max(len(row) for row in rows)

    out: dict[dt.date, dict[str, float]] = {}
    for base in range(0, width, BLOCK):
        for row in rows:
            def cell(offset: int):
                index = base + offset
                return row[index] if index < len(row) else None

            date = cell(COL_DATE)
            if not isinstance(date, dt.datetime):
                continue
            date = date.date()

            def number(value) -> float:
                return float(value) if isinstance(value, (int, float)) else 0.0

            radiation = number(cell(COL_RADIATION))
            # ハウス日射量の単位が年度で変わっている（j/cm² の年がある）。
            # 100 を超える値は j/cm² とみなして 0.01 を掛ける
            # （tools/潅水実績を作り直す.py と同じ扱い）。
            if radiation > 100.0:
                radiation *= 0.01

            out[date] = {
                "irrigation_mm": number(cell(COL_IRRIGATION)),
                "radiation_mj": radiation,
            }
    return out


def load_fertilizer_daily() -> dict[dt.date, float]:
    """施肥実績.json から中央ハウスの日別施肥N [g-N/m²] を読む。

    ★1 g-N/m² = 1 kg-N/10a なので、JSON の kg/10a をそのまま使える。
    """
    if not SHIHI_JSON.exists():
        raise FileNotFoundError(
            f"{SHIHI_JSON} が無い。先に "
            f"`py tools/施肥実績を作り直す.py` を走らせること。"
        )
    document = json.loads(SHIHI_JSON.read_text(encoding="utf-8"))
    out: dict[dt.date, float] = {}
    for doy, year, n, _p, _k in document["ハウス"][HOUSE]:
        date = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
        out[date] = out.get(date, 0.0) + float(n)
    return out


def load_yield() -> dict[int, float]:
    """data.js から中央ハウスの作期別収量 [kg/m²] を読む。"""
    if not YIELD_JS.exists():
        raise FileNotFoundError(
            f"{YIELD_JS} が無い。Dropbox の同期を確かめること。"
        )
    text = YIELD_JS.read_text(encoding="utf-8")
    document = json.loads(text[text.find("{"):].rstrip().rstrip(";"))
    out = {}
    for season in document["seasons"]:
        record = season["houses"].get(HOUSE)
        if not record:
            continue
        kilograms = sum(float(v) for v in record["kg"])
        if kilograms > 0:
            out[int(season["start"])] = kilograms / (
                float(record["area"]) * 1000.0)
    return out


def load_or_build_cache(dates: list[dt.date],
                        radiation: dict[dt.date, float]) -> dict[str, dict]:
    """advisor を日ごとに回した結果（蒸散量・糖）をキャッシュ越しに返す。

    advisor は1日に 0.2 秒かかる。750日で2〜3分なので、一度回したら
    JSON に残す。作り直したいときはファイルを消す。
    """
    cache: dict[str, dict] = {}
    if CACHE_JSON.exists():
        cache = json.loads(CACHE_JSON.read_text(encoding="utf-8")).get("日", {})

    missing = [d for d in dates if d.isoformat() not in cache]
    if missing:
        print(f"  advisor を {len(missing)} 日ぶん回す"
              f"（1日 0.2 秒ほど。{len(missing) * 0.2 / 60:.1f} 分の見込み）…")
        for index, date in enumerate(missing, start=1):
            forecast = forecast_from_date(radiation[date], date)
            lai = leaf_count_per_m2(date) * LEAF_AREA_PER_LEAF_M2
            # 潅水量はこちらで実績を使うので、上乗せ0・土壌は既定で回す。
            # 要るのは蒸散量（潜在）と糖だけ。
            result = advise(forecast, house=HOUSE, lai=lai,
                            leaching_fraction=0.0)
            cache[date.isoformat()] = {
                "蒸散": result.transpiration_l_per_m2,
                "糖": result.sugar_g_per_m2,
            }
            if index % 100 == 0:
                print(f"    {index}/{len(missing)} 日")
        CACHE_JSON.write_text(
            json.dumps({
                "説明": (
                    "advisor を日ごとに回した結果。[蒸散 L/m²/日, 糖 g/m²/日]。"
                    "重いのでキャッシュしている。作り直すときはこのファイルを消す。"
                ),
                "作成日": dt.date.today().isoformat(),
                "ハウス": HOUSE,
                "日": cache,
            }, ensure_ascii=False),
            encoding="utf-8")
        print(f"  書き出した: {CACHE_JSON}")
    else:
        print(f"  キャッシュを使う（{CACHE_JSON.name}）")
    return cache


# =============================================================================
# 2. 1作期ぶんを回す
# =============================================================================

def run_season(
    season: int,
    dates: list[dt.date],
    diary: dict[dt.date, dict[str, float]],
    fertilizer: dict[dt.date, float],
    cache: dict[str, dict],
    fresh_yield_kg_per_m2: float,
    start_n_g_per_m2: float,
    start_ec_solution: float,
) -> dict:
    """1作期ぶんの水・N・塩の収支をまとめて回す。"""
    # --- (a) 水収支。潅水は実績、蒸散はモデル ---
    water_days = [
        DailyWaterInput(
            date=date.isoformat(),
            irrigation_mm=diary[date]["irrigation_mm"],
            potential_transpiration_mm=cache[date.isoformat()]["蒸散"],
        )
        for date in dates
    ]
    water = simulate_water_balance(water_days)

    # --- (b) N需要。作期の合計を、糖の日別配分で割る ---
    #
    # ★絶対値は収量から、日ごとの形は光合成から。
    #   光合成モデルは素のままだと実収量の約1/3しか出ないので、
    #   絶対値には使わない（core/nitrogen.py の説明）。
    demand_total = demand_from_fresh_yield(
        fresh_yield_kg_per_m2).uptake_g_per_m2

    # ★★★ 配分の重みに「糖」を使ってはいけない（2026-10-01 に判明）★★★
    #
    # 当初は光合成モデルの糖を重みにする設計だった（core/nitrogen.py の
    # 「絶対値は収量から、日ごとの形は光合成から」）。これは**成り立たない**。
    #
    # 月ごとの糖の平均 [g/m²/日] と、そのときのLAI:
    #
    #     月      12     1     2     3     4     5     6
    #     糖    5.27  5.68  5.27  3.80  0.81 -0.61 -3.49
    #     LAI   2.89  4.02  5.08  6.24  7.19  7.19  6.40
    #
    # **5月・6月は平均が負**になる。LAIが大きいと暗呼吸（LAIに比例）が
    # 総光合成を上回るためだ。README の「実収量から逆算すると約3倍足りない」は
    # 作期合計の話で、**誤差が一定倍率ではなくLAI依存**だったということ。
    # 一定倍率なら分子と分母で打ち消し合うので配分には使えたが、
    # LAI依存なら打ち消し合わない。春の需要が消えてしまう。
    #
    # 【代わりに使うもの: 受光量】
    #
    #     重み = 日射 [MJ/m²] × ( 1 − exp( −k × LAI ) )
    #
    # 群落が実際に受け取った光のエネルギー。光利用効率（LUE）の考え方で、
    # 乾物生産がこれにほぼ比例することはよく確かめられている。
    # 使うのは**実測の日射**と**作業計画の葉枚数**だけなので、
    # 壊れている光合成モデルを通らない。
    #
    #     月      12     1     2     3     4     5     6
    #     日射   6.9   8.0   9.7  12.6  13.9  15.5  14.3
    #     受光率 0.87  0.94  0.97  0.99  0.99  0.99  0.99
    #     重み   6.0   7.5   9.4  12.4  13.8  15.4  14.1
    #
    # 春が冬の2.6倍になる。乾物生産の季節変化として素直な形だ。
    #
    # 【定植前はゼロ】
    # 作期は8月始まりだが定植は9月中旬（ISO週38）。作物がいないので需要ゼロ。
    # ★基肥は8月に入るので、収支そのものは8月から回す。需要だけを止める。
    weights = []
    for date in dates:
        if before_planting(date):
            weights.append(0.0)
            continue
        lai = leaf_count_per_m2(date) * LEAF_AREA_PER_LEAF_M2
        intercepted = 1.0 - math.exp(-EXTINCTION_COEFFICIENT_K * lai)
        weights.append(diary[date]["radiation_mj"] * intercepted)
    weight_total = sum(weights)
    if weight_total <= 0.0:
        raise ValueError(
            f"作期{season}: 受光量の重みの合計が0。日射の読み取りか、"
            f"定植日の判定（before_planting）を疑うこと。"
        )
    demand = [demand_total * value / weight_total for value in weights]

    # --- (c) 無機化。地力窒素を気温の重みで日割り ---
    #
    # 日平均気温は advisor が日付から推定した値を使えないので
    # （advise は日中平均しか返さない）、forecast から取り直す。
    temperatures = [
        forecast_from_date(diary[date]["radiation_mj"], date).mean_temp_c
        for date in dates
    ]
    weights = temperature_weights(temperatures)
    mineralization = [SOIL_N_SUPPLY_KG_PER_10A * w for w in weights]

    # --- (d) N収支 ---
    nitrogen_days = [
        DailyNitrogenInput(
            date=date.isoformat(),
            fertilizer_n_g_per_m2=fertilizer.get(date, 0.0),
            demand_n_g_per_m2=demand[index],
            irrigation_mm=water[index].irrigation_mm,
            drainage_mm=water[index].drainage_mm,
            mineralization_n_g_per_m2=mineralization[index],
            water_content=water[index].end_water_content,
        )
        for index, date in enumerate(dates)
    ]
    nitrogen = simulate(nitrogen_days, start_n_g_per_m2)

    # --- (e) 給液EC。その日の施肥NとN濃度から ---
    #
    # 基肥（8月の固形肥料）は液肥ではないので給液ECには乗せない。
    # 潅水がゼロの日も当然ゼロ。
    feed_ec = []
    for index, date in enumerate(dates):
        irrigation = water[index].irrigation_mm
        fertilizer_n = fertilizer.get(date, 0.0)
        if irrigation <= 0.0 or date.month == 8:
            feed_ec.append(0.0)
            continue
        # g/m² ÷ (L/m²) = g/L → ×1000 で mg/L
        n_mg_per_l = fertilizer_n / irrigation * 1000.0
        feed_ec.append(
            feed_ec_from_n_concentration(n_mg_per_l, FERTILIZER_N_FRACTION))

    # --- (f) 塩類収支 ---
    salt_days = [
        DailySaltInput(
            date=date.isoformat(),
            irrigation_mm=water[index].irrigation_mm,
            drainage_mm=water[index].drainage_mm,
            feed_ec_ds_per_m=feed_ec[index],
            n_uptake_g_per_m2=nitrogen[index].uptake_n_g_per_m2,
            water_content=water[index].end_water_content,
        )
        for index, date in enumerate(dates)
    ]
    salt = simulate_salt(salt_days, start_ec_solution, ROOT_ZONE_DEPTH_M)

    return {
        "作期": season,
        "日付": dates,
        "水": water,
        "N": nitrogen,
        "塩": salt,
        "給液EC": feed_ec,
        "需要N": demand,
        "無機化": mineralization,
        "気温": temperatures,
        "収量": fresh_yield_kg_per_m2,
        "需要N合計": demand_total,
    }


# =============================================================================
# 3. 月ごとにまとめる
# =============================================================================

def monthly(result: dict) -> dict[int, dict[str, float]]:
    """1作期ぶんの日別結果を、月ごとの平均・合計にまとめる。"""
    buckets: dict[int, dict[str, list]] = defaultdict(
        lambda: defaultdict(list))
    for index, date in enumerate(result["日付"]):
        bucket = buckets[date.month]
        bucket["施肥N"].append(result["N"][index].fertilizer_n_g_per_m2)
        bucket["無機化"].append(result["無機化"][index])
        bucket["需要N"].append(result["需要N"][index])
        bucket["吸収N"].append(result["N"][index].uptake_n_g_per_m2)
        bucket["流亡N"].append(result["N"][index].leaching_n_g_per_m2)
        bucket["根圏N"].append(result["N"][index].end_n_g_per_m2)
        bucket["潅水"].append(result["水"][index].irrigation_mm)
        bucket["蒸散"].append(result["水"][index].actual_transpiration_mm)
        bucket["潜在蒸散"].append(
            result["水"][index].potential_transpiration_mm)
        bucket["流亡水"].append(result["水"][index].drainage_mm)
        bucket["日射"].append(result["気温"][index])  # 置き場所の都合で気温
        bucket["含水率"].append(result["水"][index].end_water_content)
        bucket["給液EC"].append(result["給液EC"][index])
        bucket["根圏EC"].append(result["塩"][index].end_ec_solution)
        bucket["ECe"].append(result["塩"][index].end_ece)

    out: dict[int, dict[str, float]] = {}
    for month, bucket in buckets.items():
        days = len(bucket["潅水"])
        feed = [v for v in bucket["給液EC"] if v > 0.0]
        root_ec = sum(bucket["根圏EC"]) / days
        ece = sum(bucket["ECe"]) / days
        matric = FIELD_CAPACITY_POTENTIAL_J_PER_KG
        out[month] = {
            "日数": days,
            "施肥N": sum(bucket["施肥N"]),
            "無機化": sum(bucket["無機化"]),
            "需要N": sum(bucket["需要N"]),
            "吸収N": sum(bucket["吸収N"]),
            "流亡N": sum(bucket["流亡N"]),
            "根圏N": bucket["根圏N"][-1],
            "潅水": sum(bucket["潅水"]) / days,
            "蒸散": sum(bucket["蒸散"]) / days,
            "潜在蒸散": sum(bucket["潜在蒸散"]) / days,
            "流亡水": sum(bucket["流亡水"]) / days,
            "気温": sum(bucket["日射"]) / days,
            "含水率": sum(bucket["含水率"]) / days,
            "給液EC": sum(feed) / len(feed) if feed else 0.0,
            "根圏EC": root_ec,
            "ECe": ece,
            "浸透": osmotic_potential_j_per_kg(root_ec),
            "駆動力比": uptake_driving_force_ratio(matric, root_ec),
            "収量比": salinity_response(ece).yield_ratio,
            "糖度": salinity_response(ece).brix_change,
        }
    return out


# =============================================================================
# 4. 走らせて表示する
# =============================================================================

def main() -> None:
    theta_fc = water_content_from_potential(FIELD_CAPACITY_POTENTIAL_J_PER_KG)

    print("=" * 96)
    print(f"■ 0. 準備（{HOUSE}ハウスのみ。東は過去の記録が当てにならない）")
    print("=" * 96)

    diary = read_diary()
    fertilizer = load_fertilizer_daily()
    yields = load_yield()

    # 作期ごとに、記録のある日を並べる
    by_season: dict[int, list[dt.date]] = defaultdict(list)
    for date, record in diary.items():
        if record["radiation_mj"] <= 0.0:
            continue
        if date.month == 7:            # 栽培していない
            continue
        by_season[season_of(date)].append(date)
    for season in by_season:
        by_season[season].sort()

    usable = [s for s in SEASONS if s in by_season and s in yields]
    if not usable:
        raise ValueError(
            f"使える作期が無い。作業日誌の日付と data.js の start、"
            f"SEASONS={SEASONS} が合っているか確かめること。"
        )
    for season in usable:
        print(f"  作期{season}: {len(by_season[season])} 日  "
              f"収量 {yields[season]:.1f} kg/m²")

    all_dates = [d for s in usable for d in by_season[s]]
    radiation = {d: diary[d]["radiation_mj"] for d in all_dates}
    cache = load_or_build_cache(all_dates, radiation)

    # --- 作付け前のNとEC。土壌診断の実測から置く ---
    #
    # ★診断は6月採取＝その作期の終わり。だから「作期Nの終わりの値」を
    #   「作期N+1の始まりの値」として使う。診断年 = 作期年 + 1。
    print()
    print("=" * 96)
    print("■ 1. 作付け前のNとEC（土壌診断の実測から置く）")
    print("=" * 96)
    print("  ★診断は6月採取なので『その作期の終わり』の値。")
    print("    診断年 Y の値は、作期 Y の始まりの状態として使う")
    print("    （作期 Y−1 の終わり ＝ 作期 Y の始まり）。")
    print()
    print(f"{'診断年':>7}{'EC(1:5)':>9}{'土壌溶液EC':>12}{'無機態N':>9}"
          f"{'→ 作期':>9}")
    print(f"{'':>7}{'dS/m':>9}{'dS/m':>12}{'kg-N/10a':>9}")
    start_state: dict[int, tuple[float, float]] = {}
    for year, record in sorted(SOIL_DIAGNOSIS_CHUO.items()):
        ec5 = record["ec_1to5_ds_per_m"]
        solution = ec_1to5_to_soil_solution(ec5, theta_fc)
        nitrogen = mineral_n_to_g_per_m2(
            record["mineral_n_mg_per_100g"], SOIL_DIAGNOSIS_SAMPLE_DEPTH_M)
        start_state[year] = (nitrogen, solution)
        print(f"{year:>7}{ec5:>9.2f}{solution:>12.2f}{nitrogen:>9.1f}"
              f"{year:>9}")

    # 診断の無い作期は、ある年の平均を使う
    mean_n = sum(v[0] for v in start_state.values()) / len(start_state)
    mean_ec = sum(v[1] for v in start_state.values()) / len(start_state)
    print(f"\n  診断の無い作期には平均を使う: N {mean_n:.1f} kg/10a ／ "
          f"EC {mean_ec:.2f} dS/m")

    # --- 回す ---
    results = {}
    for season in usable:
        nitrogen0, ec0 = start_state.get(season, (mean_n, mean_ec))
        results[season] = run_season(
            season, by_season[season], diary, fertilizer, cache,
            yields[season], nitrogen0, ec0)

    # =========================================================================
    print()
    print("=" * 96)
    print("■ 2. 作期ごとの収支（N）")
    print("=" * 96)
    for season in usable:
        report = summarize(results[season]["N"])
        print(f"\n--- 作期{season}（収量 {yields[season]:.1f} kg/m²）---")
        print(report.describe())
        if abs(report.balance_residual) > 1e-6:
            raise ValueError(
                f"作期{season}: N収支が閉じていない"
                f"（残差 {report.balance_residual:.3e}）。計算が壊れている。"
            )

    # =========================================================================
    print()
    print("=" * 96)
    print("■ 3. 月ごとの姿（3作期の平均）")
    print("=" * 96)

    tables = {s: monthly(results[s]) for s in usable}
    averaged: dict[int, dict[str, float]] = {}
    for month in SEASON_MONTHS:
        rows = [tables[s][month] for s in usable if month in tables[s]]
        if not rows:
            continue
        averaged[month] = {
            key: sum(r[key] for r in rows) / len(rows) for key in rows[0]
        }

    print("  《N の側》  単位は kg-N/10a（= g-N/m²）、潅水は L/m²/日")
    print(f"{'月':>4}{'気温':>7}{'潅水':>7}{'流亡水':>8}{'施肥N':>8}"
          f"{'無機化':>8}{'需要N':>8}{'吸収N':>8}{'流亡N':>8}{'根圏N':>8}")
    for month, row in averaged.items():
        print(f"{month:>4}{row['気温']:>7.1f}{row['潅水']:>7.2f}"
              f"{row['流亡水']:>8.2f}{row['施肥N']:>8.2f}{row['無機化']:>8.2f}"
              f"{row['需要N']:>8.2f}{row['吸収N']:>8.2f}{row['流亡N']:>8.2f}"
              f"{row['根圏N']:>8.1f}")

    print()
    print()
    print("  《EC の側》  ★ここがデルフィの指導を確かめるところ")
    print("  ※根圏ECは現地の含水率のままの濃度。土が乾くと濃くなるので、")
    print("    含水率も並べてある。ECe は飽和まで薄めた基準（文献と比べる用）。")
    print(f"{'月':>4}{'含水率':>8}{'給液EC':>9}{'根圏EC':>9}{'ECe':>7}{'浸透':>8}"
          f"{'駆動力比':>10}{'収量比':>8}{'糖度':>8}")
    print(f"{'':>4}{'m³/m³':>8}{'dS/m':>9}{'dS/m':>9}{'dS/m':>7}{'J/kg':>8}"
          f"{'':>10}{'':>8}{'°Brix':>8}")
    for month, row in averaged.items():
        print(f"{month:>4}{row['含水率']:>8.3f}{row['給液EC']:>9.2f}"
              f"{row['根圏EC']:>9.2f}"
              f"{row['ECe']:>7.2f}{row['浸透']:>8.0f}{row['駆動力比']:>10.3f}"
              f"{row['収量比'] * 100:>7.1f}%{row['糖度']:>+8.2f}")

    # =========================================================================
    print()
    print("=" * 96)
    print("■ 4. 浸透ポテンシャルはマトリックポテンシャルの何倍か")
    print("=" * 96)
    print(f"  圃場容水量（pF 1.8）のマトリックポテンシャル "
          f"{FIELD_CAPACITY_POTENTIAL_J_PER_KG} J/kg")
    print()
    print(f"{'月':>4}{'根圏EC':>9}{'浸透 J/kg':>11}{'マトリックの何倍':>18}")
    for month, row in averaged.items():
        ratio = row["浸透"] / FIELD_CAPACITY_POTENTIAL_J_PER_KG
        print(f"{month:>4}{row['根圏EC']:>9.2f}{row['浸透']:>11.0f}"
              f"{ratio:>17.1f}倍")
    print()
    print("  ★十分に湿った土では、根が感じる水ポテンシャルはECで決まる。")
    print("    第9章の式9.20・9.22（含水率だけで吸水を決める）には")
    print("    この項が入っていない。")

    # =========================================================================
    print()
    print("=" * 96)
    print("■ 5. ECが吸水を削る量は、時期によって意味が違う")
    print("=" * 96)
    print("  駆動力の減り方（%）は同じでも、蒸散要求の大きさで効き方が変わる。")
    print("  蒸散要求が根の吸水能力の上限に近い時期ほど、削られた分が")
    print("  そのまま不足になる。")
    print()
    print("  判定は「ECで失う蒸散が、すでに届いていない量と比べてどれだけか」で付ける。")
    print("  すでに潜在蒸散に実蒸散が届いていない時期にさらに削るのが、いちばん痛い。")
    print()
    print(f"{'月':>4}{'潜在蒸散':>10}{'実蒸散':>9}{'届かない量':>12}"
          f"{'駆動力の減り':>13}{'ECで失う蒸散':>14}{'判定':>12}")
    print(f"{'':>4}{'mm/日':>10}{'mm/日':>9}{'mm/日':>12}{'%':>13}{'mm/日':>14}")
    for month, row in averaged.items():
        if month == 8:
            continue              # 定植前。作物がいない
        loss_ratio = 1.0 - row["駆動力比"]
        lost = row["潜在蒸散"] * loss_ratio
        gap = row["潜在蒸散"] - row["蒸散"]
        # すでに水が届いていない時期に、ECがさらに削る量が同じ桁なら重い
        if gap > 0.3 and lost >= gap * 0.5:
            verdict = "★いちばん痛い"
        elif lost >= 0.3:
            verdict = "効く"
        elif lost >= 0.15:
            verdict = "やや"
        else:
            verdict = "ほぼ無害"
        print(f"{month:>4}{row['潜在蒸散']:>10.2f}{row['蒸散']:>9.2f}"
              f"{gap:>12.2f}{loss_ratio * 100:>12.1f}%{lost:>14.2f}"
              f"{verdict:>12}")

    # =========================================================================
    print()
    print("=" * 96)
    print("■ 6. 作終わりのECを土壌診断の実測と突き合わせる")
    print("=" * 96)
    print("  ★これがこのモデルの答え合わせ。診断年 = 作期年 + 1。")
    print()
    print(f"{'作期':>6}{'診断年':>8}{'モデルの作終わりEC':>20}"
          f"{'実測から換算':>14}{'比':>8}")
    print(f"{'':>6}{'':>8}{'dS/m（土壌溶液）':>20}{'dS/m':>14}")
    for season in usable:
        model_ec = results[season]["塩"][-1].end_ec_solution
        year = season + 1
        if year in SOIL_DIAGNOSIS_CHUO:
            measured = ec_1to5_to_soil_solution(
                SOIL_DIAGNOSIS_CHUO[year]["ec_1to5_ds_per_m"], theta_fc)
            print(f"{season:>6}{year:>8}{model_ec:>20.2f}{measured:>14.2f}"
                  f"{model_ec / measured:>8.2f}")
        else:
            print(f"{season:>6}{year:>8}{model_ec:>20.2f}"
                  f"{'診断なし':>14}{'—':>8}")

    print()
    print("  定常状態の濃縮倍率（根圏EC ÷ 給液EC = 1 ÷ 流亡率）:")
    print("  ★流亡率がゼロに近いと発散する。塩が出ていく先が無いので、")
    print("    そもそも落ち着き先が存在しないという意味だ。")
    for month in SEASON_MONTHS:
        if month not in averaged:
            continue
        row = averaged[month]
        if row["潅水"] <= 0.0:
            continue
        fraction = row["流亡水"] / row["潅水"]
        if fraction < 0.02:
            print(f"    {month:2d}月  流亡率 {fraction * 100:4.1f}% → "
                  f"落ち着き先なし（塩は溜まり続ける）")
            continue
        factor = steady_state_concentration_factor(fraction)
        print(f"    {month:2d}月  流亡率 {fraction * 100:4.1f}% → "
              f"濃縮 {factor:5.1f} 倍（給液 {row['給液EC']:.2f} → 上限 "
              f"{row['給液EC'] * factor:.2f} dS/m）")

    # =========================================================================
    print()
    print("=" * 96)
    print("■ 7. まとめ — デルフィの指導は妥当か")
    print("=" * 96)
    winter = [averaged[m] for m in (12, 1, 2) if m in averaged]
    spring = [averaged[m] for m in (4, 5, 6) if m in averaged]
    if winter and spring:
        def mean(rows, key):
            return sum(r[key] for r in rows) / len(rows)

        print(f"{'':<18}{'冬（12〜2月）':>16}{'春（4〜6月）':>16}")
        for label, key, form in (
            ("給液EC dS/m", "給液EC", "{:.2f}"),
            ("根圏EC dS/m", "根圏EC", "{:.2f}"),
            ("ECe dS/m", "ECe", "{:.2f}"),
            ("潜在蒸散 mm/日", "潜在蒸散", "{:.2f}"),
            ("駆動力の減り %", "駆動力比", "{:.1f}"),
            ("ECで失う蒸散 mm/日", None, "{:.2f}"),
            ("施肥N kg/10a/月", "施肥N", "{:.2f}"),
            ("無機化 kg/10a/月", "無機化", "{:.2f}"),
            ("需要N kg/10a/月", "需要N", "{:.2f}"),
        ):
            if key == "駆動力比":
                values = [(1 - mean(rows, key)) * 100 for rows in (winter, spring)]
            elif key is None:
                values = [
                    mean(rows, "潜在蒸散") * (1 - mean(rows, "駆動力比"))
                    for rows in (winter, spring)
                ]
            else:
                values = [mean(rows, key) for rows in (winter, spring)]
            print(f"{label:<18}{form.format(values[0]):>16}"
                  f"{form.format(values[1]):>16}")

    print()
    print(f"  塩害のしきい値（ECe {SALINITY_YIELD_THRESHOLD_ECE} dS/m）を"
          f"超えた月:")
    over = [m for m, r in averaged.items()
            if r["ECe"] > SALINITY_YIELD_THRESHOLD_ECE]
    print(f"    {over if over else 'なし'}")
    print(f"  土壌診断の理想値（EC(1:5) {SOIL_DIAGNOSIS_EC_TARGET_1TO5} 未満）"
          f"に相当する土壌溶液EC: "
          f"{ec_1to5_to_soil_solution(SOIL_DIAGNOSIS_EC_TARGET_1TO5, theta_fc):.2f} dS/m")


if __name__ == "__main__":
    main()
