"""施肥Nの実績と、収量から逆算したN需要を突き合わせる。

【このツールが答える問い】
「乾物 → 器官別N濃度 → 吸収N」という組み立ては成り立っているか。

成り立っているかどうかは、**逆算した地力窒素がありうる値に収まるか**で
判断する。地力窒素（土が作期を通して出すN）は土耕では測りようがないので、
これを唯一の未知数にして実績から解く。

    地力窒素 = 吸収N需要（収量から逆算） − 施肥N × 利用率

8作期 × 2ハウス = 16点あるので、16通りの答えが出る。値がそろっていて、
かつ施設土壌としてありうる範囲なら、前提は成り立っている。
ばらばらに散ったり、負になったり、ありえない大きさになったら、
N濃度・分配率・収量のどこかが間違っている。

【データの出どころ】
    施肥N   data/施肥実績.json（tools/施肥実績を作り直す.py が作る）
    収量     Dropbox/作業記録/data.js（収量ビューアが読んでいるもの）

★収量をこのリポジトリに書き出さない。config 第8-3節と同じ扱いで、
  実収量は手元のファイルから読むだけにする。

【使い方】
    py tools/N収支を検証.py
"""
from __future__ import annotations

import datetime as dt
import io
import json
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (                                      # noqa: E402
    FERTILIZER_N_EFFICIENCY,
    FERTILIZER_N_EFFICIENCY_SENSITIVITY,
    N_PER_TONNE_FRUIT_LITERATURE_RANGE,
    SOIL_N_SUPPLY_KG_PER_10A,
    SOIL_N_SUPPLY_PLAUSIBLE_RANGE,
    YIELD_CONVERSION_SENSITIVITY,
)

#: 土耕の点滴潅水で見込まれる施肥Nの利用率の範囲（文献）。
#: config.FERTILIZER_N_EFFICIENCY_SENSITIVITY の端から取る。
EFFICIENCY_RANGE = (min(FERTILIZER_N_EFFICIENCY_SENSITIVITY),
                    max(FERTILIZER_N_EFFICIENCY_SENSITIVITY))
from core.nitrogen import (                                # noqa: E402
    NitrogenSettings,
    demand_from_fresh_yield,
    solve_soil_n_supply,
)

SHIHI_JSON = Path(__file__).resolve().parent.parent / "data" / "施肥実績.json"
YIELD_JS = Path(r"c:\Users\kimij\Dropbox\作業記録\data.js")

HOUSES = ("中央", "東")


# =============================================================================
# 1. データを読む
# =============================================================================

def load_fertilizer() -> dict:
    """施肥実績.json を読む。"""
    if not SHIHI_JSON.exists():
        raise FileNotFoundError(
            f"{SHIHI_JSON} が無い。先に "
            f"`py tools/施肥実績を作り直す.py` を走らせること。"
        )
    return json.loads(SHIHI_JSON.read_text(encoding="utf-8"))


def load_yield() -> dict[tuple[int, str], float]:
    """data.js から「(作期, ハウス) → 収量 [kg/m²]」を読む。

    data.js は `window.DATA = {...};` という形の JavaScript なので、
    最初の { から末尾のセミコロンまでを切り出して JSON として読む。
    """
    if not YIELD_JS.exists():
        raise FileNotFoundError(
            f"{YIELD_JS} が無い。Dropbox の同期が済んでいるか確かめること。"
            f"このファイルは shuryo_export.py が作る。"
        )
    text = YIELD_JS.read_text(encoding="utf-8")
    start = text.find("{")
    if start < 0:
        raise ValueError(
            f"{YIELD_JS} に JSON の始まり '{{' が見つからない。"
            f"ファイルの作りが変わった可能性がある。"
        )
    document = json.loads(text[start:].rstrip().rstrip(";"))

    out = {}
    for season in document["seasons"]:
        year = int(season["start"])          # 作期は8月始まりの年
        for house, record in season["houses"].items():
            if house not in HOUSES:
                continue
            area_m2 = float(record["area"]) * 1000.0   # 1.74 → 1740 m²
            kg = sum(float(v) for v in record["kg"])
            if kg > 0:
                out[(year, house)] = kg / area_m2
    return out


# =============================================================================
# 2. 1作期ぶんの突き合わせ
# =============================================================================

def inspect_season(
    fresh_yield_kg_per_m2: float,
    fertilizer_n_kg_per_10a: float,
    settings: NitrogenSettings,
    efficiency: float,
) -> dict[str, float]:
    """1作期・1ハウスぶんの数字をまとめる。"""
    demand = demand_from_fresh_yield(fresh_yield_kg_per_m2, settings)
    soil = solve_soil_n_supply(
        demand.uptake_kg_per_10a, fertilizer_n_kg_per_10a, efficiency)
    return {
        "収量": fresh_yield_kg_per_m2,
        "需要N": demand.uptake_kg_per_10a,
        "持ち出しN": demand.removed_g_per_m2,
        "施肥N": fertilizer_n_kg_per_10a,
        "地力窒素": soil,
        # 施肥N 1kg あたり何kgの果実が穫れたか。
        # 収量 [kg/m²] × 1000 = [kg/10a] なので ×1000 する。
        "果実kg/施肥N1kg": (
            fresh_yield_kg_per_m2 * 1000.0 / fertilizer_n_kg_per_10a
            if fertilizer_n_kg_per_10a > 0 else float("nan")
        ),
    }


# =============================================================================
# 3. 表示
# =============================================================================

def main() -> None:
    fertilizer = load_fertilizer()
    yields = load_yield()
    settings = NitrogenSettings()

    print("=" * 74)
    print("■ 1. 前提（config 第10節。すべて文献の中央値で、実測ではない）")
    print("=" * 74)
    print(settings.describe())
    low, high = N_PER_TONNE_FRUIT_LITERATURE_RANGE
    if not low <= settings.n_per_tonne_fruit_kg <= high:
        print(f"  ★文献の幅 {low}〜{high} の外にある。前提を見直すこと。")
    elif settings.n_per_tonne_fruit_kg < low + (high - low) * 0.25:
        print(f"  ※文献の幅 {low}〜{high} の下端寄り。"
              f"N需要を低めに見積もる向きで、地力窒素の逆算値も低めに出る。")

    totals = fertilizer["作期の合計"]
    seasons = sorted(int(s) for s in totals["中央"])

    print()
    print("=" * 74)
    print(f"■ 2. 作期ごとの突き合わせ（利用率 {FERTILIZER_N_EFFICIENCY:.2f} のとき）")
    print("=" * 74)
    print("  地力窒素 = 需要N − 施肥N × 利用率")
    print()
    header = (f"{'作期':>5}{'ハウス':>5}{'収量':>7}{'需要N':>8}{'施肥N':>8}"
              f"{'基肥':>7}{'追肥':>7}{'地力窒素':>9}{'果実/N':>8}")
    print(header)
    print(f"{'':>5}{'':>5}{'kg/m²':>7}" + "".join(f"{'kg/10a':>8}" for _ in range(2))
          + f"{'kg/10a':>7}{'kg/10a':>7}{'kg/10a':>9}{'kg/kg':>8}")

    # ★記録が不完全な作期は印を付け、逆算からは外す。
    #   施肥Nが実際より少なく出るので、地力窒素を過大に見積もってしまう。
    incomplete = {
        (row["作期"], row["ハウス"]) for row in fertilizer.get("★記録が不完全な作期", [])
    }

    rows = []
    excluded = []
    for season in seasons:
        for house in HOUSES:
            bucket = totals[house].get(str(season))
            harvest = yields.get((season, house))
            if bucket is None or harvest is None:
                print(f"{season:>5}{house:>5}  収量または施肥の記録が無いので飛ばす")
                continue
            result = inspect_season(
                harvest, bucket["N"], settings, FERTILIZER_N_EFFICIENCY)
            flagged = (season, house) in incomplete
            (excluded if flagged else rows).append((season, house, bucket, result))
            print(f"{season:>5}{house:>5}{result['収量']:>7.1f}"
                  f"{result['需要N']:>8.1f}{result['施肥N']:>8.1f}"
                  f"{bucket['基肥N']:>7.1f}{bucket['追肥N']:>7.1f}"
                  f"{result['地力窒素']:>9.1f}{result['果実kg/施肥N1kg']:>8.0f}"
                  + ("  ★記録が不完全" if flagged else ""))

    if not rows:
        raise ValueError(
            "突き合わせられる作期が1つも無かった。"
            "施肥実績.json の作期と data.js の start が合っているか確かめること。"
        )

    if excluded:
        print()
        print("  ★記録が不完全な作期を逆算から外した:")
        for row in fertilizer.get("★記録が不完全な作期", []):
            print(f"    作期{row['作期']} {row['ハウス']}ハウス — {row['理由']}")
        print(f"  使うのは {len(rows)} 点（外したのは {len(excluded)} 点）。")

    soil_values = [r[3]["地力窒素"] for r in rows]
    plausible_low, plausible_high = SOIL_N_SUPPLY_PLAUSIBLE_RANGE
    print()
    print(f"  逆算した地力窒素: {min(soil_values):.1f} 〜 {max(soil_values):.1f} "
          f"（中央値 {sorted(soil_values)[len(soil_values) // 2]:.1f}）kg-N/10a")
    print(f"  施設土壌でありうる範囲として置いた値: "
          f"{plausible_low}〜{plausible_high} kg-N/10a")
    outside = [v for v in soil_values if not plausible_low <= v <= plausible_high]
    if outside:
        print(f"  ★16点のうち {len(outside)} 点がこの範囲の外にある。"
              f"下の第4節を読むこと。")

    # --- 3. 施肥Nを増やすと収量は増えたのか ---
    print()
    print("=" * 74)
    print("■ 3. 施肥Nと収量の関係（施肥Nは足りていたのか）")
    print("=" * 74)
    print("  施肥Nの少ない順に並べる。収量が頭打ちなら、その水準では")
    print("  Nは律速になっていない（＝それ以上入れても収量は増えない）。")
    print()
    print(f"{'施肥N':>8}{'作期':>7}{'ハウス':>6}{'収量':>8}")
    print(f"{'kg/10a':>8}{'':>7}{'':>6}{'kg/m²':>8}")
    for season, house, bucket, result in sorted(rows, key=lambda r: r[3]["施肥N"]):
        print(f"{result['施肥N']:>8.1f}{season:>7}{house:>6}{result['収量']:>8.1f}")

    # --- 4. 感度 ---
    print()
    print("=" * 74)
    print("■ 4. 前提を振ったときの地力窒素（直近4作期・中央ハウスの平均）")
    print("=" * 74)
    print("  実測していない値を1つずつ動かし、逆算される地力窒素がどう動くかを見る。")
    print("  幅が大きいものほど、測る価値が高い。")
    print()

    recent = [r for r in rows if r[0] >= 2022 and r[1] == "中央"]
    if not recent:
        raise ValueError("直近の中央ハウスのデータが取れなかった。")

    def average_soil(custom: NitrogenSettings, efficiency: float) -> float:
        values = [
            solve_soil_n_supply(
                demand_from_fresh_yield(r[3]["収量"], custom).uptake_kg_per_10a,
                r[3]["施肥N"], efficiency)
            for r in recent
        ]
        return sum(values) / len(values)

    base = average_soil(settings, FERTILIZER_N_EFFICIENCY)
    print(f"{'動かすもの':<22}{'値':>9}{'地力窒素':>11}{'既定との差':>11}")
    print(f"{'（既定のまま）':<22}{'—':>9}{base:>11.1f}{0.0:>11.1f}")

    cases: list[tuple[str, str, NitrogenSettings, float]] = []
    for value in YIELD_CONVERSION_SENSITIVITY["FRUIT_DRY_MATTER_CONTENT"]:
        cases.append(("果実乾物率", f"{value * 100:.1f}%",
                      replace(settings, fruit_dry_matter_content=value),
                      FERTILIZER_N_EFFICIENCY))
    for value in YIELD_CONVERSION_SENSITIVITY["FRUIT_ALLOCATION"]:
        cases.append(("果実分配率", f"{value:.2f}",
                      replace(settings, fruit_allocation=value),
                      FERTILIZER_N_EFFICIENCY))
    for value in (0.035, 0.040, 0.050):
        cases.append(("葉のN濃度", f"{value * 100:.1f}%",
                      replace(settings, leaf_n=value), FERTILIZER_N_EFFICIENCY))
    for value in (0.015, 0.018, 0.022):
        cases.append(("果実のN濃度", f"{value * 100:.1f}%",
                      replace(settings, fruit_n=value), FERTILIZER_N_EFFICIENCY))
    for value in FERTILIZER_N_EFFICIENCY_SENSITIVITY + (1.0,):
        cases.append(("施肥Nの利用率", f"{value:.2f}", settings, value))

    previous_label = None
    for label, shown, custom, efficiency in cases:
        if label != previous_label:
            print(f"{'':<22}{'':>9}{'':>11}{'':>11}")
            previous_label = label
        value = average_soil(custom, efficiency)
        print(f"{label:<22}{shown:>9}{value:>11.1f}{value - base:>+11.1f}")

    # --- 5. 向きを変えて解く: 地力窒素を一定と置いて利用率を求める ---
    print()
    print("=" * 74)
    print("■ 5. 向きを変えて解く（地力窒素を一定と置き、利用率を求める）")
    print("=" * 74)
    print("  第2節は利用率を固定して地力窒素を解いた。すると施肥Nを減らした年ほど")
    print("  地力窒素が大きく出る。土が変わったのではなく、ただの引き算の結果だ。")
    print()
    print("  同じ土なら地力窒素は年によらずほぼ一定のはず。そう置くと、")
    print("  動かざるをえないのは利用率のほうになる。")
    print()
    print("      利用率 = ( 需要N − 地力窒素 ) ÷ 施肥N")
    print()

    assumed_soils = (15.0, 20.0, 25.0)
    print(f"{'作期':>5}{'ハウス':>5}{'施肥N':>8}"
          + "".join(f"{'地力' + str(int(s)):>9}" for s in assumed_soils))
    print(f"{'':>5}{'':>5}{'kg/10a':>8}"
          + "".join(f"{'の利用率':>9}" for _ in assumed_soils))
    for season, house, bucket, result in rows:
        line = f"{season:>5}{house:>5}{result['施肥N']:>8.1f}"
        for soil in assumed_soils:
            efficiency = (result["需要N"] - soil) / result["施肥N"]
            # 1 を超えたら「施肥だけでは需要に届かない」という意味。
            # 物理的にありえないので印を付ける。
            mark = "★" if efficiency > 1.0 else " "
            line += f"{efficiency:>8.2f}{mark}"
        print(line)

    print()
    print("  ★ が付いた行は利用率が1を超えている（＝施肥Nを全部吸っても足りない）。")
    print("    その地力窒素の仮定では説明がつかないということだ。")
    print()
    print("  中央ハウスだけで、施肥Nの多い年と少ない年を比べる:")
    central = sorted((r for r in rows if r[1] == "中央"),
                     key=lambda r: r[3]["施肥N"])
    for soil in assumed_soils:
        least, most = central[0], central[-1]
        eff_least = (least[3]["需要N"] - soil) / least[3]["施肥N"]
        eff_most = (most[3]["需要N"] - soil) / most[3]["施肥N"]
        print(f"    地力 {soil:.0f} のとき: "
              f"施肥 {most[3]['施肥N']:.1f}（作期{most[0]}）→ 利用率 {eff_most:.2f} ／ "
              f"施肥 {least[3]['施肥N']:.1f}（作期{least[0]}）→ 利用率 {eff_least:.2f}")
    print()
    print("  施肥Nを減らすほど利用率が上がる。これは実績のN効率")
    print("  （果実kg ÷ 施肥N kg）が第2節で 417 → 620 へ上がったことと同じ現象で、")
    print("  **多く入れた年は、増やした分が吸われずに抜けていた**と読める。")

    # --- 5-2. 利用率は1を超えられない。そこから地力窒素の下限が決まる ---
    print()
    print("-" * 74)
    print("  5-2. 地力窒素の下限（まず硬い制約を1つ置く）")
    print("-" * 74)
    print("  利用率は物理的に1を超えられない。施肥したNより多くは吸えないからだ。")
    print("  したがって、どの作期でも")
    print()
    print("      地力窒素 ≧ 需要N − 施肥N")
    print()
    print("  が成り立たなければならない。いちばん厳しい作期がこの下限を決める。")
    print()

    binding = max(rows, key=lambda r: r[3]["需要N"] - r[3]["施肥N"])
    floor = binding[3]["需要N"] - binding[3]["施肥N"]
    print(f"{'作期':>5}{'ハウス':>5}{'需要N':>8}{'施肥N':>8}{'差（下限）':>11}")
    for season, house, bucket, result in sorted(
            rows, key=lambda r: -(r[3]["需要N"] - r[3]["施肥N"]))[:5]:
        gap = result["需要N"] - result["施肥N"]
        mark = "  ← これが効く" if (season, house) == binding[:2] else ""
        print(f"{season:>5}{house:>5}{result['需要N']:>8.1f}"
              f"{result['施肥N']:>8.1f}{gap:>11.1f}{mark}")

    print()
    print(f"  → 地力窒素は少なくとも {floor:.1f} kg-N/10a/作 ないと、")
    print(f"     作期{binding[0]}の{binding[1]}ハウス"
          f"（施肥N {binding[3]['施肥N']:.1f}・収量 {binding[3]['収量']:.1f} kg/m²）"
          f"が説明できない。")
    print()

    # その下限を採ったときの利用率。1に近すぎるなら、下限は実際の値ではない。
    efficiencies = [
        (r[0], r[1], (r[3]["需要N"] - floor) / r[3]["施肥N"]) for r in rows
    ]
    print(f"  ただし、地力窒素をこの下限 {floor:.1f} に置くと利用率が "
          f"{min(e[2] for e in efficiencies):.2f}〜"
          f"{max(e[2] for e in efficiencies):.2f} になる。")
    print(f"  土耕の点滴で見込まれる {EFFICIENCY_RANGE[0]:.1f}〜"
          f"{EFFICIENCY_RANGE[1]:.1f} より高い。")
    print("  **下限は下限であって、実際の値ではない。**もっと上にあるはずだ。")

    # --- 5-3. 利用率の範囲と組み合わせて絞る ---
    print()
    print("-" * 74)
    print("  5-3. 利用率の範囲と組み合わせて絞る（★ここが答え）")
    print("-" * 74)
    print("  16点から解けるのは1本の式しかない。")
    print()
    print("      地力窒素 + 施肥N × 利用率 = 需要N")
    print()
    print("  未知数は2つ（地力窒素・利用率）なので、**組み合わせだけが決まる**。")
    print("  そこで利用率に文献の範囲を当て、地力窒素の範囲を出す。")
    print()
    print("  使うのは**中央ハウスの直近4作期**。記録が完全で、かつ施肥水準が")
    print("  いまの運用（約48 kg-N/10a）に近いからだ。")
    print()

    recent = [r for r in rows if r[0] >= 2022 and r[1] == "中央"]
    mean_demand = sum(r[3]["需要N"] for r in recent) / len(recent)
    mean_fertilizer = sum(r[3]["施肥N"] for r in recent) / len(recent)
    print(f"    平均の需要N   {mean_demand:.1f} kg-N/10a")
    print(f"    平均の施肥N   {mean_fertilizer:.1f} kg-N/10a")
    print()
    print(f"{'利用率':>8}{'地力窒素':>11}{'判定':>8}")
    for efficiency in (0.50, 0.60, 0.70, 0.80, 0.90, 1.00):
        soil = mean_demand - mean_fertilizer * efficiency
        inside = EFFICIENCY_RANGE[0] <= efficiency <= EFFICIENCY_RANGE[1]
        print(f"{efficiency:>8.2f}{soil:>11.1f}"
              + f"{'  ← 文献の範囲' if inside else '':>8}")

    soil_low = mean_demand - mean_fertilizer * EFFICIENCY_RANGE[1]
    soil_high = mean_demand - mean_fertilizer * EFFICIENCY_RANGE[0]
    soil_mid = mean_demand - mean_fertilizer * FERTILIZER_N_EFFICIENCY
    print()
    print(f"  → 利用率が {EFFICIENCY_RANGE[0]:.1f}〜{EFFICIENCY_RANGE[1]:.1f} なら、")
    print(f"     地力窒素は {soil_low:.1f}〜{soil_high:.1f} kg-N/10a/作。")
    print(f"     中ほど（利用率 {FERTILIZER_N_EFFICIENCY:.2f}）で "
          f"{soil_mid:.1f} kg-N/10a/作。")
    print()
    print(f"  下限 {floor:.1f}（利用率1のとき）よりも上にあり、矛盾しない。")
    print(f"  config.SOIL_N_SUPPLY_KG_PER_10A は {SOIL_N_SUPPLY_KG_PER_10A:.1f}。"
          + ("この逆算と合っている。"
             if abs(soil_mid - SOIL_N_SUPPLY_KG_PER_10A) < 1.0
             else f"逆算値 {soil_mid:.1f} と食い違うので見直すこと。"))
    print()
    print("  ★地力窒素と利用率は対で使うこと。片方だけ動かすと合わなくなる。")

    print()
    print("=" * 74)
    print("■ 6. 読み方")
    print("=" * 74)
    print("  ・地力窒素が負になるなら、施肥だけで需要を超えている（入れすぎ）")
    print("  ・地力窒素が 40 を超えるなら、需要の見積りが高すぎるか、")
    print("    この土が本当に多くのNを出している（堆肥を長年入れた施設土壌では")
    print("    ありうる）。土壌診断の硝酸態N・アンモニア態Nを1回読めば決着する")
    print("  ・第3節で収量が頭打ちなら、その水準の施肥Nは足りている")


if __name__ == "__main__":
    main()
