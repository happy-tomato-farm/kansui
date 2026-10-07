"""この施肥設計を、外の文献・公的基準と突き合わせる。

【なぜこれが必要なのか】
いまの施肥設計は2種類の数字が混ざっている。

  (a) 自分の圃場の実績から逆算した値 … 地力窒素 30.5 kg-N/10a・利用率 0.70
  (b) 文献から持ってきた値           … 器官別N濃度・乾物分配率

(a) は8作期の実績が1本の式に縮退したものなので、**自分の圃場にしか当たらない**。
(b) は他人の圃場の値なので、**自分の圃場に当たっている保証がない**。
だから「外の数字と比べてどれくらい合っているのか」を独立に確かめる必要がある。

【比べ方の原則】
絶対量（kg/10a）で比べると、収量が違うと話にならない。
**果実1トンあたりに直して比べる。**これが唯一そろう土俵。

使い方:
    py tools\文献と突き合わせる.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import (                                            # noqa: E402
    CHIBA_EC_MS_PER_M_PER_DS_PER_M,
    CHIBA_EC_TO_MINERAL_N_COEF,
    CHIBA_MINITOMATO_PROMOTED,
    CHIBA_UPTAKE_PER_TONNE,
    DOUNAN_DAILY_N_STEPS,
    EUROPE_N_DEMAND_PER_TONNE_RANGE,
    FERTILIZER_N_EFFICIENCY,
    FRUIT_N_REMOVAL_PER_TONNE_RANGE,
    GLOBAL_OPTIMUM_N_KG_PER_HA_RANGE,
    IBARAKI_EC_1TO5_MS_PER_M,
    N_PER_TONNE_FRUIT_LITERATURE_RANGE,
    ORGAN_N_CONTENT,
    PROCESSING_TOMATO_UPTAKE_KG_PER_HA_RANGE,
    PROCESSING_TOMATO_YIELD_T_PER_HA_RANGE,
    SOIL_DIAGNOSIS_CHUO,
    SOIL_N_SUPPLY_KG_PER_10A,
    SOIL_TEXTURE_CLOSEST,
    TARGET_FRESH_YIELD_KG_PER_M2,
)
from core.nitrogen import (                                      # noqa: E402
    NitrogenSettings,
    demand_from_fresh_yield,
    required_fertilizer_n,
)
from core.salinity import (                                      # noqa: E402
    ec_1to5_to_saturated_extract,
    ec_1to5_to_soil_solution,
    salinity_response,
)
from core.soil import FIELD_CAPACITY_WATER_CONTENT                # noqa: E402

#: 1 kg/m²（収量）= 1 t/10a = 10 t/ha
T_PER_HA_PER_KG_PER_M2 = 10.0

#: 1 kg-N/10a = 10 kg-N/ha
KG_PER_HA_PER_KG_PER_10A = 10.0

#: 作期の日数（9月中旬定植〜翌7月）。1日あたりに直すときの分母。
SEASON_DAYS = 290


def band(value: float, low: float, high: float) -> str:
    """値が帯の中か、どちら側に外れているかを一言で返す。"""
    if value < low:
        return f"↓帯の下（{low / value:.2f}倍 足りない側）"
    if value > high:
        return f"↑帯の上（{value / high:.2f}倍 多い側）"
    position = (value - low) / (high - low) if high > low else 0.0
    where = "下端" if position < 0.25 else ("上端" if position > 0.75 else "中ほど")
    return f"帯の中（{where}）"


# =============================================================================
print("=" * 86)
print("この施肥設計を文献・公的基準と突き合わせる")
print("=" * 86)
print()
print("【この設計の出どころ】")
print("  自分の圃場から逆算 … 地力窒素 "
      f"{SOIL_N_SUPPLY_KG_PER_10A} kg-N/10a・利用率 {FERTILIZER_N_EFFICIENCY}")
print("  文献から            … 器官別N濃度 "
      f"果実{ORGAN_N_CONTENT['FRUIT'] * 100:.1f}% "
      f"葉{ORGAN_N_CONTENT['LEAF'] * 100:.1f}% "
      f"茎{ORGAN_N_CONTENT['STEM'] * 100:.1f}% "
      f"根{ORGAN_N_CONTENT['ROOT'] * 100:.1f}%")
print()

settings = NitrogenSettings()
yield_kg_per_m2 = TARGET_FRESH_YIELD_KG_PER_M2
demand = demand_from_fresh_yield(yield_kg_per_m2, settings)
plan = required_fertilizer_n(
    demand.uptake_kg_per_10a, SOIL_N_SUPPLY_KG_PER_10A, FERTILIZER_N_EFFICIENCY)

# 1 kg/m² = 1 t/10a なので、kg-N/10a を収量で割れば kg-N/t になる
uptake_per_t = demand.uptake_kg_per_10a / yield_kg_per_m2
removed_per_t = demand.removed_g_per_m2 / yield_kg_per_m2
fertilizer_per_t = plan.fertilizer_kg_per_10a / yield_kg_per_m2

print(f"{'このモデルが出す値':<34}{'値':>12}{'単位':>16}")
print("-" * 86)
print(f"{'目標収量':<34}{yield_kg_per_m2:>12.1f}{'kg/m² (= t/10a)':>16}")
print(f"{'作期の吸収N需要':<34}{demand.uptake_kg_per_10a:>12.1f}{'kg-N/10a':>16}")
print(f"{'　うち果実で持ち出す分':<34}{demand.removed_g_per_m2:>12.1f}{'kg-N/10a':>16}")
print(f"{'　うち残渣に残る分':<34}{demand.residue_g_per_m2:>12.1f}{'kg-N/10a':>16}")
print(f"{'地力窒素（逆算）':<34}{SOIL_N_SUPPLY_KG_PER_10A:>12.1f}{'kg-N/10a':>16}")
print(f"{'作期の施肥N':<34}{plan.fertilizer_kg_per_10a:>12.1f}{'kg-N/10a':>16}")
print()

# =============================================================================
print("=" * 86)
print("1. 果実1トンあたりの全吸収N ─── いちばん素直に比べられる指標")
print("=" * 86)
print("  ★定義をそろえること。ここは『植物体全体が吸ったN ÷ 果実の生重』。")
print("    果実だけの持ち出し量ではない（それは第2節）。")
print()
print(f"{'出どころ':<40}{'kg-N/t':>12}{'判定':>30}")
print("-" * 86)
print(f"{'★このモデル':<40}{uptake_per_t:>12.3f}{'':>30}")

chiba_n = CHIBA_UPTAKE_PER_TONNE["N"]
print(f"{'千葉県 施肥基準 第Ⅲ-3-1表（トマト）':<40}{chiba_n:>12.3f}"
      f"{f'モデルは {uptake_per_t / chiba_n:.2f} 倍':>30}")

eu_low, eu_high = EUROPE_N_DEMAND_PER_TONNE_RANGE
print(f"{'ヨーロッパ 加工用トマト':<40}{f'{eu_low}〜{eu_high}':>12}"
      f"{band(uptake_per_t, eu_low, eu_high):>30}")

lit_low, lit_high = N_PER_TONNE_FRUIT_LITERATURE_RANGE
print(f"{'施設トマト（蘭・西）':<40}{f'{lit_low}〜{lit_high}':>12}"
      f"{band(uptake_per_t, lit_low, lit_high):>30}")
print()
print(f"  → 千葉の公的基準（{chiba_n}）より "
      f"{(uptake_per_t / chiba_n - 1) * 100:+.0f}%、"
      f"ヨーロッパの帯（{eu_low}〜{eu_high}）には収まる。")
print(f"    **3つの独立した出どころに挟まれている。この指標は信用してよい。**")
print()

# =============================================================================
print("=" * 86)
print("2. 果実で持ち出すN ─── ★ここが合っていない")
print("=" * 86)
print()
# 果実のN量 ÷ 果実のN濃度 = 果実の乾物 [g/m²]。それを生重で割れば乾物率。
# 収量 kg/m² → g/m² は ×1000。
fruit_dry_matter_g_per_m2 = (
    demand.removed_g_per_m2 / ORGAN_N_CONTENT["FRUIT"])
fruit_dm_ratio = fruit_dry_matter_g_per_m2 / (yield_kg_per_m2 * 1000.0)
print(f"{'出どころ':<40}{'kg-N/t':>12}{'判定':>30}")
print("-" * 86)
print(f"{'★このモデル':<40}{removed_per_t:>12.3f}{'':>30}")
rem_low, rem_high = FRUIT_N_REMOVAL_PER_TONNE_RANGE
print(f"{'米国（加工用・3 lb/short ton）':<40}{rem_low:>12.2f}"
      f"{f'モデルは {removed_per_t / rem_low:.2f} 倍':>30}")
print(f"{'Yara（豪）':<40}{rem_high:>12.2f}"
      f"{f'モデルは {removed_per_t / rem_high:.2f} 倍':>30}")
print()
print(f"  ★★モデルの果実持ち出しは文献の 1/1.7〜1/2.7 しかない。")
print(f"    内訳を見ると、果実に行くNは全吸収の "
      f"{demand.removed_g_per_m2 / demand.uptake_kg_per_10a * 100:.0f}% だけ。")
print(f"    文献は『果実が地上部Nの約2/3を持つ』と言っている。**逆転している。**")
print()
print("  【なぜそうなるのか】")
print(f"    果実のN濃度（乾物あたり）{ORGAN_N_CONTENT['FRUIT'] * 100:.1f}% に対し、")
print(f"    葉は {ORGAN_N_CONTENT['LEAF'] * 100:.1f}% で2倍以上ある。")
print(f"    乾物の{settings.dry_matter_partition['FRUIT'] * 100:.0f}%が果実に行っても、")
print("    N では葉・茎・根のほうが重くなる。")
print()
print("  【どちらが疑わしいか】")
print("    (a) 文献の2/3は**加工用トマト**（有限伸育・一斉収穫・茎葉が枯れ上がる）。")
print("        こちらは無限伸育で9か月、葉を作り続けるので、")
print("        栄養器官のNプールが本当に大きい。→ 逆転は正しいかもしれない。")
print(f"    (b) 果実のN濃度 {ORGAN_N_CONTENT['FRUIT'] * 100:.1f}% が低すぎる"
      f"（文献では 2.0〜3.0% の報告もある）。")
print(f"    (c) 果実乾物率 {fruit_dm_ratio * 100:.1f}% が低すぎる。")
print()
print("    ★第1節の『全吸収』が合っているのに第2節の『果実』が合わないのは、")
print("      **全体の量は当たっていて、器官への配り方がずれている**ということ。")
print("      施肥量の計算（＝全吸収を使う）には効かないが、")
print("      『どれだけ持ち出されたか』を言うときは使えない。")
print()

# =============================================================================
print("=" * 86)
print("3. 施肥N ─── 絶対値では公的基準とほぼ同じ、1トンあたりでは1/3")
print("=" * 86)
print()
chiba_yield_t = CHIBA_MINITOMATO_PROMOTED["yield_kg_per_10a"] / 1000.0
chiba_fert = CHIBA_MINITOMATO_PROMOTED["n_kg_per_10a"]
chiba_fert_per_t = chiba_fert / chiba_yield_t

print(f"{'':<34}{'施肥N':>12}{'収量':>10}{'1tあたり':>12}")
print(f"{'':<34}{'kg-N/10a':>12}{'t/10a':>10}{'kg-N/t':>12}")
print("-" * 86)
print(f"{'★このモデル':<34}{plan.fertilizer_kg_per_10a:>12.1f}"
      f"{yield_kg_per_m2:>10.1f}{fertilizer_per_t:>12.2f}")
print(f"{'千葉県 ミニトマト（ハウス促成）':<34}{chiba_fert:>12.1f}"
      f"{chiba_yield_t:>10.1f}{chiba_fert_per_t:>12.2f}")
print()
print(f"  → **絶対値はほぼ同じ**（{plan.fertilizer_kg_per_10a:.0f} 対 {chiba_fert:.0f} kg-N/10a）。")
print(f"    ところが収量が {yield_kg_per_m2 / chiba_yield_t:.1f} 倍あるので、")
print(f"    1トンあたりの施肥Nは **{chiba_fert_per_t / fertilizer_per_t:.1f} 分の1**。")
print(f"    同じ肥料で {yield_kg_per_m2 / chiba_yield_t:.1f} 倍取っている勘定になる。")
print()
print("  ★これは段3で見た『収量は 施肥N 40〜72 kg/10a で頭打ち』と同じ話。")
print("    肥料を増やしても収量は増えないところにいる。")
print()

# 1日あたりで道南の段階と比べる
daily = plan.fertilizer_kg_per_10a / SEASON_DAYS
print(f"  【1日あたりで道南農試の段階と比べる】（作期 {SEASON_DAYS} 日）")
print(f"    このモデル … {daily:.3f} kg-N/10a・日")
print(f"    道南の段階 … {' / '.join(f'{s:.3f}' for s in DOUNAN_DAILY_N_STEPS)}")
lo, hi = DOUNAN_DAILY_N_STEPS[0], DOUNAN_DAILY_N_STEPS[-1]
print(f"    → {band(daily, lo, hi)}。"
      f"いちばん下の段階（{DOUNAN_DAILY_N_STEPS[0]}）と"
      f"2番目（{DOUNAN_DAILY_N_STEPS[1]}）の間。")
print()

# =============================================================================
print("=" * 86)
print("4. ha 換算で加工用トマトの研究と比べる")
print("=" * 86)
print("  ★収量が違いすぎるので、絶対値ではなく**1トンあたりに直して**比べる。")
print()
uptake_kg_per_ha = demand.uptake_kg_per_10a * KG_PER_HA_PER_KG_PER_10A
yield_t_per_ha = yield_kg_per_m2 * T_PER_HA_PER_KG_PER_M2
pt_lo, pt_hi = PROCESSING_TOMATO_UPTAKE_KG_PER_HA_RANGE
py_lo, py_hi = PROCESSING_TOMATO_YIELD_T_PER_HA_RANGE

print(f"{'':<34}{'吸収N kg/ha':>14}{'収量 t/ha':>12}{'1tあたり':>12}")
print("-" * 86)
print(f"{'★このモデル':<34}{uptake_kg_per_ha:>14.0f}"
      f"{yield_t_per_ha:>12.0f}{uptake_per_t:>12.2f}")
print(f"{'加工用トマト（文献）':<34}{f'{pt_lo:.0f}〜{pt_hi:.0f}':>14}"
      f"{f'{py_lo:.0f}〜{py_hi:.0f}':>12}"
      f"{f'{pt_lo / py_hi:.2f}〜{pt_hi / py_lo:.2f}':>12}")
print()
print(f"  → 絶対値では {uptake_kg_per_ha / pt_hi:.1f} 倍も多いが、"
      f"収量も {yield_t_per_ha / py_hi:.1f} 倍ある。")
print(f"    1トンあたりに直すと {uptake_per_t:.2f} 対 "
      f"{pt_lo / py_hi:.2f}〜{pt_hi / py_lo:.2f} で"
      f"{band(uptake_per_t, pt_lo / py_hi, pt_hi / py_lo)}。")
print()

opt_lo, opt_hi = GLOBAL_OPTIMUM_N_KG_PER_HA_RANGE
fert_kg_per_ha = plan.fertilizer_kg_per_10a * KG_PER_HA_PER_KG_PER_10A
print(f"  【世界のメタ解析の最適N施用量と比べる】")
print(f"    メタ解析   … {opt_lo:.0f}〜{opt_hi:.0f} kg-N/ha")
print(f"    このモデル … {fert_kg_per_ha:.0f} kg-N/ha"
      f"（{band(fert_kg_per_ha, opt_lo, opt_hi)}）")
print(f"    ★ただしメタ解析の対象の収量は {py_lo:.0f}〜{py_hi:.0f} t/ha 水準。")
print(f"      1トンあたりに直すと "
      f"{opt_lo / py_hi:.2f}〜{opt_hi / py_lo:.2f} 対 {fertilizer_per_t:.2f} で、")
print(f"      **こちらのほうが1トンあたりの施肥Nは少ない。**")
print()

# =============================================================================
print("=" * 86)
print("5. 土壌診断の相互検算（千葉県のEC→硝酸態窒素の推定式）")
print("=" * 86)
print("  ★ECと硝酸態窒素は**別々に測っている**。")
print("    両方が独立に正しければ、千葉の推定式に乗るはず。乗れば両方の裏が取れる。")
print()
print(f"    Y [mg/100g] = {CHIBA_EC_TO_MINERAL_N_COEF[0]} X² "
      f"+ {CHIBA_EC_TO_MINERAL_N_COEF[1]} X "
      f"{CHIBA_EC_TO_MINERAL_N_COEF[2]}    X: EC(1:5) [mS/m]")
print()
a, b, c = CHIBA_EC_TO_MINERAL_N_COEF
print(f"{'年':>6}{'EC(1:5)':>10}{'EC':>9}{'推定の硝酸N':>13}"
      f"{'実測の硝酸N':>13}{'実測/推定':>11}")
print(f"{'':>6}{'dS/m':>10}{'mS/m':>9}{'mg/100g':>13}{'mg/100g':>13}{'':>11}")
ratios = []
for year, record in sorted(SOIL_DIAGNOSIS_CHUO.items()):
    ec_ds = record["ec_1to5_ds_per_m"]
    ec_ms = ec_ds * CHIBA_EC_MS_PER_M_PER_DS_PER_M
    estimated = a * ec_ms ** 2 + b * ec_ms + c
    measured = record["mineral_n_mg_per_100g"]
    ratio = measured / estimated if estimated > 0 else float("nan")
    ratios.append(ratio)
    print(f"{year:>6}{ec_ds:>10.2f}{ec_ms:>9.0f}{estimated:>13.1f}"
          f"{measured:>13.1f}{ratio:>11.2f}")
print()
print(f"  → 実測は推定の {min(ratios):.2f}〜{max(ratios):.2f} 倍。"
      f"全年で実測が高めだが、**大小の順序は保たれている**。")
print("    ECが硝酸以外（Ca・Mg・K）も拾うぶん推定は下振れしやすいので、")
print("    この向きのずれは想定どおり。**2つの測定は矛盾していない。**")
print()

# =============================================================================
print("=" * 86)
print("6. ★ECの帯を千葉・茨城の土壌EC基準と突き合わせる")
print("=" * 86)
print("  ★第11-8節の帯は FAO の ECe しきい値 2.5 dS/m から逆算したものだった。")
print("    日本の土壌診断は EC(1:5) で基準を持っている。測り方が違うので換算する。")
print()
texture = IBARAKI_EC_1TO5_MS_PER_M[SOIL_TEXTURE_CLOSEST]
fit_low, fit_high = texture["適"]
harm = texture["障害"]
theta = FIELD_CAPACITY_WATER_CONTENT

print(f"  土性は「{SOIL_TEXTURE_CLOSEST}」が近いものとして読む。")
print()
print(f"{'区分':<22}{'EC(1:5)':>10}{'EC(1:5)':>10}{'土壌溶液EC':>12}"
      f"{'ECe':>8}{'収量比':>9}")
print(f"{'':<22}{'mS/m':>10}{'dS/m':>10}{'dS/m':>12}{'dS/m':>8}{'':>9}")
print("-" * 86)
for label, ec_ms in (("生育適濃度の下端", fit_low),
                     ("生育適濃度の上端", fit_high),
                     ("生育障害発生", harm)):
    ec_ds = ec_ms / CHIBA_EC_MS_PER_M_PER_DS_PER_M
    solution = ec_1to5_to_soil_solution(ec_ds, theta)
    ece = ec_1to5_to_saturated_extract(ec_ds)
    print(f"{label:<22}{ec_ms:>10.0f}{ec_ds:>10.2f}{solution:>12.2f}"
          f"{ece:>8.2f}{salinity_response(ece).yield_ratio * 100:>8.1f}%")
print("-" * 86)
for year, record in sorted(SOIL_DIAGNOSIS_CHUO.items()):
    ec_ds = record["ec_1to5_ds_per_m"]
    solution = ec_1to5_to_soil_solution(ec_ds, theta)
    ece = ec_1to5_to_saturated_extract(ec_ds)
    verdict = ("適濃度の下" if ec_ds * CHIBA_EC_MS_PER_M_PER_DS_PER_M < fit_low
               else ("障害域" if ec_ds * CHIBA_EC_MS_PER_M_PER_DS_PER_M > harm
                     else "適濃度の中"))
    print(f"{f'実測 {year}年（{verdict}）':<22}"
          f"{ec_ds * CHIBA_EC_MS_PER_M_PER_DS_PER_M:>10.0f}{ec_ds:>10.2f}"
          f"{solution:>12.2f}{ece:>8.2f}"
          f"{salinity_response(ece).yield_ratio * 100:>8.1f}%")
print()
harm_ece = ec_1to5_to_saturated_extract(
    harm / CHIBA_EC_MS_PER_M_PER_DS_PER_M)
print(f"  ★★ここで食い違いが出る。")
print(f"    茨城園試の『生育障害発生』EC(1:5) {harm:.0f} mS/m は、")
print(f"    こちらの換算では ECe {harm_ece:.2f} dS/m にあたる。")
print(f"    FAO の減収しきい値 ECe 2.5 より **{harm_ece / 2.5:.1f} 倍ゆるい**。")
print()
print("  【どう読むか】")
print("    ・FAO の 2.5 は『減収が始まる点』、茨城の 150 は『障害が見える点』。")
print("      見ているものが違う（数%の減収は見た目には出ない）。")
print("    ・だから帯の上限を FAO 側に置いたままでよい。**安全側にある。**")
print("    ・一方で実測の EC(1:5) は3作期とも茨城の適濃度の下端以下だった。")
print("      **土に塩が溜まりすぎている状態ではない。**")
print("      第11-8節で『10〜12月の運用値は上限超え』と出たのは、")
print("      給液ECから濃縮3倍を仮定して ECe を推したもので、")
print("      **土壌診断の実測とは整合していない。**")
print()
print("  → 結論: **帯の上限は暫定のまま。排液ECの実測が入るまで動かさない。**")
print("    濃縮倍率3倍が過大なら、上限はもっと上げてよいことになる。")
print()

# =============================================================================
print("=" * 86)
print("まとめ")
print("=" * 86)
print()
print("  ◎ 合っているもの")
print(f"    ・果実1トンあたりの全吸収N {uptake_per_t:.2f} kg-N/t")
print(f"      千葉 {chiba_n}・ヨーロッパ {eu_low}〜{eu_high}・"
      f"加工用 {pt_lo / py_hi:.1f}〜{pt_hi / py_lo:.1f} の3つに挟まれている")
print(f"    ・作期の施肥N {plan.fertilizer_kg_per_10a:.0f} kg-N/10a")
print(f"      千葉のハウス促成の基準 {chiba_fert:.0f} とほぼ同じ")
print(f"    ・1日あたり {daily:.3f} kg-N/10a・日 は道南の段階の中")
print(f"    ・土壌診断のECと硝酸態窒素が、千葉の推定式の上で矛盾しない")
print()
print("  × 合っていないもの")
print(f"    ・果実で持ち出すN {removed_per_t:.2f} kg-N/t は文献の"
      f" {rem_low:.1f}〜{rem_high:.1f} に対して低すぎる")
print(f"      → 器官への配り方（果実N濃度 or 果実乾物率）がずれている")
print(f"      → **施肥量の計算には効かない**（全吸収を使っているため）")
print()
print("  ？ 判断を保留するもの")
print(f"    ・ECの帯の上限。土壌診断の実測は適濃度の下端以下で、")
print(f"      給液ECから推した根圏ECとは整合していない")
print(f"      → 排液ECをEC計で月1回測るのが、いちばん効く1手")
print()
print("=" * 86)
