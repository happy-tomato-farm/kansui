"""塩類濃度（EC）が栽培に効く経路を数字にする。

【このモジュールがある理由】（2026-10-01 ユーザー聞き取り）
施肥の時期配分は、N栄養の都合で決めたものではない。デルフィのコーチから

    「冬は窒素を多めに流して EC を上げる、春先は EC を下げる」

という指導を受けてこの運用になっている。狙っているのは **EC** で、
Nはその手段として動いている。だから `core/nitrogen.py`（N需要）だけでは
妥当性を判断できない。

【ECが効く3つの経路】

    (1) 浸透ポテンシャル → 吸水を抑える   ← このモジュールの中心
    (2) 果実の糖度       → 上がる
    (3) 収量             → 下がる

(2) と (3) は同じ現象の表と裏。果実への水の流入が減るので、糖が薄まらず
濃くなる（糖度↑）が、同時に果実は大きくならない（収量↓）。
どちらを取るかという選択になる。

【★なぜECがマトリックポテンシャルより効くのか】
`core/soil.py` の soil_to_leaf_driving_force にすでに書いてある。

    圃場容水量（pF 1.8）のマトリックポテンシャル   6.2 J/kg
    土壌溶液 EC 2.0 dS/m の浸透ポテンシャル       72   J/kg   ← 12倍

**十分に湿った土では、根が感じる水ポテンシャルはほぼECで決まる。**
第9章の式9.20・9.22 は含水率だけで吸水を決めているので、ECを入れないと
春の吸水不足を説明しきれない（README の「土は湿っているのに吸水が届かない」）。

【★ECの3つの基準を混同しないこと】
ECは「どれだけ水で薄めて測ったか」で値が変わる。

    EC(1:5)     土1に水5。日本の土壌診断の標準。手元の実測はこれ
    土壌溶液EC   現地の含水率のまま。**根が感じるのはこれ**
    ECe         飽和ペースト。**文献の塩害しきい値はこの基準**

このモジュールは3つを行き来する関数を持つ。引数名と戻り値の名前に
どの基準かを必ず書いてある。

【実測していないもの（すべて config 第11節に幅つきで置いてある）】
・土の仮比重 … EC基準の換算倍率にそのまま効く
・肥料のEC係数 … 給液ECの絶対値を決める。EC計で1点測れば確定する
・原水のEC … 水質分析をしていない
・作物が吸う塩の量 … 根圏の塩類収支で効く
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from config import (
    EC_EXTRACT_WATER_RATIO,
    FEED_EC_BAND_BY_MONTH,
    FERTILIZER_EC_DS_PER_M_PER_G_PER_L,
    OSMOTIC_POTENTIAL_J_PER_KG_PER_DS_PER_M,
    RAW_WATER_EC_DS_PER_M,
    SALINITY_BRIX_REFERENCE_ECE,
    SALINITY_BRIX_SLOPE_PER_DS_PER_M,
    SALINITY_YIELD_SLOPE_PER_DS_PER_M,
    SALINITY_YIELD_THRESHOLD_ECE,
    SALT_UPTAKE_PER_N_DS_PER_M_MM,
    SOIL_BULK_DENSITY_G_PER_CM3,
    SOIL_RETENTION,
)

#: m → mm
MM_PER_M = 1000.0

#: g/cm³ → kg/m³
KG_PER_M3_PER_G_PER_CM3 = 1000.0

#: mg/100g → mg/kg
MG_PER_KG_PER_MG_PER_100G = 10.0


# =============================================================================
# 1. EC の基準を行き来する
# =============================================================================

def gravimetric_water_content(
    volumetric_water_content: float,
    bulk_density_g_per_cm3: float = SOIL_BULK_DENSITY_G_PER_CM3,
) -> float:
    """体積含水率 [m³/m³] を重量含水率 [g-水/g-乾土] に直す。

        重量含水率 = 体積含水率 ÷ 仮比重

    EC の基準を換算するときに、分母として要る。
    """
    if bulk_density_g_per_cm3 <= 0.0:
        raise ValueError(
            f"土の仮比重は正の値でなければならない。"
            f"渡された値: {bulk_density_g_per_cm3} g/cm³"
        )
    if not 0.0 < volumetric_water_content < 1.0:
        raise ValueError(
            f"体積含水率は 0 より大きく 1 未満。"
            f"渡された値: {volumetric_water_content}"
        )
    return volumetric_water_content / bulk_density_g_per_cm3


def ec_1to5_to_soil_solution(
    ec_1to5_ds_per_m: float,
    volumetric_water_content: float,
    bulk_density_g_per_cm3: float = SOIL_BULK_DENSITY_G_PER_CM3,
) -> float:
    """土壌診断の EC(1:5) を、現地の土壌溶液の EC [dS/m] に直す。

        土壌溶液EC = EC(1:5) × 5 ÷ 現地の重量含水率

    1:5 は土1gに水5gを加えて測るので、現地（重量含水率 0.47 g/g くらい）
    より10倍ほど薄い。そのぶんを戻す。

    ★この線形換算は「塩が全部溶けている」という仮定の上に立つ。
      石膏など溶けきらない塩があると過大に出る。この土はCECが 42 me と
      高く腐植も10%以上あるので、1:5抽出で交換性の塩基が余分に出て
      **土壌溶液ECを過大に見積もる向き**と見ておくこと。
    """
    if ec_1to5_ds_per_m < 0.0:
        raise ValueError(
            f"EC が負: {ec_1to5_ds_per_m} dS/m。測定値を疑うこと。"
        )
    field = gravimetric_water_content(
        volumetric_water_content, bulk_density_g_per_cm3)
    return ec_1to5_ds_per_m * EC_EXTRACT_WATER_RATIO / field


def ec_1to5_to_saturated_extract(
    ec_1to5_ds_per_m: float,
    bulk_density_g_per_cm3: float = SOIL_BULK_DENSITY_G_PER_CM3,
) -> float:
    """土壌診断の EC(1:5) を、飽和抽出の ECe [dS/m] に直す。

    文献の塩害しきい値（FAO のトマト 2.5 dS/m など）はこの基準なので、
    比べるときはここを通す。飽和の含水率は全孔隙率（θs）を使う。
    """
    return ec_1to5_to_soil_solution(
        ec_1to5_ds_per_m, SOIL_RETENTION["THETA_S"], bulk_density_g_per_cm3)


def soil_solution_to_saturated_extract(
    ec_solution_ds_per_m: float,
    volumetric_water_content: float,
) -> float:
    """現地の土壌溶液 EC を、飽和抽出の ECe [dS/m] に直す。

        ECe = 土壌溶液EC × 現地の含水率 ÷ 飽和含水率

    薄めるぶんだけ下がる。仮比重は分子と分母で打ち消すので要らない。
    """
    if ec_solution_ds_per_m < 0.0:
        raise ValueError(f"EC が負: {ec_solution_ds_per_m} dS/m")
    if not 0.0 < volumetric_water_content < 1.0:
        raise ValueError(
            f"体積含水率は 0 より大きく 1 未満。"
            f"渡された値: {volumetric_water_content}"
        )
    return (ec_solution_ds_per_m * volumetric_water_content
            / SOIL_RETENTION["THETA_S"])


def mineral_n_to_g_per_m2(
    mineral_n_mg_per_100g: float,
    depth_m: float,
    bulk_density_g_per_cm3: float = SOIL_BULK_DENSITY_G_PER_CM3,
) -> float:
    """土壌診断の無機態窒素 [mg/100g] を [g-N/m²] に直す。

        土の重さ [kg/m²] = 深さ [m] × 仮比重 [kg/m³]
        N [g/m²] = (mg/100g × 10) [mg/kg] × 土の重さ ÷ 1000

    ★1 g-N/m² = 1 kg-N/10a なので、戻り値はそのまま kg/10a としても読める。

    ★診断の採土深（ふつう 0〜15 cm）より深い層へ引き伸ばすと過大になる。
      表層ほど濃いからだ。深さを変えて幅を見ること。
    """
    if mineral_n_mg_per_100g < 0.0:
        raise ValueError(
            f"無機態窒素が負: {mineral_n_mg_per_100g} mg/100g"
        )
    if depth_m <= 0.0:
        raise ValueError(f"深さは正の値。渡された値: {depth_m} m")
    soil_kg_per_m2 = depth_m * bulk_density_g_per_cm3 * KG_PER_M3_PER_G_PER_CM3
    mg_per_kg = mineral_n_mg_per_100g * MG_PER_KG_PER_MG_PER_100G
    return mg_per_kg * soil_kg_per_m2 / 1000.0


# =============================================================================
# 2. 浸透ポテンシャルと、吸水への効き
# =============================================================================

def osmotic_potential_j_per_kg(ec_solution_ds_per_m: float) -> float:
    """土壌溶液の EC から浸透ポテンシャルの大きさ [J/kg] を出す。

        |ψo| = 36 × EC [dS/m]

    （浸透圧 π[bar] ≒ 0.36 × EC、1 bar ≒ 100 J/kg より）

    ★「大きさ」を正の値で返す。符号は負なので、土の水ポテンシャルを
      合計するときは引き算になる。core/soil.py の
      potential_from_water_content も同じ約束（正の大きさ）なので、
      そのまま足し合わせればよい。

    【桁の感覚】
        EC 0.5 →  18 J/kg（pF 2.26 相当）
        EC 1.0 →  36 J/kg（pF 2.56）
        EC 2.0 →  72 J/kg（pF 2.86）
        EC 4.0 → 144 J/kg（pF 3.17）
    圃場容水量のマトリックポテンシャルは 6.2 J/kg しかない。
    """
    if ec_solution_ds_per_m < 0.0:
        raise ValueError(f"EC が負: {ec_solution_ds_per_m} dS/m")
    return OSMOTIC_POTENTIAL_J_PER_KG_PER_DS_PER_M * ec_solution_ds_per_m


def total_soil_potential_j_per_kg(
    matric_potential_j_per_kg: float,
    ec_solution_ds_per_m: float,
) -> float:
    """根が感じる土の水ポテンシャルの大きさ [J/kg]。

        |ψ土| = |ψマトリック| + |ψ浸透|

    どちらも「正の大きさ」で渡し、正の大きさで返す。
    """
    if matric_potential_j_per_kg < 0.0:
        raise ValueError(
            f"マトリックポテンシャルは正の大きさで渡すこと。"
            f"渡された値: {matric_potential_j_per_kg} J/kg"
        )
    return (matric_potential_j_per_kg
            + osmotic_potential_j_per_kg(ec_solution_ds_per_m))


def uptake_driving_force_ratio(
    matric_potential_j_per_kg: float,
    ec_solution_ds_per_m: float,
    leaf_potential_j_per_kg: float = 1000.0,
) -> float:
    """ECを入れたときの吸水の駆動力が、ECゼロのときの何倍になるか。

        駆動力 = |ψ葉| − |ψ土|
        比     = (|ψ葉| − |ψマトリック| − |ψ浸透|) ÷ (|ψ葉| − |ψマトリック|)

    【これが答える問い】
    「ECを上げると、どれだけ水が吸いにくくなるのか」

    【桁の感覚】葉内水ポテンシャル −1.0 MPa・圃場容水量のとき
        EC 0.5 → 0.983（1.7%減）
        EC 1.0 → 0.964（3.6%減）
        EC 2.0 → 0.928（7.2%減）
        EC 4.0 → 0.855（14.5%減）

    ★一見小さいが、**蒸散要求が大きい時期にはこれが律速になる**。
      5月の蒸散要求は 7.6 mm/日で、吸水の上限（根の能力 6 mm/日、
      config 第7-6節）にすでに当たっている。そこから7%削られると、
      届かない量がそのまま果実の肥大を削る。
      冬は蒸散要求が 1.5 mm/日しかないので、7%減っても上限に当たらない。
      **同じEC・同じ減り方でも、効く時期と効かない時期がある。**

    引数の leaf_potential_j_per_kg は「大きさ」で渡す（1000 なら −1.0 MPa）。
    """
    if leaf_potential_j_per_kg <= 0.0:
        raise ValueError(
            f"葉内水ポテンシャルは正の大きさで渡すこと。"
            f"渡された値: {leaf_potential_j_per_kg} J/kg"
        )
    without_salt = leaf_potential_j_per_kg - matric_potential_j_per_kg
    if without_salt <= 0.0:
        raise ValueError(
            f"マトリックポテンシャル {matric_potential_j_per_kg} J/kg が"
            f"葉内水ポテンシャル {leaf_potential_j_per_kg} J/kg に"
            f"届いている。この条件では吸水できず、比を出せない。"
        )
    with_salt = without_salt - osmotic_potential_j_per_kg(ec_solution_ds_per_m)
    return max(with_salt, 0.0) / without_salt


# =============================================================================
# 3. 収量と糖度への効き（同じ現象の表と裏）
# =============================================================================

@dataclass(frozen=True)
class SalinityResponse:
    """ある ECe のときの、収量と糖度への効き。"""

    ece_ds_per_m: float
    yield_ratio: float          # 塩害のない場合に対する収量比（0〜1）
    brix_change: float          # 基準に対する糖度の変化 [°Brix]

    def describe(self) -> str:
        return (f"ECe {self.ece_ds_per_m:.2f} dS/m → "
                f"収量 {self.yield_ratio * 100:.1f}% ／ "
                f"糖度 {self.brix_change:+.2f} °Brix")


def salinity_response(ece_ds_per_m: float) -> SalinityResponse:
    """ECe から収量比と糖度の変化を出す。

    【収量】FAO方式。しきい値を超えた分に比例して下がる。
        収量比 = 1 − 傾き × (ECe − しきい値)
        トマト: しきい値 2.5 dS/m、傾き 9.9 %/dS/m

    【糖度】1 dS/m あたり 0.4 °Brix 上がる（施設トマトの報告の中央）。
        こちらはしきい値を置かない（低ECでも薄まる方向に効く）。

    ★どちらも露地・作期平均を前提にした係数で、施設の長段どりに
      そのまま当たるとは限らない。**向きと桁を見るためのもの**。
    """
    if ece_ds_per_m < 0.0:
        raise ValueError(f"ECe が負: {ece_ds_per_m} dS/m")

    excess = max(ece_ds_per_m - SALINITY_YIELD_THRESHOLD_ECE, 0.0)
    yield_ratio = max(1.0 - SALINITY_YIELD_SLOPE_PER_DS_PER_M * excess, 0.0)
    brix = SALINITY_BRIX_SLOPE_PER_DS_PER_M * (
        ece_ds_per_m - SALINITY_BRIX_REFERENCE_ECE)
    return SalinityResponse(
        ece_ds_per_m=ece_ds_per_m,
        yield_ratio=yield_ratio,
        brix_change=brix,
    )


# =============================================================================
# 4. 肥料から給液ECへ
# =============================================================================

def feed_ec_from_n_concentration(
    n_mg_per_l: float,
    fertilizer_n_fraction: float,
    ec_coefficient: float = FERTILIZER_EC_DS_PER_M_PER_G_PER_L,
    raw_water_ec_ds_per_m: float = RAW_WATER_EC_DS_PER_M,
) -> float:
    """液肥のN濃度から給液の EC [dS/m] を出す。

        肥料濃度 [g/L] = N濃度 [mg-N/L] ÷ 1000 ÷ 肥料のN比
        給液EC = 肥料濃度 × EC係数 ＋ 原水のEC

    液肥混入機レシピが決めるのは「1日のN量」なので、ECはここを通して出す。

    Args:
        n_mg_per_l: 給液のN濃度 [mg-N/L]
        fertilizer_n_fraction: 肥料のN比（トケル養液配合1号なら 0.10）
        ec_coefficient: 肥料 1 g/L あたりのEC [dS/m]（★実測していない）
        raw_water_ec_ds_per_m: 原水のEC（★水質分析をしていない）

    【実績の桁】（直近4作期・中央。README「段4〜5で答える問い」の表）
        12月 111 mg-N/L → 肥料 1.11 g/L → 給液EC 1.59 dS/m
         5月  18 mg-N/L → 肥料 0.18 g/L → 給液EC 0.38 dS/m
    一般的な養液（2.0〜3.0 dS/m）より薄い。土耕で土が養分を出すので、
    これで足りているということだ。
    """
    if n_mg_per_l < 0.0:
        raise ValueError(f"N濃度が負: {n_mg_per_l} mg/L")
    if not 0.0 < fertilizer_n_fraction <= 1.0:
        raise ValueError(
            f"肥料のN比は 0 より大きく 1 以下。"
            f"渡された値: {fertilizer_n_fraction}。"
            f"10% なら 0.10 と書く（10 ではない）。"
        )
    if ec_coefficient <= 0.0:
        raise ValueError(f"EC係数は正の値。渡された値: {ec_coefficient}")
    if raw_water_ec_ds_per_m < 0.0:
        raise ValueError(
            f"原水のECが負: {raw_water_ec_ds_per_m} dS/m"
        )
    fertilizer_g_per_l = n_mg_per_l / 1000.0 / fertilizer_n_fraction
    return fertilizer_g_per_l * ec_coefficient + raw_water_ec_ds_per_m


def n_concentration_from_feed_ec(
    feed_ec_ds_per_m: float,
    fertilizer_n_fraction: float,
    ec_coefficient: float = FERTILIZER_EC_DS_PER_M_PER_G_PER_L,
    raw_water_ec_ds_per_m: float = RAW_WATER_EC_DS_PER_M,
) -> float:
    """★逆向き: 目標の給液ECから、液肥のN濃度 [mg-N/L] を出す。

        肥料濃度 [g/L] = ( 給液EC − 原水のEC ) ÷ EC係数
        N濃度          = 肥料濃度 × N比 × 1000

    feed_ec_from_n_concentration の逆算。朝に「ECをいくつにするか」から
    決めるときに使う。デルフィの指導はこの向きの考え方だ。

    原水のECより低い目標を渡したら、肥料を入れない（0 を返す）。
    水だけでそのECなので、それ以上は下げようがない。
    """
    if feed_ec_ds_per_m < 0.0:
        raise ValueError(f"給液ECが負: {feed_ec_ds_per_m} dS/m")
    if not 0.0 < fertilizer_n_fraction <= 1.0:
        raise ValueError(
            f"肥料のN比は 0 より大きく 1 以下。"
            f"渡された値: {fertilizer_n_fraction}。"
            f"10% なら 0.10 と書く（10 ではない）。"
        )
    if ec_coefficient <= 0.0:
        raise ValueError(f"EC係数は正の値。渡された値: {ec_coefficient}")
    fertilizer_g_per_l = max(
        feed_ec_ds_per_m - raw_water_ec_ds_per_m, 0.0) / ec_coefficient
    return fertilizer_g_per_l * fertilizer_n_fraction * 1000.0


def daily_n_from_feed_ec(
    feed_ec_ds_per_m: float,
    irrigation_l_per_m2: float,
    fertilizer_n_fraction: float,
    ec_coefficient: float = FERTILIZER_EC_DS_PER_M_PER_G_PER_L,
    raw_water_ec_ds_per_m: float = RAW_WATER_EC_DS_PER_M,
) -> float:
    """目標の給液ECと潅水量から、1日のN量 [kg-N/10a] を出す。

        N量 [g/m²] = N濃度 [mg/L] × 潅水量 [L/m²] ÷ 1000

    ★1 g-N/m² = 1 kg-N/10a なので、戻り値はそのまま液肥混入機レシピの
      「1日のN量」に入れられる。

    【この向きで決める意味】
    ECを固定すると、潅水量が増えればN量も自動で増える。液肥混入機レシピが
    「倍率を固定すればECが一定になり、潅水量を変えたぶんだけN量が増減する」
    と書いているのと同じことだ。レシピ側の作りと噛み合う。
    """
    if irrigation_l_per_m2 < 0.0:
        raise ValueError(f"潅水量が負: {irrigation_l_per_m2} L/m²")
    n_mg_per_l = n_concentration_from_feed_ec(
        feed_ec_ds_per_m, fertilizer_n_fraction,
        ec_coefficient, raw_water_ec_ds_per_m)
    return n_mg_per_l * irrigation_l_per_m2 / 1000.0


# =============================================================================
# 4-2. 給液ECの帯で挟む（量はNで決める）
# =============================================================================
# 【考え方】
# 1日のN量はN需要から決める。そのN量を潅水量で割ると給液ECが決まるので、
# それが月別の帯（config 第11-8節）から外れていたら、帯の端まで戻す。
#
#   ・上限で切られる … 濃すぎる。吸水が削られ、強光期なら尻腐果の側
#   ・下限で持ち上げ … 薄すぎる。暗い日に液肥がほとんど出なくなるのを防ぐ
#   ・帯の中          … そのまま。ECは何もしない
#
# ★「ECを目標にしてN量を決める」のではない。向きが逆。
#   ECは危ない側に出たときだけ効く括りで、ふだんは黙っている。
#   根拠は config 第11-8節（道南農試・愛知県の事例・EC×蒸散の文献）。

@dataclass(frozen=True)
class BandedDailyN:
    """帯で挟んだあとの1日のN量。

    Attributes:
        daily_n_kg_per_10a: 挟んだあとのN量 [kg-N/10a]。レシピに入れる値
        requested_n_kg_per_10a: 挟む前の、N需要から出した値 [kg-N/10a]
        feed_ec_ds_per_m: 挟んだあとの給液EC [dS/m]
        requested_feed_ec_ds_per_m: 挟む前の給液EC [dS/m]
        band_low_ds_per_m: その月の帯の下限 [dS/m]
        band_high_ds_per_m: その月の帯の上限 [dS/m]
        clamped_by: "high" / "low" / "" （挟まれなかったとき）
    """

    daily_n_kg_per_10a: float
    requested_n_kg_per_10a: float
    feed_ec_ds_per_m: float
    requested_feed_ec_ds_per_m: float
    band_low_ds_per_m: float
    band_high_ds_per_m: float
    clamped_by: str

    @property
    def is_clamped(self) -> bool:
        """帯に当たったかどうか。"""
        return self.clamped_by != ""

    def describe(self) -> str:
        """画面や検証スクリプトに出す一行。"""
        if not self.is_clamped:
            return (
                f"給液EC {self.feed_ec_ds_per_m:.2f} dS/m は帯"
                f"（{self.band_low_ds_per_m:.2f}〜"
                f"{self.band_high_ds_per_m:.2f}）の中。ECは効いていない"
            )
        side = "上限" if self.clamped_by == "high" else "下限"
        return (
            f"N需要から出した給液EC {self.requested_feed_ec_ds_per_m:.2f} dS/m が"
            f"{side}に当たったので {self.feed_ec_ds_per_m:.2f} dS/m まで戻した"
            f"（N量 {self.requested_n_kg_per_10a:.3f} → "
            f"{self.daily_n_kg_per_10a:.3f} kg-N/10a）"
        )


def clamp_daily_n_to_ec_band(
    daily_n_kg_per_10a: float,
    irrigation_l_per_m2: float,
    month: int,
    fertilizer_n_fraction: float,
    band_by_month: dict = FEED_EC_BAND_BY_MONTH,
    ec_coefficient: float = FERTILIZER_EC_DS_PER_M_PER_G_PER_L,
    raw_water_ec_ds_per_m: float = RAW_WATER_EC_DS_PER_M,
) -> BandedDailyN:
    """N需要から出した1日のN量を、その月の給液ECの帯で挟む。

    Args:
        daily_n_kg_per_10a: N需要から出した値 [kg-N/10a]
        irrigation_l_per_m2: その日の潅水量 [L/m²]
        month: 月（1〜12）。帯は月ごとに違う（config 第11-8節）
        fertilizer_n_fraction: 肥料のN比
        band_by_month: 月 → (下限, 上限) [dS/m]
        ec_coefficient: 肥料 1 g/L あたりのEC [dS/m]
        raw_water_ec_ds_per_m: 原水のEC [dS/m]

    潅水量が 0 のときは割り算ができないので、挟まずにそのまま返す。
    N量 0 とECの対応が決まらないだけで、異常ではない（定植前など）。
    """
    if daily_n_kg_per_10a < 0.0:
        raise ValueError(
            f"1日のN量が負: {daily_n_kg_per_10a} kg-N/10a。"
            f"N需要の計算（core/nitrogen.py）を疑うこと。"
        )
    if irrigation_l_per_m2 < 0.0:
        raise ValueError(f"潅水量が負: {irrigation_l_per_m2} L/m²")
    if month not in band_by_month:
        raise ValueError(
            f"月 {month} の給液ECの帯が config に無い。"
            f"持っているのは {sorted(band_by_month)}。"
        )

    low, high = band_by_month[month]
    if low > high:
        raise ValueError(
            f"{month}月の帯の下限 {low} が上限 {high} を超えている"
            f"（config 第11-8節）。"
        )

    # 潅水量が無い日は、N量とECの対応が決まらない。そのまま通す。
    if irrigation_l_per_m2 <= 0.0:
        return BandedDailyN(
            daily_n_kg_per_10a=daily_n_kg_per_10a,
            requested_n_kg_per_10a=daily_n_kg_per_10a,
            feed_ec_ds_per_m=0.0,
            requested_feed_ec_ds_per_m=0.0,
            band_low_ds_per_m=low,
            band_high_ds_per_m=high,
            clamped_by="",
        )

    # 1 kg-N/10a = 1 g-N/m² なので、潅水量で割れば mg-N/L になる
    n_mg_per_l = daily_n_kg_per_10a / irrigation_l_per_m2 * 1000.0
    requested_ec = feed_ec_from_n_concentration(
        n_mg_per_l, fertilizer_n_fraction,
        ec_coefficient, raw_water_ec_ds_per_m)

    # 栽培していない月（帯が 0〜0）は液肥を出さない
    if high <= 0.0:
        return BandedDailyN(
            daily_n_kg_per_10a=0.0,
            requested_n_kg_per_10a=daily_n_kg_per_10a,
            feed_ec_ds_per_m=0.0,
            requested_feed_ec_ds_per_m=requested_ec,
            band_low_ds_per_m=low,
            band_high_ds_per_m=high,
            clamped_by="high" if daily_n_kg_per_10a > 0.0 else "",
        )

    if requested_ec > high:
        clamped_by = "high"
        target_ec = high
    elif requested_ec < low:
        clamped_by = "low"
        target_ec = low
    else:
        return BandedDailyN(
            daily_n_kg_per_10a=daily_n_kg_per_10a,
            requested_n_kg_per_10a=daily_n_kg_per_10a,
            feed_ec_ds_per_m=requested_ec,
            requested_feed_ec_ds_per_m=requested_ec,
            band_low_ds_per_m=low,
            band_high_ds_per_m=high,
            clamped_by="",
        )

    return BandedDailyN(
        daily_n_kg_per_10a=daily_n_from_feed_ec(
            target_ec, irrigation_l_per_m2, fertilizer_n_fraction,
            ec_coefficient, raw_water_ec_ds_per_m),
        requested_n_kg_per_10a=daily_n_kg_per_10a,
        feed_ec_ds_per_m=target_ec,
        requested_feed_ec_ds_per_m=requested_ec,
        band_low_ds_per_m=low,
        band_high_ds_per_m=high,
        clamped_by=clamped_by,
    )


# =============================================================================
# 5. 根圏の塩類収支（core/soil_nitrogen.py と同じ作り）
# =============================================================================
# 塩も「よく混ざった水槽」で追う。Nとの違いは2つ。
#   (1) 塩は無機化で増えない（有機物から出るのは主にNとS）
#   (2) 作物は養分イオンを吸うので、その分は減る
#
# 塩の量は「EC × 水量」で持つ。単位は dS/m·mm。
#   EC [dS/m] = 塩の量 [dS/m·mm] ÷ 根圏の水量 [mm]
# Nを g/m² で持って濃度を割り算で出すのと、まったく同じ仕掛けだ。

@dataclass(frozen=True)
class DailySaltInput:
    """1日ぶんの入力。"""

    date: str
    irrigation_mm: float
    drainage_mm: float
    feed_ec_ds_per_m: float        # その日の給液のEC
    n_uptake_g_per_m2: float       # その日の吸収N（塩の持ち出しに使う）
    water_content: float           # その日の根圏の体積含水率

    def __post_init__(self) -> None:
        for name, value in (
            ("irrigation_mm", self.irrigation_mm),
            ("drainage_mm", self.drainage_mm),
            ("feed_ec_ds_per_m", self.feed_ec_ds_per_m),
            ("n_uptake_g_per_m2", self.n_uptake_g_per_m2),
        ):
            if value < 0.0:
                raise ValueError(f"{self.date}: {name} が負（{value}）。")
        # ★流亡が潅水を超えることは正常（溜まった水が抜ける日）。
        #   core/soil_nitrogen.py の同じ注意書きを読むこと。
        if not 0.0 < self.water_content < 1.0:
            raise ValueError(
                f"{self.date}: 体積含水率が範囲外（{self.water_content}）。"
            )


@dataclass(frozen=True)
class DailySalt:
    """1日ぶんの結果。塩の量は dS/m·mm。"""

    date: str
    start_salt: float
    added_salt: float            # 潅水で入ってきた塩
    leached_salt: float          # 流亡した塩
    removed_salt: float          # 作物が吸った塩
    end_salt: float
    root_zone_water_mm: float
    start_ec_solution: float     # 日の初めの土壌溶液EC [dS/m]
    end_ec_solution: float
    end_ece: float               # 飽和抽出に直したEC（文献と比べる用）

    @property
    def balance_residual(self) -> float:
        """収支の閉じ具合。0 に近くなければ計算が壊れている。"""
        return (self.start_salt + self.added_salt
                - self.leached_salt - self.removed_salt - self.end_salt)


def step_one_day_salt(
    start_salt: float,
    day: DailySaltInput,
    root_zone_depth_m: float,
    salt_uptake_per_n: float = SALT_UPTAKE_PER_N_DS_PER_M_MM,
    substeps: int = 24,
) -> DailySalt:
    """1日ぶんの塩類収支を進める。

    Nと同じく1日を刻んで、「入れる → 流亡させる → 吸わせる」を繰り返す。

    Args:
        start_salt: 日の初めの塩の量 [dS/m·mm]
        day: その日の入力
        root_zone_depth_m: 根群域の深さ [m]
        salt_uptake_per_n: 吸収N 1 g あたり持ち出される塩 [dS/m·mm]
        substeps: 1日の刻み数
    """
    if start_salt < 0.0:
        raise ValueError(
            f"{day.date}: 日の初めの塩の量が負（{start_salt}）。"
            f"前日の計算が壊れている。"
        )
    if substeps < 1:
        raise ValueError(f"刻み数は1以上。渡された値: {substeps}")

    water_mm = day.water_content * root_zone_depth_m * MM_PER_M
    if day.drainage_mm > water_mm:
        raise ValueError(
            f"{day.date}: 流亡 {day.drainage_mm:.1f} mm が根圏の水量 "
            f"{water_mm:.1f} mm を超えている。水収支の結果か単位を疑うこと。"
        )

    added_total = day.irrigation_mm * day.feed_ec_ds_per_m
    removal_target = day.n_uptake_g_per_m2 * salt_uptake_per_n

    added_step = added_total / substeps
    removal_step = removal_target / substeps
    retention = math.exp(-(day.drainage_mm / substeps) / water_mm)

    salt = start_salt
    leached_total = 0.0
    removed_total = 0.0
    for _ in range(substeps):
        salt += added_step
        after = salt * retention
        leached_total += salt - after
        salt = after
        # 存在する以上は持ち出せない
        removed = min(removal_step, salt)
        salt -= removed
        removed_total += removed

    start_ec = start_salt / water_mm
    end_ec = salt / water_mm
    return DailySalt(
        date=day.date,
        start_salt=start_salt,
        added_salt=added_total,
        leached_salt=leached_total,
        removed_salt=removed_total,
        end_salt=salt,
        root_zone_water_mm=water_mm,
        start_ec_solution=start_ec,
        end_ec_solution=end_ec,
        end_ece=soil_solution_to_saturated_extract(end_ec, day.water_content),
    )


def simulate_salt(
    days: Iterable[DailySaltInput],
    start_ec_solution_ds_per_m: float,
    root_zone_depth_m: float,
    salt_uptake_per_n: float = SALT_UPTAKE_PER_N_DS_PER_M_MM,
    substeps: int = 24,
) -> list[DailySalt]:
    """作期を通して塩類収支を回す。

    Args:
        days: 日ごとの入力。日付の順に並んでいること。
        start_ec_solution_ds_per_m: 作付け前の土壌溶液EC [dS/m]。
            ★土壌診断の EC(1:5) から ec_1to5_to_soil_solution で出す。
        root_zone_depth_m: 根群域の深さ [m]
        salt_uptake_per_n: 吸収N 1 g あたり持ち出される塩
        substeps: 1日の刻み数
    """
    rows = list(days)
    if not rows:
        raise ValueError("入力の日が1日もない。")
    if start_ec_solution_ds_per_m < 0.0:
        raise ValueError(
            f"作付け前のECが負: {start_ec_solution_ds_per_m} dS/m"
        )

    first_water_mm = rows[0].water_content * root_zone_depth_m * MM_PER_M
    salt = start_ec_solution_ds_per_m * first_water_mm

    results = []
    for day in rows:
        result = step_one_day_salt(
            salt, day, root_zone_depth_m, salt_uptake_per_n, substeps)
        results.append(result)
        salt = result.end_salt
    return results


def steady_state_concentration_factor(leaching_fraction: float) -> float:
    """定常状態での「根圏EC ÷ 給液EC」。

        濃縮倍率 = 1 ÷ 流亡率

    【なぜこうなるか】
    入ってくる塩 = 潅水量 × 給液EC。出ていく塩 = 流亡量 × 根圏EC。
    落ち着いた状態では両者が等しいので

        潅水 × 給液EC = 流亡 × 根圏EC
        根圏EC ÷ 給液EC = 潅水 ÷ 流亡 = 1 ÷ 流亡率

    【桁の感覚】
        流亡率 50%（春の運用）→ 2.0 倍
        流亡率 20%           → 5.0 倍
        流亡率 10%           → 10.0 倍

    ★作物が塩を吸うぶんを数えていないので、これは**上限**。
      実際にはこれより低く落ち着く。向きを押さえるための式。
    """
    if not 0.0 < leaching_fraction <= 1.0:
        raise ValueError(
            f"流亡率は 0 より大きく 1 以下。"
            f"渡された値: {leaching_fraction}。"
            f"流亡率ゼロでは塩が出ていく先がなく、定常にならない。"
        )
    return 1.0 / leaching_fraction
