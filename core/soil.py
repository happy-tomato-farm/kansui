"""
土壌の水分特性と、根からの吸水の上限を計算する。
（生物環境物理学の基礎 第9章「蒸散と植物による水分吸収」）

【このモジュールが答えようとしている問い】
1. 潅水したうち、どれだけが根群域より下へ抜けていくのか
2. 土が乾いてくると、蒸散はどこから抑えられ始めるのか
3. 「春先に圃場容水量以上の水を入れ続ける」のは物理的に意味があるのか

【全体の流れ】

    含水率 θ  ──(水分特性曲線)──→  水ポテンシャル ψ
        │                                  │
        │                                  └─→ pF、空気率
        │
        ├──(Mualem の式)──→ 透水係数 K ──→ 重力排水（地下への流亡）
        │
        └──(式9.20)──→ 有効水分度 Aw ──(式9.22)──→ 可能吸水速度 Up

【符号の約束】
水ポテンシャルは本来マイナスの値をとる（土は水を引き止めている）。
このモジュールでは、式が読みやすくなるよう **大きさ（正の値）** で扱う。
「ψ = 33 J/kg」と書いたら、実際の水ポテンシャルは −33 J/kg のことだ。

【単位】
  含水率 θ        : m³/m³（体積含水率。土1 m³ 中の水の体積）
  水ポテンシャル ψ : J/kg（単位質量あたり。1 J/kg = 1 kPa = 10.2 cm 水柱）
  透水係数 K      : mm/日（水深に直した速さ。潅水量と直接比べられる）

【この章の重要な前提（P147）】
教科書は「吸水に関わる抵抗のほとんどすべては根に存在する
（すなわち土壌の抵抗は無視できる）」という立場をとる。
したがってここでは、土の中を水が根に向かって流れる抵抗（Gardner の
円筒流など）は計算しない。土壌水分は「根の周りの水ポテンシャルを
決めるもの」としてだけ効く。
"""

from __future__ import annotations

from dataclasses import dataclass

from config import (
    CAMPBELL_UPTAKE_COEF,
    CAMPBELL_UPTAKE_EXPONENT,
    FIELD_CAPACITY_POTENTIAL_J_PER_KG,
    MIN_AIR_FILLED_POROSITY,
    MUALEM_TORTUOSITY_L,
    NEGLIGIBLE_DRAINAGE_MM_PER_DAY,
    SOIL_RETENTION,
    SOIL_SATURATED_CONDUCTIVITY_MM_PER_DAY,
    WILTING_POINT_POTENTIAL_J_PER_KG,
)

# 1 J/kg が何 cm 水柱にあたるか。pF への換算に使う。
#   ψ[J/kg] × 水の密度 1000 kg/m³ = 圧力[Pa]
#   1000 Pa ÷ (1000 kg/m³ × 9.80665 m/s²) = 0.10197 m = 10.197 cm
#
# 逆向き（cm → J/kg）は 0.0980665 倍。土性の調査.xlsx の ψm 列は
# 表記が J/kg でも中身は cm なので、読み込むときはこの逆換算が要る。
# 詳しくは config.py の第7節の注意書きを参照。
CM_WATER_PER_J_PER_KG = 10.197
J_PER_KG_PER_CM_WATER = 1.0 / CM_WATER_PER_J_PER_KG

# 水ポテンシャルの探索範囲 [J/kg]。含水率から逆算するときに使う。
_PSI_SEARCH_MIN = 1.0e-4
_PSI_SEARCH_MAX = 1.0e7


# =============================================================================
# 1. 水分特性曲線（含水率 ⇔ 水ポテンシャル）
# =============================================================================

def _effective_saturation_components(potential_j_per_kg: float) -> tuple[float, float]:
    """二峰モデルの、大きい孔隙・小さい孔隙それぞれの飽和度を返す。

    透水係数の式でも同じ2つが必要になるので、共通の下請けにしてある。
    """
    psi = max(potential_j_per_kg, 0.0)

    n_coarse = SOIL_RETENTION["N_COARSE"]
    n_fine = SOIL_RETENTION["N_FINE"]
    # van Genuchten の慣例で m は n から決まる（Mualem の式が解けるようにするため）
    m_coarse = 1.0 - 1.0 / n_coarse
    m_fine = 1.0 - 1.0 / n_fine

    se_coarse = (1.0 + (SOIL_RETENTION["ALPHA_COARSE"] * psi) ** n_coarse) ** (-m_coarse)
    se_fine = (1.0 + (SOIL_RETENTION["ALPHA_FINE"] * psi) ** n_fine) ** (-m_fine)
    return se_coarse, se_fine


def effective_saturation(potential_j_per_kg: float) -> float:
    """有効飽和度 Se（0〜1）を水ポテンシャルから求める。

    Se = 1 なら孔隙が水で満たされている、Se = 0 なら残留含水率まで乾いている。
    2つの孔隙系の重み付き平均をとる。
    """
    se_coarse, se_fine = _effective_saturation_components(potential_j_per_kg)
    weight = SOIL_RETENTION["WEIGHT_COARSE"]
    return weight * se_coarse + (1.0 - weight) * se_fine


def water_content_from_potential(potential_j_per_kg: float) -> float:
    """水ポテンシャル [J/kg] → 体積含水率 [m³/m³]。

    θ = θr + (θs − θr) · Se
    """
    theta_r = SOIL_RETENTION["THETA_R"]
    theta_s = SOIL_RETENTION["THETA_S"]
    return theta_r + (theta_s - theta_r) * effective_saturation(potential_j_per_kg)


def potential_from_water_content(water_content: float) -> float:
    """体積含水率 [m³/m³] → 水ポテンシャルの大きさ [J/kg]。

    水分特性曲線には逆関数の式がないので、二分法で探す。
    θ が大きいほど ψ は小さい（単調減少）ので、解は必ず1つに決まる。

    飽和以上を渡されたら 0（＝水を引き止める力がない）を返す。
    残留含水率以下は物理的にありえないのでエラーにする。
    """
    theta_r = SOIL_RETENTION["THETA_R"]
    theta_s = SOIL_RETENTION["THETA_S"]

    if water_content >= theta_s:
        return 0.0
    if water_content <= theta_r:
        raise ValueError(
            f"体積含水率 {water_content:.4f} が残留含水率 {theta_r:.4f} 以下。"
            f"この土ではここまで乾くことはないので、入力値か水収支の計算を疑うことだ。"
        )

    low, high = _PSI_SEARCH_MIN, _PSI_SEARCH_MAX
    # 60回も割れば幅は 10^7 / 2^60 まで縮むので、精度は十分すぎるほど足りる
    for _ in range(100):
        mid = (low + high) / 2.0
        if water_content_from_potential(mid) > water_content:
            # まだ湿りすぎ → もっと強く吸う側（ψ を大きく）へ
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def pf_from_potential(potential_j_per_kg: float) -> float:
    """水ポテンシャルの大きさ [J/kg] → pF値。

    pF は「水柱の高さ [cm] の常用対数」。日本の土壌診断で使われる表し方。
    """
    if potential_j_per_kg <= 0.0:
        return 0.0
    import math

    return math.log10(potential_j_per_kg * CM_WATER_PER_J_PER_KG)


def potential_from_pf(pf: float) -> float:
    """pF値 → 水ポテンシャルの大きさ [J/kg]。"""
    return 10.0 ** pf / CM_WATER_PER_J_PER_KG


def air_filled_porosity(water_content: float) -> float:
    """空気で満たされている孔隙の割合 [m³/m³]。

    根は土の中の空気から酸素をとって呼吸している。
    これが 10% を下回ると根の活性が落ち、さらに低いと根腐れに向かう。
    水をやりすぎたときに何が起きるかは、この値を見れば分かる。
    """
    return max(SOIL_RETENTION["THETA_S"] - water_content, 0.0)


# =============================================================================
# 2. 透水係数と重力排水（地下への流亡）
# =============================================================================

def relative_conductivity(potential_j_per_kg: float) -> float:
    """相対透水係数 Kr（0〜1）。飽和時を1としたときの水の通しやすさ。

    【考え方（Mualem の式）】
    水は太い管ほど速く流れる。土が乾くと太い孔隙から順に空になるので、
    残った細い孔隙だけで水を運ぶことになり、透水係数は急激に落ちる。
    Mualem はこれを水分特性曲線の形から計算する方法を示した。

    【二峰モデル版（Durner 1994）】
    大きい孔隙系と小さい孔隙系で別々に計算して、孔隙の大きさ（α に比例）で
    重みをつけて足し合わせる。

        Kr = Se^L · [ Σ wi·αi·(1 − (1 − Sei^(1/mi))^mi) ]² / [ Σ wi·αi ]²

    Se^L の L は「水の通り道の曲がりくねり具合」を表す補正で、慣用値 0.5。
    """
    se_coarse, se_fine = _effective_saturation_components(potential_j_per_kg)

    weight = SOIL_RETENTION["WEIGHT_COARSE"]
    alpha_coarse = SOIL_RETENTION["ALPHA_COARSE"]
    alpha_fine = SOIL_RETENTION["ALPHA_FINE"]
    m_coarse = 1.0 - 1.0 / SOIL_RETENTION["N_COARSE"]
    m_fine = 1.0 - 1.0 / SOIL_RETENTION["N_FINE"]

    se_total = weight * se_coarse + (1.0 - weight) * se_fine

    # Mualem の積分を解いた結果の項。各孔隙系がどれだけ水を通すかを表す。
    term_coarse = 1.0 - (1.0 - se_coarse ** (1.0 / m_coarse)) ** m_coarse
    term_fine = 1.0 - (1.0 - se_fine ** (1.0 / m_fine)) ** m_fine

    numerator = weight * alpha_coarse * term_coarse + (1.0 - weight) * alpha_fine * term_fine
    denominator = weight * alpha_coarse + (1.0 - weight) * alpha_fine

    return se_total ** MUALEM_TORTUOSITY_L * (numerator / denominator) ** 2


def hydraulic_conductivity_mm_per_day(water_content: float) -> float:
    """体積含水率 [m³/m³] → 不飽和透水係数 [mm/日]。"""
    potential = potential_from_water_content(water_content)
    return SOIL_SATURATED_CONDUCTIVITY_MM_PER_DAY * relative_conductivity(potential)


def drainage_flux_mm_per_day(water_content: float) -> float:
    """重力による下向きの排水速度 [mm/日]。＝ 地下への流亡速度。

    【単位勾配の仮定】
    土の中の水の流れは
        流束 = −K · (dψ/dz + 重力)
    で決まる。潅水の後しばらくすると、深さ方向の水分分布がならされて
    dψ/dz ≒ 0 になり、重力の項だけが残る。すると

        流束 = K（その含水率での透水係数）

    という単純な形になる。これを単位勾配の仮定と呼ぶ。
    根群域より下の排水を見積もるときの標準的なやり方で、
    潅水の何時間か後からはよく成り立つ。

    【この仮定が崩れる場面】
    ・潅水した直後の数十分（湿潤前線が下りている最中）
    ・地下水位が根群域に近いとき（下から押し戻される）
    """
    return hydraulic_conductivity_mm_per_day(water_content)


def _solve_water_content_for_flux(target_flux_mm_per_day: float) -> float:
    """排水速度が目標値になる含水率を二分法で探す。圃場容水量の逆算に使う。"""
    low_theta = SOIL_RETENTION["THETA_R"] + 1e-6   # 乾いている側（排水はほぼ0）
    high_theta = SOIL_RETENTION["THETA_S"] - 1e-9  # 湿っている側（排水が速い）

    for _ in range(100):
        mid = (low_theta + high_theta) / 2.0
        if hydraulic_conductivity_mm_per_day(mid) > target_flux_mm_per_day:
            high_theta = mid
        else:
            low_theta = mid
    return (low_theta + high_theta) / 2.0


# =============================================================================
# 3. 圃場容水量・しおれ点（このモジュールを読み込んだ時点で確定する定数）
# =============================================================================

#: 圃場容水量の体積含水率 [m³/m³]。ψ = 33 J/kg（pF 2.53）に対応。
FIELD_CAPACITY_WATER_CONTENT = water_content_from_potential(
    FIELD_CAPACITY_POTENTIAL_J_PER_KG
)

#: 永久しおれ点の体積含水率 [m³/m³]。ψ = 1500 J/kg（pF 4.18）に対応。
WILTING_POINT_WATER_CONTENT = water_content_from_potential(
    WILTING_POINT_POTENTIAL_J_PER_KG
)

#: 排水が実質止まる含水率 [m³/m³]。上の圃場容水量が妥当かを点検するための値。
DRAINAGE_STOP_WATER_CONTENT = _solve_water_content_for_flux(
    NEGLIGIBLE_DRAINAGE_MM_PER_DAY
)


# =============================================================================
# 4. 根からの吸水（式9.20・式9.22）
# =============================================================================

def available_water_fraction(water_content: float) -> float:
    """【式9.20】有効水分度 Aw（0〜1）。

        Aw = (θ − θpwp) / (θfc − θpwp)

    植物が使える水を 0〜1 に直したものだ。
      Aw = 1 … 圃場容水量（使える水が満タン）
      Aw = 0 … 永久しおれ点（もう吸えない）

    圃場容水量を超えて水を入れても Aw は 1 で頭打ちになる。
    これが「余分な潅水は蒸散を増やさない」ことの、式の上での理由。
    """
    span = FIELD_CAPACITY_WATER_CONTENT - WILTING_POINT_WATER_CONTENT
    fraction = (water_content - WILTING_POINT_WATER_CONTENT) / span
    return min(max(fraction, 0.0), 1.0)


def potential_uptake_ratio(available_water: float) -> float:
    """【式9.22】無次元の可能吸水速度 Up*（0〜1）。

        Up* = 1 − (1 + 1.37·Aw)^(−5)

    根がその土壌水分のもとで吸い上げられる水の、最大値に対する割合。

    【この式の形が示していること】
    Aw が 1 から減っても、しばらくは Up* がほとんど落ちない。

        Aw = 1.0 → 0.987     Aw = 0.5 → 0.926     Aw = 0.3 → 0.821
        Aw = 0.2 → 0.702     Aw = 0.1 → 0.474     Aw = 0.0 → 0

    つまり **土の水分を圃場容水量いっぱいに保っても、半分まで使わせた
    場合と比べて吸水能力は6%しか変わらない**。
    一方で空気率は大きく変わる。ここが潅水設計の勘所。
    """
    ratio = 1.0 - (1.0 + CAMPBELL_UPTAKE_COEF * available_water) ** (-CAMPBELL_UPTAKE_EXPONENT)
    return min(max(ratio, 0.0), 1.0)


def uptake_limit_ratio_from_water_content(water_content: float) -> float:
    """体積含水率から、可能吸水速度の割合 Up*（0〜1）を一気に求める。"""
    return potential_uptake_ratio(available_water_fraction(water_content))


def soil_to_leaf_driving_force(
    water_content: float,
    leaf_potential_j_per_kg: float = 1000.0,
) -> float:
    """土から葉へ水を引き上げる駆動力 [J/kg]。

        駆動力 = |ψ葉| − |ψ土|

    【これが答える問い】
    「土を湿らせておけば、水ポテンシャルの差が大きくなって蒸散が増えるのでは？」

    答えは、差はほとんど変わらない、。

    日中のトマトの葉内水ポテンシャルは −1.0 MPa（＝ 1000 J/kg）前後ある。
    それに対して土の側は、飽和でも 0、圃場容水量（pF 1.8）でわずか 6 J/kg。
    土をどれだけ湿らせても、差の分母にあたる 1000 のうち数しか動かない。

        飽和            ψ土 =    0  → 駆動力 1000（100.0%）
        圃場容水量pF1.8 ψ土 =  6.2  → 駆動力  994（ 99.4%）
        pF 2.5          ψ土 =   33  → 駆動力  967（ 96.7%）
        pF 3.0          ψ土 =   98  → 駆動力  902（ 90.2%）

    つまり圃場容水量を超えて水を入れても、駆動力は 0.6% しか増えない。
    一方で空気率は 8.2% から 0% まで落ちる。ここが判断の分かれ目。

    【果実の肥大では話が少し変わる】
    果実の水ポテンシャルは −0.35〜−0.75 MPa 程度で、葉ほど低くない。
    分母が小さいぶん土の側の変化が相対的に効き、圃場容水量から飽和まで
    湿らせたときの駆動力の増加は 1.5〜2% 程度になる。
    葉より効くとはいえ、やはり小さい。

    果実への水の入りやすさを本当に左右するのは土壌溶液の塩類濃度（EC）で、
    EC 1 dS/m がおよそ 36 J/kg（= pF 2.6 相当）に効く。マトリック
    ポテンシャルより一桁大きい。潅水で塩を洗い流すほうが、
    水分を高く保つより効き目がある、ということだ。

    引数の leaf_potential_j_per_kg は「大きさ」で渡す（1000 なら −1.0 MPa）。
    """
    soil_potential = potential_from_water_content(water_content)
    return leaf_potential_j_per_kg - soil_potential


def stomatal_water_potential_factor(leaf_potential_j_per_kg: float) -> float:
    """葉内水ポテンシャルによる気孔の閉じ具合（0〜1）。

        係数 = 1 / [1 + (ψL / ψc)^n]

    ★式番号を原典で確認できていない（config.py の注意書きを参照）。
      日単位の収支では 蒸散 = min(E_pmax, U_p) が同じ働きをするので、
      現在の計算経路では使っていない。時間単位に踏み込むときに使う。
    """
    from config import LEAF_WATER_POTENTIAL_PARAMS

    psi_leaf = max(leaf_potential_j_per_kg, 0.0)
    half = LEAF_WATER_POTENTIAL_PARAMS["HALF_CLOSURE_J_PER_KG"]
    shape = LEAF_WATER_POTENTIAL_PARAMS["SHAPE_EXPONENT"]
    return 1.0 / (1.0 + (psi_leaf / half) ** shape)


# =============================================================================
# 5. まとめて見るための入れ物
# =============================================================================

@dataclass(frozen=True)
class SoilState:
    """ある瞬間の土壌の状態。すべて計算で導かれる値をまとめて持つ。

    frozen=True にしてあるので、あとから書き換えられない。
    水収支の各日の記録として積み上げても、過去の値が壊れない。
    """

    water_content: float              # 体積含水率 [m³/m³]
    potential_j_per_kg: float         # 水ポテンシャルの大きさ [J/kg]
    pf: float                         # pF値
    air_filled_porosity: float        # 空気率 [m³/m³]
    conductivity_mm_per_day: float    # 不飽和透水係数 [mm/日]
    available_water: float            # 有効水分度 Aw（0〜1）
    uptake_ratio: float               # 可能吸水速度の割合 Up*（0〜1）

    @property
    def is_anoxic(self) -> bool:
        """根が酸欠になる領域かどうか。"""
        return self.air_filled_porosity < MIN_AIR_FILLED_POROSITY

    @property
    def is_above_field_capacity(self) -> bool:
        """圃場容水量を超えているかどうか（＝重力排水が続いている状態）。"""
        return self.water_content > FIELD_CAPACITY_WATER_CONTENT


def describe_soil_state(water_content: float) -> SoilState:
    """体積含水率から、土壌の状態を一式まとめて計算する。"""
    potential = potential_from_water_content(water_content)
    available = available_water_fraction(water_content)
    return SoilState(
        water_content=water_content,
        potential_j_per_kg=potential,
        pf=pf_from_potential(potential),
        air_filled_porosity=air_filled_porosity(water_content),
        conductivity_mm_per_day=SOIL_SATURATED_CONDUCTIVITY_MM_PER_DAY
        * relative_conductivity(potential),
        available_water=available,
        uptake_ratio=potential_uptake_ratio(available),
    )
