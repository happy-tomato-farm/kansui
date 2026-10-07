"""乾物の生産量から、作物が吸った窒素（N）の量を出す。

【このモジュールの役割】
窒素は乾物の構成要素として体に入る。タンパク質（とくに Rubisco）・核酸・
クロロフィルがN化合物だから、乾物を作れば必ずNが要る。

    N需要 = Σ( 器官の乾物 × その器官のN濃度 )

器官ごとに分けるのが肝心だ。葉のN濃度は果実の2倍以上ある（Rubisco が
葉タンパク質の2〜3割を占める）。果実は糖と水が主体でNは低い。
だから **乾物を果実に多く回す群落ほど、乾物あたりのN需要は低い**。
全身を1つのN濃度でまとめると、この違いが消える。

【★入口は2つある。どちらを使うかが設計の分かれ目】

    (A) 光合成モデルの糖から   … 日ごとの形は出るが、絶対値が信用できない
    (B) 実収量から逆算          … 絶対値は確かだが、作期の合計しか出ない

光合成モデルは実収量から逆算すると約3倍足りない（README「光合成について」）。
(A) を素で使うとN需要が1/3で出る。これは「少しずれる」ではなく、
施肥設計として使えない。

絶対値には (B) を使う。

【★日ごとの配分に糖を使ってはいけない（2026-10-01 に訂正）】
当初は「絶対値は (B)、日ごとの配分は (A)」という分担を考えていた。

    その日のN需要 = 作期のN需要(B) × ( その日の糖(A) ÷ 作期の糖の合計(A) )

較正係数が分子と分母で打ち消し合うから成り立つ——と書いていたが、
**これは間違いだった**。打ち消し合うのは誤差が一定倍率のときだけで、
実際の誤差は LAI に依存する。

月ごとの糖の平均 [g/m²/日] と、そのときの LAI:

    月      12     1     2     3     4     5     6
    糖    5.27  5.68  5.27  3.80  0.81 -0.61 -3.49
    LAI   2.89  4.02  5.08  6.24  7.19  7.19  6.40

**5月・6月は平均が負**になる。暗呼吸は LAI に比例するのに総光合成は
過小なので、LAI が大きいと呼吸が勝ってしまう。配分に使うと春の需要が
消えてしまい、話が逆になる。

README の「実収量から逆算すると約3倍足りない」は作期合計の話で、
この LAI 依存性を隠していた。

【代わりに使うもの: 受光量】

    重み = 日射 [MJ/m²] × ( 1 − exp( −k × LAI ) )

群落が実際に受け取った光のエネルギー。光利用効率（LUE）の考え方で、
乾物生産がこれにほぼ比例することはよく確かめられている。
使うのは**実測の日射と作業計画の葉枚数だけ**で、光合成モデルを通らない。
tools/時期別の妥当性を検証.py がこれを使っている。

【★単位の便利な一致】
    1 g-N/m² = 1 kg-N/10a
施肥の実績は kg/10a で記録されているので、モデル側を g/m² で持てば
換算なしで突き合わせられる。このモジュールは g-N/m² で通す。

【前提が崩れるところ】
・果実分配率を固定している。定植直後は実際には果実へほとんど回らないので、
  作期の序盤ではN需要を過大に見積もる。作期合計では問題にならない。
・日単位では合わない。植物は硝酸を液胞に貯めるので吸収と消費の時刻が
  ずれる。7〜10日の移動平均で使うこと。
"""
from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass
from typing import Iterable

from config import (
    EXTINCTION_COEFFICIENT_K,
    FERTILIZER_N_EFFICIENCY,
    LEAF_AREA_PER_LEAF_M2,
    ORGAN_N_CONTENT,
    SOIL_N_SUPPLY_KG_PER_10A,
    VEGETATIVE_PARTITION,
    YIELD_CONVERSION,
)
from core.advisor import normal_radiation_mj
from core.canopy import leaf_count_per_m2

#: 定植の ISO 週。config.LEAF_COUNT_PLAN_PER_M2 が値を持ち始める週。
#: これより前は作物がいないので、N需要をゼロとして扱う。
PLANTING_ISO_WEEK = 38

#: 1 g/m² を kg/10a に直す係数。1 なのだが、式の中で意図を示すために置く。
G_PER_M2_TO_KG_PER_10A = 1.0

#: 器官の名前。表示と内訳の順番をそろえるために持つ。
ORGANS = ("FRUIT", "LEAF", "STEM", "ROOT")

ORGAN_LABELS = {
    "FRUIT": "果実",
    "LEAF": "葉",
    "STEM": "茎",
    "ROOT": "根",
}


# =============================================================================
# 1. 前提（N濃度と分配）
# =============================================================================

@dataclass(frozen=True)
class NitrogenSettings:
    """N需要を出すときの前提。

    どれも実測していない文献値なので、感度を見られるように外から
    差し替えられる形にしてある。frozen=True なので変えたいときは
    dataclasses.replace() で新しい設定を作る。
    """

    #: 器官別のN濃度 [g-N/g-乾物]
    fruit_n: float = ORGAN_N_CONTENT["FRUIT"]
    leaf_n: float = ORGAN_N_CONTENT["LEAF"]
    stem_n: float = ORGAN_N_CONTENT["STEM"]
    root_n: float = ORGAN_N_CONTENT["ROOT"]

    #: 全乾物のうち果実に回る割合（第8節 FRUIT_ALLOCATION と同じ量）
    fruit_allocation: float = YIELD_CONVERSION["FRUIT_ALLOCATION"]

    #: 栄養器官（果実以外）の内訳。合計1になること。
    leaf_share: float = VEGETATIVE_PARTITION["LEAF"]
    stem_share: float = VEGETATIVE_PARTITION["STEM"]
    root_share: float = VEGETATIVE_PARTITION["ROOT"]

    #: 糖から乾物への変換効率（第8節と同じ量）
    sugar_to_dry_matter: float = YIELD_CONVERSION["SUGAR_TO_DRY_MATTER"]

    #: 果実の乾物率（第8節と同じ量）
    fruit_dry_matter_content: float = YIELD_CONVERSION["FRUIT_DRY_MATTER_CONTENT"]

    def __post_init__(self) -> None:
        for name, value in (
            ("fruit_n", self.fruit_n), ("leaf_n", self.leaf_n),
            ("stem_n", self.stem_n), ("root_n", self.root_n),
        ):
            if not 0.0 < value < 0.20:
                raise ValueError(
                    f"{name}（N濃度）は 0 より大きく 0.20 未満でなければならない。"
                    f"渡された値: {value}。乾物あたりの割合なので、"
                    f"4% なら 0.04 と書く（4 ではない）。"
                )
        if not 0.0 < self.fruit_allocation <= 1.0:
            raise ValueError(
                f"fruit_allocation は 0 より大きく 1 以下。"
                f"渡された値: {self.fruit_allocation}"
            )
        if not 0.0 < self.fruit_dry_matter_content <= 1.0:
            raise ValueError(
                f"fruit_dry_matter_content は 0 より大きく 1 以下。"
                f"渡された値: {self.fruit_dry_matter_content}"
            )
        if not 0.0 < self.sugar_to_dry_matter <= 1.0:
            raise ValueError(
                f"sugar_to_dry_matter は 0 より大きく 1 以下。"
                f"渡された値: {self.sugar_to_dry_matter}"
            )
        share_total = self.leaf_share + self.stem_share + self.root_share
        if abs(share_total - 1.0) > 1e-6:
            raise ValueError(
                f"栄養器官の内訳（葉・茎・根）の合計が1にならない: "
                f"{self.leaf_share} + {self.stem_share} + {self.root_share} "
                f"= {share_total}。割合なので合計1にすること。"
            )

    # --- 乾物の分配率（全乾物を1としたときの各器官の取り分）---

    @property
    def dry_matter_partition(self) -> dict[str, float]:
        """全乾物に対する器官別の分配率。合計1になる。"""
        vegetative = 1.0 - self.fruit_allocation
        return {
            "FRUIT": self.fruit_allocation,
            "LEAF": vegetative * self.leaf_share,
            "STEM": vegetative * self.stem_share,
            "ROOT": vegetative * self.root_share,
        }

    @property
    def organ_n_content(self) -> dict[str, float]:
        """器官別のN濃度 [g-N/g-乾物]。"""
        return {
            "FRUIT": self.fruit_n,
            "LEAF": self.leaf_n,
            "STEM": self.stem_n,
            "ROOT": self.root_n,
        }

    @property
    def whole_plant_n_content(self) -> float:
        """全乾物あたりの平均N濃度 [g-N/g-乾物]。

        器官別の分配率とN濃度の積の和。これ1つで
        「乾物1kg作るのにN何g要るか」が決まる。
        """
        partition = self.dry_matter_partition
        content = self.organ_n_content
        return sum(partition[organ] * content[organ] for organ in ORGANS)

    @property
    def n_per_tonne_fruit_kg(self) -> float:
        """果実1トン（生重）あたりの吸収N量 [kg-N/t]。

        ★独立した答え合わせに使う指標。施設トマトで最もよく引用される。
        文献の幅は 2.0〜3.0 kg-N/t。

            全乾物 = 果実生重 × 果実乾物率 ÷ 果実分配率
            吸収N  = 全乾物 × 全身平均N濃度
        """
        dry_matter_per_kg_fresh = (
            self.fruit_dry_matter_content / self.fruit_allocation)
        return self.whole_plant_n_content * dry_matter_per_kg_fresh * 1000.0

    def describe(self) -> str:
        """前提を日本語で説明する。"""
        partition = self.dry_matter_partition
        content = self.organ_n_content
        lines = [
            f"{'器官':<6}{'乾物の分配':>12}{'N濃度':>10}{'全乾物あたりN':>16}",
        ]
        for organ in ORGANS:
            share = partition[organ] * content[organ]
            lines.append(
                f"{ORGAN_LABELS[organ]:<6}{partition[organ] * 100:>11.1f}%"
                f"{content[organ] * 100:>9.1f}%{share * 100:>15.3f}%"
            )
        lines.append(
            f"{'合計':<6}{100.0:>11.1f}%{'—':>10}"
            f"{self.whole_plant_n_content * 100:>15.3f}%"
        )
        lines.append("")
        lines.append(
            f"果実1トンあたりの吸収N  {self.n_per_tonne_fruit_kg:.2f} kg-N/t"
            f"（文献 2.0〜3.0）"
        )
        return "\n".join(lines)


# =============================================================================
# 2. N需要（器官別の内訳つき）
# =============================================================================

@dataclass(frozen=True)
class NitrogenDemand:
    """吸収N量とその内訳。すべて床面積あたり [g-N/m²]。

    ★1 g-N/m² = 1 kg-N/10a なので、施肥実績とそのまま比べられる。
    """

    total_dry_matter_g_per_m2: float      # 全乾物 [g-DM/m²]
    by_organ_g_per_m2: dict[str, float]   # 器官別のN [g-N/m²]
    settings: NitrogenSettings

    @property
    def uptake_g_per_m2(self) -> float:
        """作物が吸ったNの合計 [g-N/m²]。"""
        return sum(self.by_organ_g_per_m2.values())

    @property
    def uptake_kg_per_10a(self) -> float:
        """同じ量を kg-N/10a で。単位が一致するので値は変わらない。"""
        return self.uptake_g_per_m2 * G_PER_M2_TO_KG_PER_10A

    @property
    def removed_g_per_m2(self) -> float:
        """果実として圃場から持ち出されるN [g-N/m²]。

        茎葉根のNは作終わりに残渣として土へ戻るので、持ち出しではない。
        次作の地力窒素になる。**吸収量と持ち出し量を混同しないこと。**
        """
        return self.by_organ_g_per_m2["FRUIT"]

    @property
    def residue_g_per_m2(self) -> float:
        """残渣として土に戻るN [g-N/m²]（茎・葉・根）。"""
        return self.uptake_g_per_m2 - self.removed_g_per_m2

    def describe(self) -> str:
        """結果を日本語の表で返す。"""
        lines = [
            f"全乾物              {self.total_dry_matter_g_per_m2 / 1000.0:8.2f} kg-DM/m²",
            "",
            f"{'器官':<6}{'N [g/m²]':>12}{'割合':>10}",
        ]
        total = self.uptake_g_per_m2
        for organ in ORGANS:
            value = self.by_organ_g_per_m2[organ]
            share = value / total * 100.0 if total > 0 else float("nan")
            lines.append(f"{ORGAN_LABELS[organ]:<6}{value:>12.2f}{share:>9.1f}%")
        lines += [
            f"{'合計':<6}{total:>12.2f}{100.0:>9.1f}%",
            "",
            f"吸収N              {self.uptake_g_per_m2:8.2f} g-N/m² "
            f"（= {self.uptake_kg_per_10a:.2f} kg-N/10a）",
            f"うち果実で持ち出し  {self.removed_g_per_m2:8.2f} g-N/m²",
            f"うち残渣で土へ戻る  {self.residue_g_per_m2:8.2f} g-N/m²",
        ]
        return "\n".join(lines)


def demand_from_dry_matter(
    total_dry_matter_g_per_m2: float,
    settings: NitrogenSettings | None = None,
) -> NitrogenDemand:
    """全乾物から吸収N量を出す。これが計算の本体。

    Args:
        total_dry_matter_g_per_m2: 全乾物 [g-DM/m²]
        settings: 前提。省くと config の既定値。
    """
    if settings is None:
        settings = NitrogenSettings()
    if total_dry_matter_g_per_m2 < 0.0:
        raise ValueError(
            f"全乾物が負: {total_dry_matter_g_per_m2} g/m²。"
            f"光合成モデルの出力か積算のしかたを疑うこと。"
        )

    partition = settings.dry_matter_partition
    content = settings.organ_n_content
    by_organ = {
        organ: total_dry_matter_g_per_m2 * partition[organ] * content[organ]
        for organ in ORGANS
    }
    return NitrogenDemand(
        total_dry_matter_g_per_m2=total_dry_matter_g_per_m2,
        by_organ_g_per_m2=by_organ,
        settings=settings,
    )


def demand_from_sugar(
    sugar_g_per_m2: float,
    settings: NitrogenSettings | None = None,
    calibration_factor: float = 1.0,
) -> NitrogenDemand:
    """光合成モデルの糖から吸収N量を出す（入口A）。

    ★★★ 実務では使わないこと ★★★
      素の値（calibration_factor=1.0）は実収量の約1/3になる。しかも
      誤差は一定倍率ではなく LAI に依存し、LAI 7 では糖が負になる
      （上の説明）。較正係数を掛けても春が救えない。
      **絶対値も日ごとの配分も demand_from_fresh_yield 側で作ること。**
      この関数は、光合成モデルを直したときに改善を測るために残してある。

    Args:
        sugar_g_per_m2: 糖の生産量 [g-糖/m²]（advisor の sugar_g_per_m2）
        settings: 前提
        calibration_factor: 光合成の較正係数（config 第8-2節）
    """
    if settings is None:
        settings = NitrogenSettings()
    if calibration_factor <= 0.0:
        raise ValueError(
            f"較正係数は正の値でなければならない。"
            f"渡された値: {calibration_factor}"
        )
    dry_matter = sugar_g_per_m2 * calibration_factor * settings.sugar_to_dry_matter
    return demand_from_dry_matter(dry_matter, settings)


def demand_from_fresh_yield(
    fresh_yield_kg_per_m2: float,
    settings: NitrogenSettings | None = None,
) -> NitrogenDemand:
    """実収量（果実生重）から吸収N量を逆算する（入口B）。★絶対値はこちらを使う。

    果実生重 → 果実乾物 → 全乾物 → 器官別N と遡る。
    光合成モデルを通らないので、較正係数3倍の不確かさが入らない。

    Args:
        fresh_yield_kg_per_m2: 作期の果実生重 [kg/m²]
        settings: 前提
    """
    if settings is None:
        settings = NitrogenSettings()
    if fresh_yield_kg_per_m2 < 0.0:
        raise ValueError(
            f"収量が負: {fresh_yield_kg_per_m2} kg/m²。"
        )
    fruit_dry_g = (
        fresh_yield_kg_per_m2 * 1000.0 * settings.fruit_dry_matter_content)
    total_dry_g = fruit_dry_g / settings.fruit_allocation
    return demand_from_dry_matter(total_dry_g, settings)


# =============================================================================
# 3. 施肥N量へ（土壌の供給を差し引き、利用率で割り戻す）
# =============================================================================

@dataclass(frozen=True)
class FertilizerPlan:
    """必要な施肥N量。すべて kg-N/10a。"""

    uptake_kg_per_10a: float          # 作物が吸う必要のあるN
    soil_supply_kg_per_10a: float     # 土が出してくれるN（地力窒素）
    shortfall_kg_per_10a: float       # 足りない分（根の前まで届けばよい量）
    efficiency: float                 # 施肥Nの利用率
    fertilizer_kg_per_10a: float      # 実際に施肥すべきN量

    def describe(self) -> str:
        return "\n".join([
            f"吸収N需要        {self.uptake_kg_per_10a:7.2f} kg-N/10a",
            f"地力窒素（土）   {self.soil_supply_kg_per_10a:7.2f} kg-N/10a",
            f"足りない分       {self.shortfall_kg_per_10a:7.2f} kg-N/10a",
            f"÷ 利用率 {self.efficiency:.2f}",
            f"施肥N量          {self.fertilizer_kg_per_10a:7.2f} kg-N/10a",
        ])


def required_fertilizer_n(
    uptake_kg_per_10a: float,
    soil_supply_kg_per_10a: float,
    efficiency: float,
) -> FertilizerPlan:
    """吸収N需要から、施肥すべきN量を出す。

    ★流亡をここで二重に数えないこと。利用率に流亡を込める設計（設計A）と、
      根圏のN収支で流亡を明示的に出す設計（設計B）のどちらかを選び、
      選んだほうだけを使う。config 第10-5節に書いてある。

    Args:
        uptake_kg_per_10a: 作物が吸う必要のあるN [kg-N/10a]
        soil_supply_kg_per_10a: 地力窒素 [kg-N/10a]
        efficiency: 施肥Nの利用率（0〜1）
    """
    if not 0.0 < efficiency <= 1.0:
        raise ValueError(
            f"利用率は 0 より大きく 1 以下でなければならない。"
            f"渡された値: {efficiency}"
        )
    shortfall = uptake_kg_per_10a - soil_supply_kg_per_10a
    return FertilizerPlan(
        uptake_kg_per_10a=uptake_kg_per_10a,
        soil_supply_kg_per_10a=soil_supply_kg_per_10a,
        shortfall_kg_per_10a=shortfall,
        efficiency=efficiency,
        # 土の供給が需要を上回れば施肥は要らない。負の施肥量は出さない。
        fertilizer_kg_per_10a=max(shortfall, 0.0) / efficiency,
    )


def _light_interception(date: dt.date) -> float:
    """その日の「葉が光を受け取る割合」（0〜1）。

        受光率 = 1 − exp( −k × LAI )

    Beer 則。LAI が大きいほど1に近づく（それ以上増えても受け取れない）。
    """
    lai = leaf_count_per_m2(date) * LEAF_AREA_PER_LEAF_M2
    return 1.0 - math.exp(-EXTINCTION_COEFFICIENT_K * lai)


def daily_light_capture(date: dt.date, radiation_mj: float) -> float:
    """その日の受光量 [MJ/m²]。乾物生産の重みに使う。

        受光量 = 日射 × 受光率

    ★光合成モデルの糖は使えない（このファイルの冒頭の説明）。
      日射は実測、葉枚数は作業計画の値なので、壊れているモデルを通らない。
    """
    if radiation_mj < 0.0:
        raise ValueError(f"日射が負: {radiation_mj} MJ/m²")
    return radiation_mj * _light_interception(date)


def season_light_capture(season_year: int) -> float:
    """作期を通した平年の受光量の合計 [MJ/m²]。

    平年の日射（`core.advisor.normal_radiation_mj`）と作業計画の葉枚数から
    組み立てる。その日の受光量をこれで割れば「作期のうち何割の日か」が出る。

    作期は8月始まり。定植前（9月中旬より前）は葉がないので数えない。

    Args:
        season_year: 作期の年。2025 なら 2025年8月〜2026年7月。
    """
    total = 0.0
    date = dt.date(season_year, 8, 1)
    end = dt.date(season_year + 1, 7, 31)
    while date <= end:
        if not _before_planting(date):
            total += daily_light_capture(
                date, normal_radiation_mj(date.timetuple().tm_yday))
        date += dt.timedelta(days=1)
    if total <= 0.0:
        raise ValueError(
            f"作期{season_year}の受光量の合計が0になった。"
            f"平年日射か葉枚数の計画を疑うこと。"
        )
    return total


def _before_planting(date: dt.date) -> bool:
    """定植前かどうか。

    作期は8月始まりだが定植は9月中旬（ISO週38）。
    config.LEAF_COUNT_PLAN_PER_M2 が値を持ち始める週に合わせている。
    """
    return date.month in (7, 8) or (
        date.month == 9 and date.isocalendar().week < PLANTING_ISO_WEEK)


@dataclass(frozen=True)
class DailyFertilizerAdvice:
    """その日の施肥のめやす。すべて kg-N/10a（= g-N/m²）。"""

    date: dt.date
    radiation_mj: float               # 入れた日射予測（ハウスセンサー基準）
    light_capture_mj: float           # その日の受光量
    season_light_capture_mj: float    # 作期の平年受光量の合計
    share: float                      # その日の取り分（作期合計に対する割合）
    season_uptake_kg_per_10a: float   # 作期の吸収N需要
    season_fertilizer_kg_per_10a: float  # 作期に施肥すべきN
    daily_n_kg_per_10a: float         # ★液肥混入機レシピに入れる数字
    before_planting: bool             # 定植前なら True（N量はゼロ）

    @property
    def daily_n_g_per_m2(self) -> float:
        """同じ量を g-N/m² で。単位が一致するので値は変わらない。"""
        return self.daily_n_kg_per_10a


def advise_fertilizer_n(
    date: dt.date,
    radiation_mj: float,
    fresh_yield_target_kg_per_m2: float,
    soil_n_supply_kg_per_10a: float = SOIL_N_SUPPLY_KG_PER_10A,
    efficiency: float = FERTILIZER_N_EFFICIENCY,
    settings: NitrogenSettings | None = None,
) -> DailyFertilizerAdvice:
    """★朝に使う: その日入れるべきN量 [kg-N/10a] を出す。

    道すじは3段。

        (1) 目標収量 → 作期の吸収N需要        （demand_from_fresh_yield）
        (2) 地力窒素を引き、利用率で割り戻す  （required_fertilizer_n）
        (3) その日の受光量の取り分で割る

    (3) の分母は**平年の受光量の合計**。その日の日射が平年より明るければ
    取り分が増え、暗ければ減る。日射に比例して潅水量を決めるのと同じ考え方だ。

    ★出口はこの1つの数字で、液肥混入機レシピの「1日のN量」にそのまま入れる。
      潅水アプリが「10MJあたり潅水量」を手渡しているのと同じ形。

    【★日ごとに合わせる必要はない】
    Nの吸収はその日の光合成と日単位では一致しない。植物は硝酸を液胞に貯める。
    7〜10日の移動平均で動かすのが現実的（このファイルの冒頭の説明）。

    Args:
        date: 日付
        radiation_mj: その日の日射予測 [MJ/m²]（ハウスセンサー基準）
        fresh_yield_target_kg_per_m2: 目標収量 [kg/m²]（作期を通した果実生重）
        soil_n_supply_kg_per_10a: 地力窒素（config 第10-4節）
        efficiency: 施肥Nの利用率（config 第10-5節。★地力窒素と対で使う）
        settings: N濃度・分配の前提
    """
    if radiation_mj < 0.0:
        raise ValueError(f"日射が負: {radiation_mj} MJ/m²")

    demand = demand_from_fresh_yield(fresh_yield_target_kg_per_m2, settings)
    plan = required_fertilizer_n(
        demand.uptake_kg_per_10a, soil_n_supply_kg_per_10a, efficiency)

    season_year = date.year if date.month >= 8 else date.year - 1
    season_total = season_light_capture(season_year)

    before = _before_planting(date)
    capture = 0.0 if before else daily_light_capture(date, radiation_mj)
    share = capture / season_total

    return DailyFertilizerAdvice(
        date=date,
        radiation_mj=radiation_mj,
        light_capture_mj=capture,
        season_light_capture_mj=season_total,
        share=share,
        season_uptake_kg_per_10a=demand.uptake_kg_per_10a,
        season_fertilizer_kg_per_10a=plan.fertilizer_kg_per_10a,
        daily_n_kg_per_10a=plan.fertilizer_kg_per_10a * share,
        before_planting=before,
    )


@dataclass(frozen=True)
class WeeklyFertilizerAdvice:
    """1週間で使うN量のめやす。すべて kg-N/10a（= g-N/m²）。

    【なぜ週で見るのか】
    Nの吸収はその日の光合成と日単位では一致しない（植物は硝酸を液胞に貯める）。
    7〜10日でならして見るのが現実的で、液肥混入機レシピ側も
    「1週間で使いたいN量」を主の欄のひとつに持っている。

    ★レシピ側の定義に合わせてある。
      レシピは `nday = nweek ÷ 7` と割るだけなので、
      こちらも**快晴が7日つづいた場合**の合計を渡す。
      そうしないと向こうで割った値と、こちらの日別の値がずれる。

    Attributes:
        start_date: 期間の初日
        end_date: 期間の最終日
        days: 日数（ふつう7）
        growing_days: そのうち定植後の日数（定植前は液肥を出さない）
        week_light_capture_mj: 期間の受光量の合計
        season_light_capture_mj: 作期の平年受光量の合計
        share: 期間の取り分（作期合計に対する割合）
        season_uptake_kg_per_10a: 作期の吸収N需要
        season_fertilizer_kg_per_10a: 作期に施肥すべきN
        week_n_kg_per_10a: ★レシピの「1週間で使いたいN量」に入れる数字
        daily_n_kg_per_10a: week_n ÷ 7。レシピ側の割り算と同じ
    """

    start_date: dt.date
    end_date: dt.date
    days: int
    growing_days: int
    week_light_capture_mj: float
    season_light_capture_mj: float
    share: float
    season_uptake_kg_per_10a: float
    season_fertilizer_kg_per_10a: float
    week_n_kg_per_10a: float
    daily_n_kg_per_10a: float

    def describe(self) -> str:
        """画面や検証スクリプトに出す一行。"""
        if self.growing_days == 0:
            return (
                f"{self.start_date}〜{self.end_date} は全日が定植前なので"
                f"液肥を出さない"
            )
        note = (
            "" if self.growing_days == self.days
            else f"（うち定植後 {self.growing_days} 日）"
        )
        return (
            f"{self.start_date}〜{self.end_date} の{self.days}日{note}で "
            f"{self.week_n_kg_per_10a:.3f} kg-N/10a"
            f"（1日あたり {self.daily_n_kg_per_10a:.4f}）"
        )


def advise_fertilizer_n_week(
    start_date: dt.date,
    radiations_mj: Iterable[float],
    fresh_yield_target_kg_per_m2: float,
    soil_n_supply_kg_per_10a: float = SOIL_N_SUPPLY_KG_PER_10A,
    efficiency: float = FERTILIZER_N_EFFICIENCY,
    settings: NitrogenSettings | None = None,
) -> WeeklyFertilizerAdvice:
    """★1週間で使うN量 [kg-N/10a] を出す。

    `advise_fertilizer_n` を日数分呼んで足すだけ。式を二重に書かないため、
    日別の計算はそちらに任せている。

    【快晴7日ぶんを渡すときの使い方】
    `radiations_mj` にその7日ぶんの**快晴日射**を入れる。
    快晴日射は日付で変わる（冬は少なく夏は多い）ので、
    1つの値を7回ではなく、日ごとの値を並べて渡すこと。

        radiations = [sensor_basis_radiation_mj(start + dt.timedelta(days=i), ...)
                      for i in range(7)]

    Args:
        start_date: 期間の初日
        radiations_mj: 初日から順に並べた日射 [MJ/m²]。長さが日数になる
        fresh_yield_target_kg_per_m2: 目標収量 [kg/m²]
        soil_n_supply_kg_per_10a: 地力窒素（config 第10-4節）
        efficiency: 施肥Nの利用率（config 第10-5節。★地力窒素と対で使う）
        settings: N濃度・分配の前提

    Raises:
        ValueError: 日射の並びが空のとき
    """
    radiations = list(radiations_mj)
    if not radiations:
        raise ValueError(
            f"日射の並びが空。{start_date} から何日ぶんかを渡すこと。"
        )

    total_n = 0.0
    total_capture = 0.0
    growing_days = 0
    first: DailyFertilizerAdvice | None = None
    for offset, radiation_mj in enumerate(radiations):
        date = start_date + dt.timedelta(days=offset)
        daily = advise_fertilizer_n(
            date, radiation_mj, fresh_yield_target_kg_per_m2,
            soil_n_supply_kg_per_10a, efficiency, settings)
        total_n += daily.daily_n_kg_per_10a
        total_capture += daily.light_capture_mj
        if not daily.before_planting:
            growing_days += 1
        if first is None:
            first = daily

    days = len(radiations)
    return WeeklyFertilizerAdvice(
        start_date=start_date,
        end_date=start_date + dt.timedelta(days=days - 1),
        days=days,
        growing_days=growing_days,
        week_light_capture_mj=total_capture,
        season_light_capture_mj=first.season_light_capture_mj,
        share=total_capture / first.season_light_capture_mj,
        season_uptake_kg_per_10a=first.season_uptake_kg_per_10a,
        season_fertilizer_kg_per_10a=first.season_fertilizer_kg_per_10a,
        week_n_kg_per_10a=total_n,
        # ★レシピ側は nweek を 7 で割る。日数が7でなくても合うように
        #   「日数で割る」ではなく「レシピと同じ割り算」をしておく。
        daily_n_kg_per_10a=total_n / days,
    )


def solve_soil_n_supply(
    uptake_kg_per_10a: float,
    fertilizer_kg_per_10a: float,
    efficiency: float,
) -> float:
    """★逆算: 実績が成り立つために土が出していたはずのN量 [kg-N/10a]。

        地力窒素 = 吸収N需要 − 施肥N × 利用率

    これが第10-4節の「唯一の未知数」。逆算した値が文献の範囲
    （施設土壌で 10〜30 kg-N/10a/作）に収まれば、前提は成り立っている。
    範囲から外れたら、N濃度・分配率・収量のどれかを疑うことになる。

    較正係数と同じ役割を果たす数字だ。
    """
    if not 0.0 < efficiency <= 1.0:
        raise ValueError(
            f"利用率は 0 より大きく 1 以下でなければならない。"
            f"渡された値: {efficiency}"
        )
    return uptake_kg_per_10a - fertilizer_kg_per_10a * efficiency
