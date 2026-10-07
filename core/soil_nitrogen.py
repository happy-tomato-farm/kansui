"""根群域の無機態窒素の収支。施肥したNがどこへ行くかを追う。

【このモジュールの役割】
`core/water_balance.py` と**まったく同じ作り**をしている。
貯水槽に水を入れて出すのと、硝酸を入れて出すのとで構造が同じだからだ。

    根圏の無機態N [g-N/m²]
      ＋ 施肥N（その日入れた分）
      ＋ 無機化N（有機物から出てくる分。地力窒素を日割りしたもの）
      ＋ 潅水の水に入っていたN
      − 作物の吸収N
      − 流亡N（潅水で根群域の下へ押し出された分）
      ＝ 翌日の残量

これを作期を通して回せば、**どの時期に根圏のNが薄くなるか**が出る。
施肥量の合計だけ合っていても、要る時期に足りなければ意味がない。

【★単位】
    1 g-N/m² = 1 kg-N/10a
施肥の実績が kg/10a で記録されているので、g/m² で通せば換算が要らない。
濃度だけは mg-N/L で持つ（液肥の濃さと直接くらべられるようにするため）。

【流亡をどう計算するか】
硝酸は土に吸着しない（陰イオンなので、むしろ土粒子に反発される）。
だから**水と一緒に動く**と考えてよい。根圏を「よく混ざった水槽」とみると、

    土壌溶液の濃度 C [g-N/L] = 無機態N量 [g/m²] ÷ 根圏の水量 [L/m²]
    流亡N = 流亡水量 × C

水が抜けるにつれて濃度も下がるので、1日ぶんをまとめて解くと

    残量 = 流亡前の残量 × exp( −流亡水量 ÷ 根圏の水量 )

という形になる（よく混ざった水槽を一定流量で流し続けたときの厳密解）。

★この「よく混ざった水槽」という見方が、このモジュールでいちばん強い仮定だ。
  実際の土では
    ・大きな孔隙を水が素通りする（バイパス流）→ 流亡は計算より少ない
    ・古い溶液が押し出される（ピストン流）  → 流亡は計算より多い
  の両方が起きる。この土は二峰性の孔隙分布（config 第7-1節）なので
  バイパス流が起きやすい側で、**流亡を過大に見積もる向き**と見ておくこと。

【無機化をどう日割りするか】
地力窒素（作期の合計）は実績から逆算した値（config 第10-4節）。
これを日ごとに配るのに、温度で重みを付ける。

    重み = Q10 ^ ( (土温 − 基準温度) ÷ 10 )

微生物の働きは温度で決まるので、**夏と秋に多く、冬に少なく**出てくる。
これは運用（冬に濃く・春に薄く施肥する）とは逆向きなので、
時期別の妥当性を見るうえで効く。

★土温は実測していない。日平均気温を代わりに使う。土は気温より
  振幅が小さく遅れるので、冬の無機化はこの計算より多め、
  夏は少なめになるはずだ。幅を見るときは Q10 を振ること。

【前提が崩れるところ】
・脱窒（嫌気条件で硝酸が窒素ガスになって抜ける）を数えていない。
  空気率が10%を下回る時間が長い日は、実際にはここで失われる。
  `core/water_balance.py` が hours_anoxic を出しているので、
  必要になったらそれを使って足せる。
・アンモニア態Nを硝酸と同じ扱いにしている。アンモニアは土に吸着して
  動きにくいが、施設の土では数日で硝酸に変わるので日単位では差が小さい。
・根群域の深さと含水率で「水槽の大きさ」が決まる。どちらも実測していない。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

from config import (
    FIELD_CAPACITY_POTENTIAL_J_PER_KG,
    ROOT_ZONE_DEPTH_M,
    SOIL_N_SUPPLY_KG_PER_10A,
)
from core.soil import water_content_from_potential

#: m → mm
MM_PER_M = 1000.0

#: g/L → mg/L
MG_PER_G = 1000.0

#: 1日を何ステップに分けるか。
#:
#: 【なぜ刻むのか】
#: 1日をまとめて計算すると、「施肥・吸収・流亡のどれを先に処理するか」で
#: 答えが変わる。先に流亡させれば入れたばかりの肥料が抜け、先に吸収させれば
#: 抜ける前に吸えてしまう。どちらを選ぶ理由もない。
#: 刻めば順番の影響が消える（24ステップで1%未満になる。
#: tests/test_soil_nitrogen.py 第3節で測ってある）。
SUBSTEPS_PER_DAY = 24


# =============================================================================
# 1. 設定
# =============================================================================

@dataclass(frozen=True)
class SoilNitrogenSettings:
    """根圏のN収支を計算するときの前提。"""

    #: 根群域の深さ [m]。水収支と同じ値を使う。
    root_zone_depth_m: float = ROOT_ZONE_DEPTH_M

    #: 根圏の体積含水率 [m³/m³]。日ごとの値を渡さないときに使う既定。
    #: 圃場容水量（pF 1.8）のときの含水率。
    water_content: float = water_content_from_potential(
        FIELD_CAPACITY_POTENTIAL_J_PER_KG)

    #: 作期を通した無機化N（地力窒素）[g-N/m² = kg-N/10a]。
    #: config 第10-4節で実績から逆算した値。
    mineralization_total_g_per_m2: float = SOIL_N_SUPPLY_KG_PER_10A

    #: 無機化の温度依存。10℃あたり何倍になるか。
    #: 土壌微生物の一般値は 2〜3。中央を採る。
    mineralization_q10: float = 2.5

    #: 温度の重みの基準 [℃]。重みの絶対値は正規化で消えるので、
    #: どこに置いても結果は変わらない。読みやすさのために20℃にしてある。
    mineralization_reference_temp_c: float = 20.0

    #: 潅水の水に入っているN [mg-N/L]。井戸水だと硝酸が入っていることがある。
    #: ★水質分析をしていないので既定は0。測ったら入れること。
    irrigation_water_n_mg_per_l: float = 0.0

    #: 吸収が半分になる土壌溶液のN濃度 [mg-N/L]。
    #: トマトの硝酸吸収の Km は 0.14〜0.7 mg-N/L と非常に低いので、
    #: 数 mg/L あればほぼ制限されない。実質「切っていない」状態の値を置く。
    #: 濃度による制限を見たいときだけ上げる。
    uptake_half_saturation_mg_per_l: float = 1.0

    def __post_init__(self) -> None:
        if self.root_zone_depth_m <= 0.0:
            raise ValueError(
                f"根群域の深さは正の値でなければならない。"
                f"渡された値: {self.root_zone_depth_m} m"
            )
        if not 0.0 < self.water_content < 1.0:
            raise ValueError(
                f"体積含水率は 0 より大きく 1 未満。"
                f"渡された値: {self.water_content}"
            )
        if self.mineralization_total_g_per_m2 < 0.0:
            raise ValueError(
                f"無機化Nが負: {self.mineralization_total_g_per_m2} g/m²。"
                f"土がNを吸い込むことはこのモデルでは扱わない。"
            )
        if self.mineralization_q10 <= 1.0:
            raise ValueError(
                f"Q10 は 1 より大きくなければならない（温度が上がれば"
                f"無機化は速くなる）。渡された値: {self.mineralization_q10}"
            )
        if self.irrigation_water_n_mg_per_l < 0.0:
            raise ValueError(
                f"潅水の水のN濃度が負: "
                f"{self.irrigation_water_n_mg_per_l} mg/L"
            )
        if self.uptake_half_saturation_mg_per_l <= 0.0:
            raise ValueError(
                f"吸収の半飽和濃度は正の値でなければならない。"
                f"渡された値: {self.uptake_half_saturation_mg_per_l}"
            )

    def root_zone_water_mm(self, water_content: float | None = None) -> float:
        """根圏が抱えている水の量 [mm] （= L/m²）。

        これが「水槽の大きさ」。濃度の分母になる。
        既定の設定（深さ0.35 m・含水率0.42）では 147 L/m²。
        """
        theta = self.water_content if water_content is None else water_content
        if not 0.0 < theta < 1.0:
            raise ValueError(
                f"体積含水率は 0 より大きく 1 未満。渡された値: {theta}"
            )
        return theta * self.root_zone_depth_m * MM_PER_M


def concentration_mg_per_l(
    nitrogen_g_per_m2: float,
    water_mm: float,
) -> float:
    """根圏の無機態N量を、土壌溶液の濃度 [mg-N/L] に直す。

    液肥の濃さ（mg-N/L）と直接くらべられる形にしておく。
    一般的な養液は 170〜220 mg-N/L。
    """
    if water_mm <= 0.0:
        raise ValueError(
            f"根圏の水量が 0 以下: {water_mm} mm。濃度を割り算で出せない。"
        )
    # g/m² ÷ (L/m²) = g/L → ×1000 で mg/L
    return nitrogen_g_per_m2 / water_mm * MG_PER_G


# =============================================================================
# 2. 無機化の日割り
# =============================================================================

def temperature_weights(
    temperatures_c: Iterable[float],
    settings: SoilNitrogenSettings | None = None,
) -> list[float]:
    """日平均気温の並びから、無機化の日ごとの重みを作る（合計1）。

        重み = Q10 ^ ( (温度 − 基準) ÷ 10 )  を正規化したもの

    暖かい日ほど大きい。基準温度をどこに置いても、正規化で消えるので
    結果は変わらない。
    """
    if settings is None:
        settings = SoilNitrogenSettings()
    raw = [
        settings.mineralization_q10
        ** ((temp - settings.mineralization_reference_temp_c) / 10.0)
        for temp in temperatures_c
    ]
    total = sum(raw)
    if total <= 0.0:
        raise ValueError(
            "温度の重みの合計が0以下になった。気温の並びが空か、"
            "値が壊れている。"
        )
    return [value / total for value in raw]


# =============================================================================
# 3. 入力と結果
# =============================================================================

@dataclass(frozen=True)
class DailyNitrogenInput:
    """1日ぶんの入力。すべて床面積あたり。"""

    date: str
    fertilizer_n_g_per_m2: float      # その日の施肥N（実績 or 設計値）
    demand_n_g_per_m2: float          # その日のN需要（core/nitrogen.py から）
    irrigation_mm: float              # 潅水量（水収支と同じ値）
    drainage_mm: float                # 流亡水量（水収支が出した値）
    mineralization_n_g_per_m2: float  # その日の無機化N（温度で日割りした値）
    water_content: float | None = None  # その日の根圏の含水率。省けば既定値

    def __post_init__(self) -> None:
        for name, value in (
            ("fertilizer_n_g_per_m2", self.fertilizer_n_g_per_m2),
            ("demand_n_g_per_m2", self.demand_n_g_per_m2),
            ("irrigation_mm", self.irrigation_mm),
            ("drainage_mm", self.drainage_mm),
            ("mineralization_n_g_per_m2", self.mineralization_n_g_per_m2),
        ):
            if value < 0.0:
                raise ValueError(
                    f"{self.date}: {name} が負（{value}）。"
                    f"入力の作り方を疑うこと。"
                )
        # ★流亡が潅水を超えることは正常。
        #   前の日に土に溜まった水が、潅水しない日に抜けていくからだ。
        #   水収支（core/water_balance.py）は貯留変化を持っているので、
        #   潅水ゼロの日にも排水が出る。ここで弾いてはいけない。
        #   上限の検査は step_one_day で根圏の水量と比べて行う。


@dataclass(frozen=True)
class DailyNitrogenBalance:
    """1日ぶんの結果。すべて床面積あたり [g-N/m² = kg-N/10a]。"""

    date: str
    start_n_g_per_m2: float           # 日の初めの根圏の無機態N
    fertilizer_n_g_per_m2: float      # 入れた施肥N
    mineralization_n_g_per_m2: float  # 土から出てきたN
    irrigation_water_n_g_per_m2: float  # 潅水の水に入っていたN
    demand_n_g_per_m2: float          # 作物が要求したN
    uptake_n_g_per_m2: float          # 実際に吸えたN
    leaching_n_g_per_m2: float        # 流亡したN
    end_n_g_per_m2: float             # 日の終わりの残量
    start_concentration_mg_per_l: float
    end_concentration_mg_per_l: float
    min_concentration_mg_per_l: float

    @property
    def shortfall_n_g_per_m2(self) -> float:
        """需要に届かなかった分 [g-N/m²]。0 より大きければN不足の日。"""
        return max(self.demand_n_g_per_m2 - self.uptake_n_g_per_m2, 0.0)

    @property
    def is_short(self) -> bool:
        """根圏のNが足りず、需要を満たせなかった日かどうか。"""
        return self.shortfall_n_g_per_m2 > 1e-9

    @property
    def leaching_fraction(self) -> float:
        """入ってきたNのうち流亡した割合。入力がゼロの日は 0。"""
        supplied = (self.fertilizer_n_g_per_m2
                    + self.mineralization_n_g_per_m2
                    + self.irrigation_water_n_g_per_m2)
        if supplied <= 0.0:
            return 0.0
        return self.leaching_n_g_per_m2 / supplied

    @property
    def balance_residual_g_per_m2(self) -> float:
        """収支の閉じ具合の点検用。0 に近くなければ計算が壊れている。"""
        return (
            self.start_n_g_per_m2
            + self.fertilizer_n_g_per_m2
            + self.mineralization_n_g_per_m2
            + self.irrigation_water_n_g_per_m2
            - self.uptake_n_g_per_m2
            - self.leaching_n_g_per_m2
            - self.end_n_g_per_m2
        )


# =============================================================================
# 4. 1日を進める
# =============================================================================

def step_one_day(
    start_n_g_per_m2: float,
    day: DailyNitrogenInput,
    settings: SoilNitrogenSettings | None = None,
    substeps: int = SUBSTEPS_PER_DAY,
) -> DailyNitrogenBalance:
    """1日ぶんのN収支を進める。

    1日を substeps に刻み、各ステップで
    「入れる → 流亡させる → 吸わせる」を繰り返す。
    刻むので、この順番の影響はほとんど残らない（第3節のテストで確認）。

    Args:
        start_n_g_per_m2: 日の初めの根圏の無機態N [g-N/m²]
        day: その日の入力
        settings: 前提
        substeps: 1日の刻み数
    """
    if settings is None:
        settings = SoilNitrogenSettings()
    if start_n_g_per_m2 < 0.0:
        raise ValueError(
            f"{day.date}: 日の初めの無機態Nが負（{start_n_g_per_m2}）。"
            f"前日の計算が壊れている。"
        )
    if substeps < 1:
        raise ValueError(f"刻み数は1以上。渡された値: {substeps}")

    water_mm = settings.root_zone_water_mm(day.water_content)

    # 1日で根圏の水を何回も入れ替えるような流亡はありえない。
    # 流亡が潅水を超えるのは正常（溜まった水が抜ける日）だが、
    # 根圏の水量を超えるなら水収支の結果か単位を疑う。
    if day.drainage_mm > water_mm:
        raise ValueError(
            f"{day.date}: 流亡水量 {day.drainage_mm:.1f} mm が"
            f"根圏の水量 {water_mm:.1f} mm を超えている。"
            f"1日で根圏の水が全部入れ替わることはない。"
            f"水収支の結果か単位を疑うこと。"
        )

    # 潅水の水に入っていたN。mg/L × L/m² ÷ 1000 = g/m²
    irrigation_water_n = (
        day.irrigation_mm * settings.irrigation_water_n_mg_per_l / MG_PER_G)

    # 1ステップぶんに割る
    fertilizer_step = day.fertilizer_n_g_per_m2 / substeps
    mineralization_step = day.mineralization_n_g_per_m2 / substeps
    water_n_step = irrigation_water_n / substeps
    demand_step = day.demand_n_g_per_m2 / substeps
    drainage_step = day.drainage_mm / substeps

    # 流亡の係数。よく混ざった水槽を一定流量で流したときの残存率。
    #   残量 = 流亡前 × exp( −流亡水量 ÷ 水槽の水量 )
    retention = math.exp(-drainage_step / water_mm)

    nitrogen = start_n_g_per_m2
    uptake_total = 0.0
    leaching_total = 0.0
    minimum = concentration_mg_per_l(nitrogen, water_mm)

    for _ in range(substeps):
        # (1) 入れる
        nitrogen += fertilizer_step + mineralization_step + water_n_step

        # (2) 流亡させる
        after_leaching = nitrogen * retention
        leaching_total += nitrogen - after_leaching
        nitrogen = after_leaching

        # (3) 吸わせる。
        #     濃度が低いと吸いにくくなる（Michaelis 型）。既定の半飽和濃度は
        #     1 mg/L と低いので、ふつうは制限がかからない。
        concentration = concentration_mg_per_l(nitrogen, water_mm)
        ratio = concentration / (
            concentration + settings.uptake_half_saturation_mg_per_l)
        # 存在する量より多くは吸えない
        uptake = min(demand_step * ratio, nitrogen)
        nitrogen -= uptake
        uptake_total += uptake

        minimum = min(minimum, concentration_mg_per_l(nitrogen, water_mm))

    return DailyNitrogenBalance(
        date=day.date,
        start_n_g_per_m2=start_n_g_per_m2,
        fertilizer_n_g_per_m2=day.fertilizer_n_g_per_m2,
        mineralization_n_g_per_m2=day.mineralization_n_g_per_m2,
        irrigation_water_n_g_per_m2=irrigation_water_n,
        demand_n_g_per_m2=day.demand_n_g_per_m2,
        uptake_n_g_per_m2=uptake_total,
        leaching_n_g_per_m2=leaching_total,
        end_n_g_per_m2=nitrogen,
        start_concentration_mg_per_l=concentration_mg_per_l(
            start_n_g_per_m2, water_mm),
        end_concentration_mg_per_l=concentration_mg_per_l(nitrogen, water_mm),
        min_concentration_mg_per_l=minimum,
    )


def simulate(
    days: Iterable[DailyNitrogenInput],
    start_n_g_per_m2: float,
    settings: SoilNitrogenSettings | None = None,
    substeps: int = SUBSTEPS_PER_DAY,
) -> list[DailyNitrogenBalance]:
    """作期を通してN収支を回す。

    Args:
        days: 日ごとの入力。日付の順に並んでいること。
        start_n_g_per_m2: 作付け前の根圏の無機態N [g-N/m²]。
            ★土壌診断の硝酸態N＋アンモニア態Nから出すのが本筋。
              測っていないあいだは、前作の残りとして見当をつける。
        settings: 前提
        substeps: 1日の刻み数
    """
    if settings is None:
        settings = SoilNitrogenSettings()
    results = []
    nitrogen = start_n_g_per_m2
    for day in days:
        result = step_one_day(nitrogen, day, settings, substeps)
        results.append(result)
        nitrogen = result.end_n_g_per_m2
    return results


# =============================================================================
# 5. まとめ
# =============================================================================

@dataclass(frozen=True)
class NitrogenBalanceSummary:
    """作期を通した集計。すべて [g-N/m² = kg-N/10a]。"""

    days: int
    fertilizer_n: float
    mineralization_n: float
    irrigation_water_n: float
    demand_n: float
    uptake_n: float
    leaching_n: float
    start_n: float
    end_n: float
    short_days: int                  # 需要を満たせなかった日数
    shortfall_n: float               # 足りなかったNの合計
    min_concentration_mg_per_l: float
    min_concentration_date: str

    @property
    def leaching_fraction(self) -> float:
        """供給したNのうち流亡した割合。"""
        supplied = (self.fertilizer_n + self.mineralization_n
                    + self.irrigation_water_n)
        if supplied <= 0.0:
            return 0.0
        return self.leaching_n / supplied

    @property
    def uptake_ratio(self) -> float:
        """需要のうち実際に吸えた割合。1.0 なら不足なし。"""
        if self.demand_n <= 0.0:
            return float("nan")
        return self.uptake_n / self.demand_n

    @property
    def apparent_efficiency(self) -> float:
        """施肥Nの見かけの利用率（吸収N ÷ 施肥N）。

        ★地力窒素も吸われているので、この値は1を超えうる。
          config 第10-5節の「利用率」とは別の量なので混同しないこと。
        """
        if self.fertilizer_n <= 0.0:
            return float("nan")
        return self.uptake_n / self.fertilizer_n

    @property
    def balance_residual(self) -> float:
        """収支の閉じ具合。0 に近くなければ計算が壊れている。"""
        return (
            self.start_n + self.fertilizer_n + self.mineralization_n
            + self.irrigation_water_n
            - self.uptake_n - self.leaching_n - self.end_n
        )

    def describe(self) -> str:
        """結果を日本語の表で返す。"""
        return "\n".join([
            f"日数                {self.days:8d} 日",
            "",
            f"{'入ってきたN':<20}",
            f"  作付け前の残り    {self.start_n:8.2f} kg-N/10a",
            f"  施肥              {self.fertilizer_n:8.2f}",
            f"  無機化（地力）    {self.mineralization_n:8.2f}",
            f"  潅水の水          {self.irrigation_water_n:8.2f}",
            f"{'出ていったN':<20}",
            f"  作物の吸収        {self.uptake_n:8.2f}",
            f"  流亡              {self.leaching_n:8.2f}"
            f"（供給の {self.leaching_fraction * 100:.0f}%）",
            f"  作終わりの残り    {self.end_n:8.2f}",
            "",
            f"N需要               {self.demand_n:8.2f} kg-N/10a",
            f"吸えた割合          {self.uptake_ratio * 100:8.1f}%",
            f"足りなかった日数    {self.short_days:8d} 日"
            f"（計 {self.shortfall_n:.2f} kg-N/10a）",
            "",
            f"いちばん薄かった濃度 {self.min_concentration_mg_per_l:7.1f} mg-N/L"
            f"（{self.min_concentration_date}）",
        ])


def summarize(
    results: Iterable[DailyNitrogenBalance],
) -> NitrogenBalanceSummary:
    """日ごとの結果を作期の集計にまとめる。"""
    rows = list(results)
    if not rows:
        raise ValueError(
            "集計する日が1日もない。入力の作り方を疑うこと。"
        )
    lowest = min(rows, key=lambda r: r.min_concentration_mg_per_l)
    return NitrogenBalanceSummary(
        days=len(rows),
        fertilizer_n=sum(r.fertilizer_n_g_per_m2 for r in rows),
        mineralization_n=sum(r.mineralization_n_g_per_m2 for r in rows),
        irrigation_water_n=sum(r.irrigation_water_n_g_per_m2 for r in rows),
        demand_n=sum(r.demand_n_g_per_m2 for r in rows),
        uptake_n=sum(r.uptake_n_g_per_m2 for r in rows),
        leaching_n=sum(r.leaching_n_g_per_m2 for r in rows),
        start_n=rows[0].start_n_g_per_m2,
        end_n=rows[-1].end_n_g_per_m2,
        short_days=sum(1 for r in rows if r.is_short),
        shortfall_n=sum(r.shortfall_n_g_per_m2 for r in rows),
        min_concentration_mg_per_l=lowest.min_concentration_mg_per_l,
        min_concentration_date=lowest.date,
    )
