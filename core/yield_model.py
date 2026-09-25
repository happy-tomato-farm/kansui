"""
光合成で作った糖を果実の生重に直し、実収量と突き合わせる。

【このモジュールの役割】
モデルの出力（糖 [kg/m²]）と、実際に穫れた量（果実生重 [kg/m²]）は
必ずずれる。そのずれを式の中に紛れ込ませず、**較正係数というひとつの
名前の付いた数字**に閉じ込めるのがここの仕事。

    案A（較正係数 = 1.0）
        物理と文献だけで組み立てた素の見積り。
        実収量は一切使っていない。残差はそのまま表に出す。

    案B（較正係数 = 実収量から逆算）
        実務の見積りに使う値。

出力は常に両方を並べる。こうしておくと、あとで過小の原因が分かって
モデルを直したとき、較正係数がその分だけ 1.0 に近づく。
つまり **較正係数そのものが「まだ説明できていない量」の目盛り**になる。

【換算の道すじ】

    糖 [kg/m²]
      │  × 糖→乾物の変換効率（生長呼吸で失われる分を引く）
      ↓
    全乾物 [kg/m²]
      │  × 果実への分配率（葉・茎・根に回る分を引く）
      ↓
    果実乾物 [kg/m²]
      │  ÷ 果実の乾物率（果実の95%は水）
      ↓
    果実生重 [kg/m²]  ← 売るものだ

【いちばん効くのはどこか】
果実の乾物率。1ポイント（5%→6%）動くと生重は約20%動く。
水収支では果実の水は蒸散の6%しかない脇役だが、収量では95%の主役になる。
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from config import (
    ACTUAL_YIELD_RECORDS,
    PHOTOSYNTHESIS_CALIBRATION_FACTOR,
    YIELD_CONVERSION,
)


# =============================================================================
# 1. 換算の設定
# =============================================================================

@dataclass(frozen=True)
class YieldConversionSettings:
    """糖から果実生重へ換算するときの前提。

    どれも実測していない値なので、感度を見られるように外から差し替えられる
    形にしてある。frozen=True なので、変えたいときは replace() で
    新しい設定を作る（元の設定は壊さない）。
    """

    sugar_to_dry_matter: float = YIELD_CONVERSION["SUGAR_TO_DRY_MATTER"]
    fruit_allocation: float = YIELD_CONVERSION["FRUIT_ALLOCATION"]
    fruit_dry_matter_content: float = YIELD_CONVERSION["FRUIT_DRY_MATTER_CONTENT"]

    #: 光合成の較正係数。1.0 が案A（素）、>1.0 が案B（実収量に合わせた値）。
    calibration_factor: float = PHOTOSYNTHESIS_CALIBRATION_FACTOR

    def __post_init__(self) -> None:
        checks = [
            ("sugar_to_dry_matter", self.sugar_to_dry_matter, 0.0, 1.0),
            ("fruit_allocation", self.fruit_allocation, 0.0, 1.0),
            ("fruit_dry_matter_content", self.fruit_dry_matter_content, 0.0, 1.0),
        ]
        for name, value, low, high in checks:
            if not low < value <= high:
                raise ValueError(
                    f"{name} は {low} より大きく {high} 以下でなければならない。"
                    f"渡された値: {value}"
                )
        if self.calibration_factor <= 0.0:
            raise ValueError(
                f"較正係数は正の値でなければならない。"
                f"渡された値: {self.calibration_factor}"
            )

    @property
    def is_calibrated(self) -> bool:
        """実収量に合わせた設定かどうか（案B なら True）。"""
        return abs(self.calibration_factor - 1.0) > 1e-9

    def describe(self) -> str:
        """設定内容を日本語で説明する。"""
        mode = "案B（実収量に較正）" if self.is_calibrated else "案A（素の値）"
        return (
            f"{mode} / 較正係数 {self.calibration_factor:.2f} / "
            f"糖→乾物 {self.sugar_to_dry_matter:.2f} / "
            f"果実分配 {self.fruit_allocation:.2f} / "
            f"果実乾物率 {self.fruit_dry_matter_content * 100:.1f}%"
        )


# =============================================================================
# 2. 換算の結果
# =============================================================================

@dataclass(frozen=True)
class YieldEstimate:
    """糖から見積もった収量。すべて床面積あたり。"""

    sugar_kg_per_m2: float             # モデルが出した糖（較正前）
    calibrated_sugar_kg_per_m2: float  # 較正係数を掛けたあとの糖
    total_dry_matter_kg_per_m2: float  # 全乾物
    fruit_dry_matter_kg_per_m2: float  # 果実の乾物
    fresh_yield_kg_per_m2: float       # 果実生重
    settings: YieldConversionSettings

    def fresh_yield_t_per_house(self, floor_area_m2: float) -> float:
        """ハウス1棟あたりの果実生重 [t]。"""
        return self.fresh_yield_kg_per_m2 * floor_area_m2 / 1000.0


def estimate_yield(
    sugar_kg_per_m2: float,
    settings: YieldConversionSettings | None = None,
) -> YieldEstimate:
    """糖の生産量から果実生重を見積もる。

    Args:
        sugar_kg_per_m2: 作期を通した糖の生産量 [kg-糖/m²]
        settings: 換算の前提。省くと config の既定（＝案A）を使う。
    """
    if settings is None:
        settings = YieldConversionSettings()
    if sugar_kg_per_m2 < 0.0:
        raise ValueError(
            f"糖の生産量が負: {sugar_kg_per_m2} kg/m²。"
            f"モデルの出力か積算のしかたを疑うことだ。"
        )

    calibrated = sugar_kg_per_m2 * settings.calibration_factor
    total_dry = calibrated * settings.sugar_to_dry_matter
    fruit_dry = total_dry * settings.fruit_allocation
    fresh = fruit_dry / settings.fruit_dry_matter_content

    return YieldEstimate(
        sugar_kg_per_m2=sugar_kg_per_m2,
        calibrated_sugar_kg_per_m2=calibrated,
        total_dry_matter_kg_per_m2=total_dry,
        fruit_dry_matter_kg_per_m2=fruit_dry,
        fresh_yield_kg_per_m2=fresh,
        settings=settings,
    )


# =============================================================================
# 3. 逆向き: 実収量から必要な同化量を出す
# =============================================================================

def required_sugar_kg_per_m2(
    fresh_yield_kg_per_m2: float,
    settings: YieldConversionSettings | None = None,
) -> float:
    """その収量を上げるのに必要だった糖の量 [kg/m²]。

    estimate_yield の逆算。較正係数は掛けない（＝実際に要る量を返す）。
    """
    if settings is None:
        settings = YieldConversionSettings()
    fruit_dry = fresh_yield_kg_per_m2 * settings.fruit_dry_matter_content
    total_dry = fruit_dry / settings.fruit_allocation
    return total_dry / settings.sugar_to_dry_matter


def solve_calibration_factor(
    model_sugar_kg_per_m2: float,
    actual_fresh_yield_kg_per_m2: float,
    settings: YieldConversionSettings | None = None,
) -> float:
    """実収量に合わせるための較正係数を逆算する（案B）。

        較正係数 = 必要だった糖 ÷ モデルが出した糖

    この値が 1.0 に近いほど、モデルが実物に近いということだ。
    """
    if model_sugar_kg_per_m2 <= 0.0:
        raise ValueError(
            f"モデルの糖生産量が 0 以下: {model_sugar_kg_per_m2} kg/m²。"
            f"較正係数を割り算で求められない。"
        )
    needed = required_sugar_kg_per_m2(actual_fresh_yield_kg_per_m2, settings)
    return needed / model_sugar_kg_per_m2


# =============================================================================
# 4. 案Aと案Bを並べて示す
# =============================================================================

@dataclass(frozen=True)
class YieldComparison:
    """案A・案B・実績を並べたものだ。"""

    plain: YieldEstimate               # 案A（較正係数 1.0）
    calibrated: YieldEstimate          # 案B（実収量に合わせたもの）
    actual_fresh_kg_per_m2: float      # 実績
    calibration_factor: float          # 案Bで使った較正係数
    floor_area_m2: float

    @property
    def shortfall_ratio(self) -> float:
        """案Aが実績の何倍にあたるか。1.0 なら残差なし。"""
        if self.actual_fresh_kg_per_m2 <= 0.0:
            return float("nan")
        return self.plain.fresh_yield_kg_per_m2 / self.actual_fresh_kg_per_m2

    def describe(self) -> str:
        """結果を日本語の表で返す。"""
        lines = [
            f"糖の生産量（モデル）      {self.plain.sugar_kg_per_m2:8.2f} kg/m²",
            "",
            f"{'':24}{'案A（素）':>14}{'案B（較正）':>14}{'実績':>12}",
            f"{'較正係数':<24}{1.0:>14.2f}{self.calibration_factor:>14.2f}{'—':>12}",
            f"{'全乾物 [kg/m²]':<24}{self.plain.total_dry_matter_kg_per_m2:>14.2f}"
            f"{self.calibrated.total_dry_matter_kg_per_m2:>14.2f}{'—':>12}",
            f"{'果実乾物 [kg/m²]':<24}{self.plain.fruit_dry_matter_kg_per_m2:>14.2f}"
            f"{self.calibrated.fruit_dry_matter_kg_per_m2:>14.2f}"
            f"{self.actual_fresh_kg_per_m2 * self.plain.settings.fruit_dry_matter_content:>12.2f}",
            f"{'果実生重 [kg/m²]':<24}{self.plain.fresh_yield_kg_per_m2:>14.2f}"
            f"{self.calibrated.fresh_yield_kg_per_m2:>14.2f}"
            f"{self.actual_fresh_kg_per_m2:>12.2f}",
            f"{'果実生重 [t/棟]':<24}"
            f"{self.plain.fresh_yield_t_per_house(self.floor_area_m2):>14.1f}"
            f"{self.calibrated.fresh_yield_t_per_house(self.floor_area_m2):>14.1f}"
            f"{self.actual_fresh_kg_per_m2 * self.floor_area_m2 / 1000.0:>12.1f}",
            "",
            f"案Aは実績の {self.shortfall_ratio * 100:.0f}%。"
            f"説明できていない不足は {self.calibration_factor:.2f} 倍。",
        ]
        return "\n".join(lines)


def compare_with_actual(
    model_sugar_kg_per_m2: float,
    house: str = "中央",
    year: int = 2026,
    settings: YieldConversionSettings | None = None,
) -> YieldComparison:
    """モデルの同化量を、記録に残っている実収量と突き合わせる。

    Args:
        model_sugar_kg_per_m2: 作期を通した糖の生産量 [kg-糖/m²]
        house: ハウス名（config.ACTUAL_YIELD_RECORDS のキー）
        year: 年度
        settings: 換算の前提。較正係数は無視され、案A・案Bの両方を組み立てる。
    """
    key = (house, year)
    if key not in ACTUAL_YIELD_RECORDS:
        raise KeyError(
            f"実収量の記録がない: {key}。"
            f"使えるのは {list(ACTUAL_YIELD_RECORDS.keys())} 。"
        )
    record = ACTUAL_YIELD_RECORDS[key]
    floor_area = record["floor_area_m2"]
    actual_fresh = record["fresh_weight_g"] / 1000.0 / floor_area

    base = settings or YieldConversionSettings()
    plain_settings = replace(base, calibration_factor=1.0)
    factor = solve_calibration_factor(model_sugar_kg_per_m2, actual_fresh, plain_settings)
    calibrated_settings = replace(base, calibration_factor=factor)

    return YieldComparison(
        plain=estimate_yield(model_sugar_kg_per_m2, plain_settings),
        calibrated=estimate_yield(model_sugar_kg_per_m2, calibrated_settings),
        actual_fresh_kg_per_m2=actual_fresh,
        calibration_factor=factor,
        floor_area_m2=floor_area,
    )
