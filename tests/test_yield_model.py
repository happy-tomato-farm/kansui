"""
収量換算と較正係数の検証。

【何を確かめるか】
1. 換算が往復して元に戻るか（糖→生重→糖）
2. 較正係数が定義どおりに効くか（案Aの何倍かがそのまま出るか）
3. 案Bが実収量にぴたりと合うか
4. 入力の検査が効くか
5. どの前提がいちばん結果を動かすか（感度）
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dataclasses import replace

from config import ACTUAL_YIELD_RECORDS, YIELD_CONVERSION_SENSITIVITY
from core.yield_model import (
    YieldConversionSettings,
    compare_with_actual,
    estimate_yield,
    required_sugar_kg_per_m2,
    solve_calibration_factor,
)

failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    mark = "OK " if condition else "NG "
    if not condition:
        failures.append(label)
    print(f"  [{mark}] {label}" + (f"  … {detail}" if detail else ""))


# =============================================================================
print("=" * 78)
print("1. 換算が往復して元に戻るか")
print("=" * 78)

settings = YieldConversionSettings()
print(f"  設定: {settings.describe()}")

worst = 0.0
for sugar in (0.5, 1.0, 1.36, 2.0, 4.0):
    estimate = estimate_yield(sugar, settings)
    back = required_sugar_kg_per_m2(estimate.fresh_yield_kg_per_m2, settings)
    worst = max(worst, abs(back - sugar))
    print(f"  糖 {sugar:.2f} kg/m² → 生重 {estimate.fresh_yield_kg_per_m2:6.2f} kg/m² "
          f"→ 逆算 {back:.4f} kg/m²")
check(worst < 1e-9, "糖 → 生重 → 糖 が元に戻る", f"最大差 {worst:.2e}")

# =============================================================================
print()
print("=" * 78)
print("2. 較正係数が定義どおりに効くか")
print("=" * 78)

base = estimate_yield(1.36, YieldConversionSettings(calibration_factor=1.0))
print(f"  較正係数  果実生重 [kg/m²]  案Aの何倍")
ratios_ok = True
for factor in (1.0, 1.5, 2.0, 2.5):
    scaled = estimate_yield(1.36, YieldConversionSettings(calibration_factor=factor))
    ratio = scaled.fresh_yield_kg_per_m2 / base.fresh_yield_kg_per_m2
    print(f"  {factor:8.2f}  {scaled.fresh_yield_kg_per_m2:14.2f}  {ratio:8.2f}")
    if abs(ratio - factor) > 1e-9:
        ratios_ok = False
check(ratios_ok, "較正係数を n 倍すると生重もちょうど n 倍になる")

check(not YieldConversionSettings(calibration_factor=1.0).is_calibrated,
      "較正係数 1.0 は案A と判定される")
check(YieldConversionSettings(calibration_factor=2.4).is_calibrated,
      "較正係数 2.4 は案B と判定される")

# =============================================================================
print()
print("=" * 78)
print("3. 案Bが実収量に合うか")
print("=" * 78)

MODEL_SUGAR = 1.36   # 実データ257日＋欠測日の外挿（scratchpad/yield_gap.py）
comparison = compare_with_actual(MODEL_SUGAR, house="中央", year=2026)
print()
print(comparison.describe())
print()

check(abs(comparison.calibrated.fresh_yield_kg_per_m2
          - comparison.actual_fresh_kg_per_m2) < 1e-9,
      "案Bの生重が実収量にぴたりと一致する",
      f"{comparison.calibrated.fresh_yield_kg_per_m2:.4f} vs "
      f"{comparison.actual_fresh_kg_per_m2:.4f}")
check(comparison.calibration_factor > 1.0,
      "較正係数は 1 より大きい（モデルは過小）",
      f"{comparison.calibration_factor:.2f} 倍")
check(abs(comparison.shortfall_ratio * comparison.calibration_factor - 1.0) < 1e-9,
      "「案Aは実績の何割」と「較正係数」が互いの逆数になっている")

# 記録にないハウス・年度はエラーにする
try:
    compare_with_actual(MODEL_SUGAR, house="西", year=2026)
    check(False, "記録にないハウスはエラーになる")
except KeyError as error:
    check("実収量の記録がない" in str(error), "記録にないハウスはエラーになる")

# =============================================================================
print()
print("=" * 78)
print("4. 入力の検査")
print("=" * 78)

for bad, label in [
    ({"fruit_dry_matter_content": 0.0}, "果実乾物率 0 を弾く"),
    ({"fruit_dry_matter_content": 1.5}, "果実乾物率 1.5 を弾く"),
    ({"fruit_allocation": -0.1}, "分配率が負なら弾く"),
    ({"calibration_factor": 0.0}, "較正係数 0 を弾く"),
]:
    try:
        YieldConversionSettings(**bad)
        check(False, label)
    except ValueError:
        check(True, label)

try:
    estimate_yield(-1.0)
    check(False, "糖が負ならエラーになる")
except ValueError as error:
    check("負" in str(error), "糖が負ならエラーになる")

try:
    solve_calibration_factor(0.0, 31.45)
    check(False, "モデルの糖が 0 なら較正係数を求めずエラーにする")
except ValueError:
    check(True, "モデルの糖が 0 なら較正係数を求めずエラーにする")

# =============================================================================
print()
print("=" * 78)
print("5. 感度: どの前提がいちばん較正係数を動かすか")
print("=" * 78)
print("  （較正係数が小さいほど、モデルが実物に近いということ）")

def sensitivity_against_actual(record: dict) -> None:
    """実収量と突き合わせて、どの前提が較正係数をいちばん動かすかを見る。

    実収量は private_data.py に置いてあり、公開リポジトリには入っていない
    （config.py 第8-3節）。記録がある環境でだけ呼ぶ。
    """
    actual = record["fresh_weight_g"] / 1000.0 / record["floor_area_m2"]

    for key, values in YIELD_CONVERSION_SENSITIVITY.items():
        print(f"\n  {key}")
        factors = []
        for value in values:
            trial = replace(YieldConversionSettings(), **{key.lower(): value},
                            calibration_factor=1.0)
            factor = solve_calibration_factor(MODEL_SUGAR, actual, trial)
            factors.append(factor)
            print(f"    {value:<8} → 較正係数 {factor:.2f} 倍")
        print(f"    幅 {min(factors):.2f}〜{max(factors):.2f}"
              f"（{max(factors) - min(factors):.2f}）")

    spreads = {}
    for key, values in YIELD_CONVERSION_SENSITIVITY.items():
        factors = [
            solve_calibration_factor(
                MODEL_SUGAR, actual,
                replace(YieldConversionSettings(), **{key.lower(): v},
                        calibration_factor=1.0))
            for v in values
        ]
        spreads[key] = max(factors) - min(factors)

    print("\n  前提ごとの幅（較正係数）:")
    for key in sorted(spreads, key=spreads.get, reverse=True):
        print(f"    {key:<26} {spreads[key]:.2f}")

    # 3つとも無視できない幅を持つ。どれが「最大」かは仮定した振り幅次第で
    # 入れ替わるので、順位ではなく「どれも効く」ことを確かめる。
    check(all(spread > 0.3 for spread in spreads.values()),
          "3つの前提いずれも較正係数を 0.3 以上動かす（どれも無視できない）",
          " / ".join(f"{k}={v:.2f}" for k, v in spreads.items()))

    # 較正係数の振れ幅そのもの（全前提を最良・最悪に振った場合）
    best = solve_calibration_factor(MODEL_SUGAR, actual, replace(
        YieldConversionSettings(), sugar_to_dry_matter=0.75, fruit_allocation=0.65,
        fruit_dry_matter_content=0.040, calibration_factor=1.0))
    worst_case = solve_calibration_factor(MODEL_SUGAR, actual, replace(
        YieldConversionSettings(), sugar_to_dry_matter=0.65, fruit_allocation=0.45,
        fruit_dry_matter_content=0.055, calibration_factor=1.0))
    print(f"\n  全前提をモデルに最も有利に振った場合 : 較正係数 {best:.2f} 倍")
    print(f"  全前提をモデルに最も不利に振った場合 : 較正係数 {worst_case:.2f} 倍")
    check(best > 1.5,
          "前提をどれだけ甘くしても、なお 1.5 倍以上の不足が残る"
          "（＝換算の前提だけでは説明できない）",
          f"最良でも {best:.2f} 倍")


if ("中央", 2026) in ACTUAL_YIELD_RECORDS:
    sensitivity_against_actual(ACTUAL_YIELD_RECORDS[("中央", 2026)])
else:
    print("  実収量の記録がないので、この節は飛ばした。")
    print("  （private_data.py は公開リポジトリに入れていない。"
          "手元で動かすと較正係数の感度が出る）")

# =============================================================================
print()
print("=" * 78)
if failures:
    print(f"NG が {len(failures)} 件ある:")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("すべて通過")
print("=" * 78)
