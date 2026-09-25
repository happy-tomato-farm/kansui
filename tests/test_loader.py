"""
データ読み込みの検証。

【確認すること】
1. ベクトル版（loader）とスカラー版（psychrometry）の計算結果が一致するか
   … 速度のために2通りの実装があるので、ずれていないことを確かめる
2. 12か月分すべてを読み込めるか、異常値がどれだけあるか
3. 復元した湿度・VPDが生理的に妥当な値か

【実行のしかた】
    cd "c:\\Users\\kimij\\OneDrive\\Desktop\\claude作業場\\蒸散光合成モデル"
    python tests\\test_loader.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (  # noqa: E402
    LOCAL_PRESSURE_KPA,
    MAX_PLAUSIBLE_SOLAR_W_PER_M2,
)
from core.loader import (  # noqa: E402
    _saturation_vapor_pressure_array,
    filter_reliable_days,
    load_and_prepare,
    summarize_daily_quality,
)
from core.psychrometry import (  # noqa: E402
    relative_humidity_from_wet_bulb,
    saturation_vapor_pressure_kpa,
    vapor_pressure_deficit_kpa,
)

DATA_DIR = Path(
    r"c:\Users\kimij\OneDrive\Desktop\栽培記録アプリ\saibai-kiroku\data"
)

ALL_MONTHS = [
    "202507", "202508", "202509", "202510", "202511", "202512",
    "202601", "202602", "202603", "202604", "202605", "202606",
]

# 東ハウスは2025年9月以降しかファイルがない
EAST_MONTHS = [
    "202509", "202510", "202511", "202512",
    "202601", "202602", "202603", "202604", "202605", "202606",
]


def test_vector_matches_scalar() -> bool:
    """
    ベクトル版とスカラー版の計算が一致することを確かめる。

    loader は速度のために numpy でまとめて計算し、
    psychrometry は1点ずつ精密に計算する。実装が2つあるので、
    同じ入力に対して同じ答えを返すことを確認しておく。
    """
    print("=" * 78)
    print("ベクトル版とスカラー版の一致確認")
    print("=" * 78)

    # 農業で現れる範囲の温度を広くとる
    temps = np.arange(-5.0, 45.1, 0.5)

    vector_result = _saturation_vapor_pressure_array(temps)
    scalar_result = np.array([saturation_vapor_pressure_kpa(t) for t in temps])

    max_error = float(np.max(np.abs(vector_result - scalar_result)))
    passed = max_error < 1e-12

    print(f"  飽和水蒸気圧 es(T)  温度範囲 -5〜45℃ を 0.5℃ 刻みで比較")
    print(f"  最大誤差: {max_error:.3e} kPa  → {'OK' if passed else '不一致'}")

    # 飽差についても、乾湿球の組み合わせをいくつか作って比べる
    errors = []
    for dry in [5.0, 15.0, 25.0, 35.0]:
        for depression in [0.0, 1.0, 3.0, 6.0, 10.0]:
            wet = dry - depression
            scalar_vpd = vapor_pressure_deficit_kpa(dry, wet, LOCAL_PRESSURE_KPA)

            # loader と同じ手順でベクトル計算する
            es_dry = _saturation_vapor_pressure_array(np.array([dry]))[0]
            es_wet = _saturation_vapor_pressure_array(np.array([wet]))[0]
            from config import PSYCHROMETRIC_GAMMA_PER_K
            ea = es_wet - PSYCHROMETRIC_GAMMA_PER_K * LOCAL_PRESSURE_KPA * depression
            ea = min(max(ea, 0.0), es_dry)
            vector_vpd = es_dry - ea

            errors.append(abs(scalar_vpd - vector_vpd))

    max_vpd_error = max(errors)
    vpd_passed = max_vpd_error < 1e-12
    print(f"  飽差 VPD  20通りの乾湿球の組み合わせで比較")
    print(f"  最大誤差: {max_vpd_error:.3e} kPa  → {'OK' if vpd_passed else '不一致'}")
    print()

    return passed and vpd_passed


def test_load_all_months() -> bool:
    """12か月分すべてを読み込み、品質を点検する。"""
    print("=" * 78)
    print("実データの読み込みと品質点検")
    print("=" * 78)

    all_ok = True

    for house, months in [("中央", ALL_MONTHS), ("東", EAST_MONTHS)]:
        print(f"\n【{house}ハウス】{months[0]} 〜 {months[-1]}")
        try:
            df, report = load_and_prepare(DATA_DIR, house, months)
        except Exception as e:
            print(f"  読み込みに失敗しました: {type(e).__name__}: {e}")
            all_ok = False
            continue

        print(report.describe())

        # 日射のリセット回数は日数と一致するはず
        expected_days = len(df["timestamp"].dt.date.unique())
        # リセットは「前日から0に戻る」ので、初日を除いた日数になる
        reset_ok = abs(report.radiation_resets - (expected_days - 1)) <= 1
        print(
            f"  リセット回数の妥当性: {report.radiation_resets} 回 vs "
            f"日数 {expected_days} 日 → {'OK' if reset_ok else '要確認'}"
        )

        # クリップ後の瞬時日射が上限に収まっているか
        radiation_ok = report.max_radiation_w_per_m2 <= MAX_PLAUSIBLE_SOLAR_W_PER_M2
        print(
            f"  日射の上限チェック  : 最大 {report.max_radiation_w_per_m2:.0f} W/m² "
            f"→ {'OK' if radiation_ok else '上限を超えている'}"
        )

        # 湿球が使える日が十分に残っているか
        # 栽培期間の主要部分（11月以降）が使えれば実用上は足りる
        reliable_ok = report.reliable_days >= 200
        print(
            f"  使える日数の確保    : {report.reliable_days} 日 "
            f"→ {'OK' if reliable_ok else '不足'}"
        )

        all_ok = all_ok and reset_ok and radiation_ok and reliable_ok

        # 月ごとの内訳を出す
        monthly = (
            df.assign(month=df["timestamp"].dt.strftime("%Y-%m"))
            .groupby("month")
            .agg(
                異常率=("was_clipped", "mean"),
                使える日=("wet_bulb_reliable", lambda s: s.iloc[::288].sum()),
            )
        )
        print(f"  月別の湿球異常率:")
        for month, row in monthly.iterrows():
            mark = "×" if row["異常率"] > 0.05 else "○"
            print(
                f"    {mark} {month}  異常率 {row['異常率'] * 100:5.2f}%"
            )

    print()
    return all_ok


def test_restored_humidity_is_plausible() -> bool:
    """
    復元した湿度・VPDが生理的に妥当な範囲かを確かめる。

    数値計算として正しくても、値が現実離れしていれば
    どこかで前提を取り違えている。日変化のパターンで確認する。
    """
    print("=" * 78)
    print("復元した湿度・飽差の妥当性")
    print("=" * 78)

    df, _ = load_and_prepare(DATA_DIR, "中央", ["202511"])
    df["hour"] = df["timestamp"].dt.hour

    hourly = df.groupby("hour").agg(
        気温=("dry_bulb_c", "mean"),
        湿球=("wet_bulb_used_c", "mean"),
        相対湿度=("relative_humidity", "mean"),
        飽差kPa=("vpd_kpa", "mean"),
        飽差g_m3=("vpd_g_per_m3", "mean"),
        日射=("solar_w_per_m2", "mean"),
        CO2=("co2_ppm", "mean"),
    )

    print("\n2025年11月・中央ハウスの平均日変化")
    print("-" * 78)
    print(f"{'時':>3} {'気温℃':>7} {'湿球℃':>7} {'RH%':>6} "
          f"{'VPD kPa':>9} {'VPD g/m³':>9} {'日射W/m²':>9} {'CO2ppm':>8}")
    print("-" * 78)
    for hour in range(0, 24, 2):
        row = hourly.loc[hour]
        print(
            f"{hour:>3} {row['気温']:>7.1f} {row['湿球']:>7.1f} "
            f"{row['相対湿度'] * 100:>6.1f} {row['飽差kPa']:>9.3f} "
            f"{row['飽差g_m3']:>9.2f} {row['日射']:>9.0f} {row['CO2']:>8.0f}"
        )
    print("-" * 78)

    checks = []

    # 相対湿度が 0〜100% に収まっているか
    rh_min = df["relative_humidity"].min()
    rh_max = df["relative_humidity"].max()
    checks.append((
        "相対湿度が 0〜100% の範囲に収まる",
        0.0 <= rh_min and rh_max <= 1.0,
        f"最小 {rh_min * 100:.1f}% / 最大 {rh_max * 100:.1f}%",
    ))

    # 夜間のほうが日中より湿度が高いか（蒸散と気温上昇の当然の帰結）
    night_rh = df[df["hour"].isin([2, 3, 4])]["relative_humidity"].mean()
    day_rh = df[df["hour"].isin([12, 13, 14])]["relative_humidity"].mean()
    checks.append((
        "夜間の湿度 > 日中の湿度",
        night_rh > day_rh,
        f"夜間 {night_rh * 100:.1f}% > 日中 {day_rh * 100:.1f}%",
    ))

    # 飽差が負になっていないか
    checks.append((
        "飽差が負にならない",
        df["vpd_kpa"].min() >= 0.0,
        f"最小 {df['vpd_kpa'].min():.4f} kPa",
    ))

    # 日中の飽差がトマトの適正域に近いか
    # （0.2〜3.0 kPa を「ありうる範囲」とする。適正は 0.5〜1.5）
    day_vpd = df[df["hour"].isin([11, 12, 13, 14])]["vpd_kpa"].mean()
    checks.append((
        "日中の平均飽差がありうる範囲（0.2〜3.0 kPa）",
        0.2 <= day_vpd <= 3.0,
        f"{day_vpd:.3f} kPa（トマトの適正域は 0.5〜1.5 kPa）",
    ))

    # 夜間は日射がゼロか
    night_solar = df[df["hour"].isin([0, 1, 2, 3, 22, 23])]["solar_w_per_m2"].max()
    checks.append((
        "夜間の日射がゼロ",
        night_solar == 0.0,
        f"夜間の最大 {night_solar:.1f} W/m²",
    ))

    # 日中にCO2が下がるか（光合成の証拠）
    night_co2 = df[df["hour"].isin([2, 3, 4])]["co2_ppm"].mean()
    day_co2 = df[df["hour"].isin([11, 12, 13])]["co2_ppm"].mean()
    checks.append((
        "日中のCO2 < 夜間のCO2（光合成による消費）",
        day_co2 < night_co2,
        f"夜間 {night_co2:.0f} ppm → 日中 {day_co2:.0f} ppm "
        f"（{night_co2 - day_co2:.0f} ppm の低下）",
    ))

    print()
    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '失敗'}  {label}")
        print(f"        {detail}")

    print()
    return all_passed


def main() -> int:
    results = [
        test_vector_matches_scalar(),
        test_load_all_months(),
        test_restored_humidity_is_plausible(),
    ]

    print("=" * 78)
    if all(results):
        print("第2段階（データ読み込み）は完了です。")
        print("=" * 78)
        return 0

    print("失敗した項目があります。先に進む前に修正が必要です。")
    print("=" * 78)
    return 1


if __name__ == "__main__":
    sys.exit(main())
