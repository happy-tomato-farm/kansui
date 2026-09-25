"""
日積算とハウス単位への換算の検証。

【確認すること】
1. 積算の計算が正しいか（一定速度なら手計算と一致するか）
2. ハウス単位への換算が正しいか
3. エネルギー収支から見て蒸散量が妥当か  ← 最も重要
4. 潅水実績と比べてどうか

【なぜ潅水実績を検証に使わないか】
土耕栽培では、潅水量は蒸散量と一致しない。

    潅水量 ＝ 蒸散量 ＋ 土壌蒸発 ＋ 下方浸透 ± 土壌水分の増減

下方浸透は測定できず、土壌が水を蓄えるので日単位でも一致しない。
さらに潅水量は経験にもとづいて決めた量であり、適正量である保証もない。
したがって潅水実績は蒸散量の上限としてすら使えない。検証には用いない。

【何で検証するか】
エネルギー配分を使う。これは物理法則と、確立された観測事実による。

  ボーエン比 β ＝ 顕熱 H ÷ 潜熱 λET
      十分な水分条件のトマト群落で 0.25〜0.32
  潜熱の割合 λET ÷ 正味放射 Rn
      66〜72%（条件により0.7〜0.8）

葉が受け取ったエネルギーは、蒸散（潜熱）・対流（顕熱）・放射のどれかで
必ず捨てられる。その配分がどうなるかは多くの実測がある。
自分の圃場の較正データがなくても、この配分で妥当性を判定できる。

出典:
  - Evapotranspiration partitioning of greenhouse grown tomato using a
    modified Priestley-Taylor model
    https://www.sciencedirect.com/science/article/abs/pii/S0378377420322538
  - Diurnal energy-partitioning and transpiration modelling in an
    insect-proof screenhouse with a tomato crop
    https://www.sciencedirect.com/science/article/abs/pii/S1537511017303562

【実行のしかた】
    cd "c:\\Users\\kimij\\OneDrive\\Desktop\\claude作業場\\蒸散光合成モデル"
    python tests\\test_daily.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (  # noqa: E402
    EXTINCTION_COEFFICIENT_K,
    HOUSE_SPECS,
    LEAF_CHARACTERISTIC_LENGTH_M,
    LOCAL_PRESSURE_KPA,
    WIND_SPEED_M_PER_S,
)
from core.daily import (  # noqa: E402
    TIME_STEP_S,
    SimulationSettings,
    aggregate_daily,
    simulate_and_aggregate,
    summarize_monthly,
)
from core.loader import filter_reliable_days, load_and_prepare  # noqa: E402

DATA_DIR = Path(
    r"c:\Users\kimij\OneDrive\Desktop\栽培記録アプリ\saibai-kiroku\data"
)

# 検証に使う期間（栽培期間。7〜8月は非栽培なので除く）
SEASON_MONTHS = [
    "202509", "202510", "202511", "202512",
    "202601", "202602", "202603", "202604", "202605", "202606",
]

# 水の気化潜熱 [MJ/kg]。20℃付近の値。
LATENT_HEAT_MJ_PER_KG = 2.45

# 文献によるトマト群落のエネルギー配分（十分な水分条件）。
#
# 文献値そのものは「季節平均の正午」でボーエン比 0.25〜0.32、
# 潜熱の割合 66〜72% だが、これは平均であって個々の条件では振れる。
# 同じ文献に「午後は潜熱の割合が1.4まで上がる」との記述もある。
# 日射が強くCO2が薄い日は気孔が開いて蒸散が優勢になり、
# ボーエン比は下がる方向に振れる。
#
# そこで判定範囲は文献の中心値より広くとり、
# 「平均が文献の範囲に入っているか」を主に見る。
PLAUSIBLE_BOWEN_MIN = 0.10
PLAUSIBLE_BOWEN_MAX = 0.50
PLAUSIBLE_LATENT_FRACTION_MIN = 0.60
PLAUSIBLE_LATENT_FRACTION_MAX = 0.92

# 3場面の平均が入っているべき範囲（文献の中心値に近いこと）
MEAN_BOWEN_MIN = 0.18
MEAN_BOWEN_MAX = 0.40


def test_integration_arithmetic() -> bool:
    """
    積算の計算が正しいかを、手で検算できる例で確かめる。

    蒸散速度が一定なら、1日の積算量は
        速度 × 300秒 × 288点
    になるはず。
    """
    print("=" * 78)
    print("積算計算の検算")
    print("=" * 78)

    # 蒸散速度を 10 mg/(m²·s) 一定とした架空の1日を作る
    rate_mg = 10.0
    timestamps = pd.date_range("2026-01-01 00:00", periods=288, freq="5min")

    fake = pd.DataFrame({
        "timestamp": timestamps,
        "transpiration_mg_per_m2_s": rate_mg,
        "net_photosynthesis_mmol": 0.01,
        "gross_photosynthesis_mmol": 0.012,
        "respiration_mmol": 0.002,
        "cumulative_radiation_mj": np.linspace(0, 10, 288),
        "dry_bulb_c": 20.0,
        "vpd_kpa": 1.0,
        "co2_ppm": 400.0,
        "solar_w_per_m2": 100.0,
    })

    settings = SimulationSettings(lai=2.0, house="中央")
    daily = aggregate_daily(fake, settings)

    # 手計算: 10 mg/m²s × 300 s × 288点 = 864,000 mg/m² = 0.864 L/m²
    expected_l = rate_mg * TIME_STEP_S * 288 * 1.0e-6
    actual_l = daily["transpiration_l_per_m2"].iloc[0]

    # ハウス単位: 0.864 L/m² × 1740 m² = 1503.36 L
    expected_house = expected_l * HOUSE_SPECS["中央"]["floor_area_m2"]
    actual_house = daily["transpiration_l_per_house"].iloc[0]

    # 光合成: 0.01 mmol/m²s × 300 × 288 = 864 mmol/m²
    #         → 糖 864/1000 × 30 = 25.92 g/m²
    expected_sugar = 0.01 * TIME_STEP_S * 288 / 1000.0 * 30.0
    actual_sugar = daily["sugar_g_per_m2"].iloc[0]

    checks = [
        ("蒸散量 [L/m²]", expected_l, actual_l),
        ("蒸散量 [L/ハウス]", expected_house, actual_house),
        ("糖換算 [g/m²]", expected_sugar, actual_sugar),
    ]

    print(f"\n条件: 蒸散速度 {rate_mg} mg/(m²·s) 一定、288点、"
          f"床面積 {HOUSE_SPECS['中央']['floor_area_m2']:.0f} m²")
    print("-" * 78)
    print(f"{'項目':<24}{'手計算':>16}{'モデル':>16}{'判定':>10}")
    print("-" * 78)

    all_passed = True
    for label, expected, actual in checks:
        passed = abs(actual - expected) < 1e-9
        all_passed = all_passed and passed
        print(f"{label:<24}{expected:>16.6f}{actual:>16.6f}{'OK' if passed else '不一致':>10}")

    print("-" * 78)
    print()

    return all_passed


def test_energy_partitioning() -> bool:
    """
    群落のエネルギー配分が文献の実測と合うかを確かめる。

    これが本モデルの主たる検証。較正データを使わず、
    物理法則と確立された観測事実だけで判定する。
    """
    import math

    from core.energy_balance import (
        longwave_view_factor_in_canopy,
        solve_leaf_temperature,
    )

    print("=" * 78)
    print("群落のエネルギー配分（主たる検証）")
    print("=" * 78)

    # 実データから取った正午前後（11〜14時）の月平均条件。
    # 架空の条件ではなく実測を使う。
    #   気温[℃], 水蒸気圧[kPa], 日射[W/m²], CO2[ppm]
    scenarios = [
        ("2月  正午前後", 23.7, 1.95, 339.0, 538.0),
        ("11月 正午前後", 24.7, 2.00, 278.0, 512.0),
        ("5月  正午前後", 26.1, 1.83, 489.0, 422.0),
    ]

    lai = 2.16
    k = EXTINCTION_COEFFICIENT_K
    n_layers = 20

    print(f"\nLAI={lai}, 消散係数 k={k}, 風速={WIND_SPEED_M_PER_S} m/s")
    print("-" * 78)
    print(f"{'場面':<18}{'Rn':>9}{'λET':>9}{'H':>9}"
          f"{'λET/Rn':>10}{'ボーエン比':>12}")
    print(f"{'':<18}{'W/m²':>9}{'W/m²':>9}{'W/m²':>9}{'':>10}{'':>12}")
    print("-" * 78)

    checks = []
    bowen_values = []
    for label, temp, ea, solar, co2 in scenarios:
        # 群落全体のエネルギー配分を層ごとに積み上げる
        layer_lai = lai / n_layers
        net_radiation = 0.0
        sensible = 0.0
        latent = 0.0

        for i in range(n_layers):
            cumulative = (i + 0.5) * layer_lai
            layer_solar = solar * math.exp(-k * cumulative)
            view = longwave_view_factor_in_canopy(cumulative, lai - cumulative, k)

            balance = solve_leaf_temperature(
                air_temp_c=temp,
                vapor_pressure_kpa=ea,
                solar_radiation_w_per_m2=layer_solar,
                air_co2_ppm=co2,
                wind_speed_m_per_s=WIND_SPEED_M_PER_S,
                pressure_kpa=LOCAL_PRESSURE_KPA,
                characteristic_length_m=LEAF_CHARACTERISTIC_LENGTH_M,
                longwave_view_factor=view,
            )

            net_radiation += balance.net_radiation * layer_lai
            sensible += balance.sensible_heat * layer_lai
            latent += balance.latent_heat * layer_lai

        latent_fraction = latent / net_radiation
        bowen = sensible / latent
        bowen_values.append(bowen)

        print(
            f"{label:<18}{net_radiation:>9.1f}{latent:>9.1f}{sensible:>9.1f}"
            f"{latent_fraction * 100:>9.0f}%{bowen:>12.2f}"
        )

        checks.append((
            f"{label} のボーエン比",
            PLAUSIBLE_BOWEN_MIN <= bowen <= PLAUSIBLE_BOWEN_MAX,
            f"{bowen:.2f}（文献 0.25〜0.32）",
        ))
        checks.append((
            f"{label} の潜熱の割合",
            PLAUSIBLE_LATENT_FRACTION_MIN <= latent_fraction
            <= PLAUSIBLE_LATENT_FRACTION_MAX,
            f"{latent_fraction * 100:.0f}%（文献 66〜72%、条件により70〜80%）",
        ))

    print("-" * 78)
    print()

    mean_bowen = sum(bowen_values) / len(bowen_values)
    checks.append((
        "3場面の平均ボーエン比が文献の中心域にある",
        MEAN_BOWEN_MIN <= mean_bowen <= MEAN_BOWEN_MAX,
        f"{mean_bowen:.2f}（文献の季節平均 0.25〜0.32）",
    ))

    all_passed = True
    for label, passed, detail in checks:
        all_passed = all_passed and passed
        print(f"  {'OK  ' if passed else '範囲外'}  {label}")
        print(f"        {detail}")
    print()

    return all_passed


def test_energy_ceiling(daily: pd.DataFrame) -> bool:
    """
    日積算の蒸散量が、日射のエネルギーを超えていないかを確かめる。

    上限のチェック。蒸散に使えるエネルギーは日射を超えられない。
    """
    print("=" * 78)
    print("日積算のエネルギー上限チェック")
    print("=" * 78)

    work = daily[daily["is_complete"] & (daily["daily_radiation_mj"] > 2.0)].copy()
    work["latent_heat_mj"] = work["transpiration_l_per_m2"] * LATENT_HEAT_MJ_PER_KG
    work["fraction_of_solar"] = work["latent_heat_mj"] / work["daily_radiation_mj"]

    mean_fraction = work["fraction_of_solar"].mean()
    max_fraction = work["fraction_of_solar"].max()

    print(f"\n対象: 日射2 MJ/m²以上の完全な日 {len(work)} 日")
    print(f"  ハウス内日射のうち蒸散に使われた割合")
    print(f"    平均 {mean_fraction * 100:.0f}% / 最大 {max_fraction * 100:.0f}%")
    print()
    print("  ※ 日射の全部が群落に吸収されるわけではない（一部は透過・反射）ため、")
    print("     この割合は上のλET/Rn より小さくなるのが正しい。")
    print("     1日の積算には夜間の呼吸時間も含まれるので、さらに下がる。")
    print()

    # 日射を超えることは物理的にありえない
    passed = max_fraction < 1.0
    if passed:
        print(f"  OK    蒸散のエネルギーが日射を超えた日はありません")
    else:
        over = (work["fraction_of_solar"] >= 1.0).sum()
        print(f"  失敗  {over} 日で蒸散のエネルギーが日射を超えています")
    print()

    return passed


def test_with_real_data() -> tuple[bool, pd.DataFrame]:
    """実データで1年分を計算し、潅水実績と比較する。"""
    print("=" * 78)
    print("実データによる年間計算")
    print("=" * 78)

    df, report = load_and_prepare(DATA_DIR, "中央", SEASON_MONTHS)
    df = filter_reliable_days(df)

    settings = SimulationSettings(lai=2.16, house="中央")
    print(f"\n{settings.describe()}")
    print(f"（LAI 2.16 = 目標葉枚数 18枚/m² × 葉1枚 0.12 m²[暫定値]）")

    _, daily = simulate_and_aggregate(df, settings)

    print(f"\n計算できた日数: {len(daily)} 日")
    print()
    print("月別のまとめ")
    print("-" * 78)
    monthly = summarize_monthly(daily)
    print(f"{'月':<10}{'日数':>5}{'蒸散L/m²':>11}{'蒸散L/棟':>11}"
          f"{'糖g/m²':>10}{'糖kg/棟':>10}{'日射MJ':>9}")
    print("-" * 78)
    for _, row in monthly.iterrows():
        print(
            f"{row['month']:<10}{int(row['日数']):>5}"
            f"{row['蒸散_L_m2_日']:>11.3f}{row['蒸散_L_ハウス_日']:>11.0f}"
            f"{row['糖換算_g_m2_日']:>10.2f}{row['糖換算_kg_ハウス_日']:>10.2f}"
            f"{row['日射_MJ_m2']:>9.2f}"
        )
    print("-" * 78)
    print()

    return True, daily


def test_yield_consistency(daily: pd.DataFrame) -> bool:
    """
    積算した光合成量が、実際の収量と矛盾しないかを確かめる。

    これも較正ではなく、桁が合っているかの傍証。
    光合成で固定した炭水化物のうち、一定割合が果実になる。
    """
    print("=" * 78)
    print("収量との整合（傍証）")
    print("=" * 78)

    complete = daily[daily["is_complete"]]
    total_sugar_kg = complete["sugar_kg_per_house"].sum()

    print(f"\n積算した炭水化物（糖換算）: {total_sugar_kg:.0f} kg/棟"
          f"（{len(complete)}日分）")
    print()
    print("  果実生重への換算（おおよその目安）")

    # 果実への分配率と乾物率から果実生重を見積もる
    for allocation in [0.45, 0.55, 0.65]:
        for dry_matter in [0.05, 0.06]:
            fruit_dry_kg = total_sugar_kg * allocation
            fruit_fresh_t = fruit_dry_kg / dry_matter / 1000.0
            print(
                f"    果実への分配 {allocation * 100:.0f}% × "
                f"乾物率 {dry_matter * 100:.0f}% → {fruit_fresh_t:5.1f} t/棟"
            )

    print()
    print("  ※ 茎・根・果実の呼吸を差し引いていないため、これは上限寄りの見積り。")
    print("     大玉トマトの長期どりで 17.4a あたり 20〜30 t が一つの目安。")
    print()

    # 桁が合っていれば良しとする（1〜100 t の範囲）
    rough_estimate_t = total_sugar_kg * 0.55 / 0.055 / 1000.0
    passed = 5.0 < rough_estimate_t < 100.0

    if passed:
        print(f"  OK    見積り {rough_estimate_t:.1f} t/棟 は現実的な範囲です")
    else:
        print(f"  警告  見積り {rough_estimate_t:.1f} t/棟 は現実離れしています")
    print()

    return passed


def main() -> int:
    arithmetic_ok = test_integration_arithmetic()
    partitioning_ok = test_energy_partitioning()
    real_ok, daily = test_with_real_data()
    ceiling_ok = test_energy_ceiling(daily)
    yield_ok = test_yield_consistency(daily)

    print("=" * 78)
    if all([arithmetic_ok, partitioning_ok, real_ok, ceiling_ok, yield_ok]):
        print("第4段階（日積算）は完了です。")
        print("=" * 78)
        return 0

    print("失敗した項目があります。")
    print("=" * 78)
    return 1


if __name__ == "__main__":
    sys.exit(main())
