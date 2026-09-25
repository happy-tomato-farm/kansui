"""
「今日の日射予測」から、その日の潅水量のめやすを出す。

【このモジュールの立ち位置】
これまでの core/ は「実測データを流し込んで、起きたことを再現する」ための
ものだった。ここは逆向きで、**これから起きることを見積もる**ために使う。

    今日の日射予測 [MJ/m²]
        ↓（1日の日射の形に展開して、5分ごとの計算にかける）
    予測蒸散量 [L/m²/日]
        ↓
    10MJあたり潅水量 [L/m²]  ← 液肥混入機レシピアプリに入れる数字
        ↓（土壌の水収支）
    流亡の見込み・土のpF・空気率

【なぜ「10MJあたり潅水量」を出すのか】
既に使っている液肥混入機レシピのアプリが、この係数を手入力で受け取って
潅水時間を決めている。

    潅水量 ＝ 実質日射 ÷ 10 × 「10MJあたり潅水量」

デルフィのコンサルが 1.5〜3.0 L/m² と言ったのはこの数字。
これを勘ではなく物理から出すのが、このモジュールの仕事。

【日射の読みかえ】
ハウスのセンサーは屋根の下にあるので、外より小さく出る。
既存アプリと同じ約束にそろえてある。

    外の日射   ＝ センサー値 ÷ 透過率（既定 0.65）
    実質日射   ＝ 外の日射 × フィルム劣化係数（既定 1.00）

蒸散の計算に使うのは**センサー値そのもの**（作物が実際に浴びる光）で、
「10MJあたり」の分母に使うのは**実質日射**。既存アプリの約束に合わせる。

どちらの係数も advise() の引数で変えられる。画面からも動かせる。

【1日の日射の形】
日積算しか分からないので、正弦の半周期で日の出から日没まで配分する。
実測の5分値と比べると、晴天日はよく合い、曇天日は山が低くなだらかになる。
日積算が同じなら蒸散量の差は数%にとどまる（tests/test_advisor.py で確認）。
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass, field

from config import (
    MOLAR_MASS_CO2_G_PER_MOL,
    DAYTIME_WINDOW_HOUR,
    DIURNAL_SMOOTHING_CORRECTION,
    DIURNAL_TEMP_AMPLITUDE_C,
    DIURNAL_VAPOR_PRESSURE_AMPLITUDE_KPA,
    EXTINCTION_COEFFICIENT_K,
    FILM_DEGRADATION_FACTOR,
    HOUSE_SPECS,
    LEAF_CHARACTERISTIC_LENGTH_M,
    LOCAL_PRESSURE_KPA,
    RECIPE_RADIATION_COEF,
    WIND_SPEED_M_PER_S,
)
from core.canopy import calculate_canopy_gas_exchange
from core.psychrometry import saturation_vapor_pressure_kpa
from core.soil import (
    FIELD_CAPACITY_WATER_CONTENT,
    air_filled_porosity,
    describe_soil_state,
)
from core.solar import solar_day
from core.water_balance import (
    DailyWaterInput,
    WaterBalanceSettings,
    step_one_day,
)

SECONDS_PER_STEP = 300.0          # 5分
GLUCOSE_G_PER_MOL = 180.0         # グルコースの分子量
STEPS_PER_DAY = 288
MJ_PER_J = 1.0e-6
GRAMS_PER_LITER = 1000.0

#: 日射がこの値を下回る時間帯は夜とみなす [W/m²]
NIGHT_THRESHOLD_W_PER_M2 = 5.0


# =============================================================================
# 1. 1日ぶんの気象を組み立てる
# =============================================================================

@dataclass(frozen=True)
class DayForecast:
    """その日の見込み。分かる範囲だけ入れればよい。

    日射以外を省くと、その月の実績平均（MONTHLY_NORMS）を使う。
    """

    #: ハウスセンサー基準の日積算日射 [MJ/m²]。既存アプリの「1日の日射予測」。
    sensor_radiation_mj: float

    #: 日中の平均気温 [℃]
    mean_temp_c: float

    #: 日中の平均相対湿度 [0〜1]
    mean_relative_humidity: float

    #: 日中の平均CO2濃度 [ppm]
    mean_co2_ppm: float

    #: 日射を配る時間帯 [時]。
    #:
    #: ★ forecast_from_date / forecast_from_month で作れば、
    #:   その日の**実際の日の出・日没**が入る。ふだんはそちらを使うこと。
    #:   ここの既定値 (6, 17) は、日付が分からないまま手で組んだときの控え。
    #:
    #: 【なぜ日の出〜日没でないといけないか】（2026-09-24）
    #: 1日ぶんの日射をこの時間帯に配るので、窓が狭いと真昼が強くなりすぎる。
    #: 5月の実際の日長は14.2時間なのに11時間に詰め込むと、
    #: 明るい日の蒸散が12%も過大になっていた。
    #:
    #:   窓          比の中央値  ±20%以内  平均絶対誤差  16〜30MJの日
    #:   6〜17時        1.008      93%     0.159 L/m²      1.125
    #:   日の出〜日没    1.021      95%     0.071 L/m²      1.028  ★
    #:
    #: ★窓を変えたら、気温・湿度の推定式も同じ窓で当てはめ直すこと。
    #:   平均の定義が変わるので、片方だけ直すとかえってずれる。
    daylight_hours: tuple[float, float] = DAYTIME_WINDOW_HOUR

    #: 気温の日較差 [℃]。日射の形に合わせて気温も揺らす。
    #: 実測の値と、日射に比例させてはいけない理由は config.py 第9-1節。
    temp_amplitude_c: float = DIURNAL_TEMP_AMPLITUDE_C

    #: 日中の平均水蒸気圧 [kPa]。
    #:
    #: ★蒸散の計算はこちらを使う。mean_relative_humidity は表示用。
    #:   相対湿度は気温の従属変数（RH = ea / es(気温)）なので、気温を
    #:   揺らすと水蒸気が出入りしなくても動いてしまう。状態量として
    #:   持つべきは水蒸気圧のほう。理由は config.py 第9-2節。
    #:
    #: None のときは平均気温と相対湿度から es(気温) × RH で作る。
    #: 手で DayForecast を組んだときの受け皿。
    mean_vapor_pressure_kpa: float | None = None

    #: 夜（日射がゼロの時刻）の気温 [℃]。
    #:
    #: None なら日中の平均気温をそのまま使う。以前は日中の振幅から
    #: 外挿していたが、暖房の効いたハウスの夜と合わず、夜間蒸散が
    #: 3倍ほど過大になっていた（config.py 第9-4節）。
    night_temp_c: float | None = None

    #: 夜の水蒸気圧 [kPa]。None なら日中の値を使う。
    night_vapor_pressure_kpa: float | None = None

    #: 夜のCO2濃度 [ppm]。None なら日中の値を使う。
    #: 夜は呼吸で溜まるので日中より 100〜200 ppm 高い。
    night_co2_ppm: float | None = None


#: 月別の実績平均。いまは参考表で、計算には使っていない
#: （推定は下の HARMONIC_COEFFICIENTS で行う）。値の桁を確かめるのに使う。
#:
#: 出典: 中央ハウス 2025-09〜2026-06 の5分値（湿球が信頼できる日のみ）。
#:
#: ★★★ これは「日の出〜日没の平均」であって24時間平均ではない ★★★
#: 蒸散が起きるのは日中だけなので、日中の値でないと合わない。
#: 24時間平均を使うと気温が3〜4℃低く、湿度が5〜7ポイント高く出る。
#: どちらも VPD を小さくする向きなので、蒸散量が3割ほど過小になる
#: （実際にそうなって直した）。
#:
#: 2026-09-24 に窓を 6〜17時 から 日の出〜日没 に直した。理由は
#: HARMONIC_COEFFICIENTS の上のコメントに書いてある。
#: 9月の値が 30.3℃ → 23.2℃ と大きく変わったのは、窓のせいではなく
#: 「健全な日が1日しかない月を調和関数でならした」ため。こちらが妥当。
MONTHLY_NORMS: dict[int, dict[str, float]] = {
    1: {"temp_c": 20.1, "rh": 0.781, "co2_ppm": 672.0},
    2: {"temp_c": 20.1, "rh": 0.766, "co2_ppm": 627.0},
    3: {"temp_c": 20.1, "rh": 0.710, "co2_ppm": 558.0},
    4: {"temp_c": 20.6, "rh": 0.648, "co2_ppm": 484.0},
    5: {"temp_c": 22.0, "rh": 0.642, "co2_ppm": 441.0},
    6: {"temp_c": 23.8, "rh": 0.705, "co2_ppm": 425.0},
    # 7・8月は栽培していないため実測がない。調和関数の外挿値。
    7: {"temp_c": 25.0, "rh": 0.787, "co2_ppm": 425.0},
    8: {"temp_c": 24.8, "rh": 0.829, "co2_ppm": 435.0},
    9: {"temp_c": 23.2, "rh": 0.805, "co2_ppm": 465.0},
    10: {"temp_c": 21.4, "rh": 0.756, "co2_ppm": 521.0},
    11: {"temp_c": 20.2, "rh": 0.733, "co2_ppm": 598.0},
    12: {"temp_c": 20.0, "rh": 0.753, "co2_ppm": 658.0},
}

#: 実測が薄い月。めやすを出すときに一言添えるために使う。
MONTHS_WITH_WEAK_DATA = {7, 8, 9}

#: 月ごとの日射の平均 [MJ/m²/日]（ハウスセンサー基準）。
#: 日付ベースの推定（下の調和関数）が使えないときの控えとして残してある。
MONTHLY_MEAN_RADIATION_MJ: dict[int, float] = {
    1: 7.46, 2: 8.27, 3: 10.48, 4: 11.96, 5: 14.60, 6: 11.78,
    7: 12.0, 8: 12.0, 9: 11.0, 10: 8.5, 11: 6.49, 12: 5.83,
}


# =============================================================================
# 1-2. 日付ベースの気象推定（調和関数）
# =============================================================================
#
# 【なぜ月区切りをやめたか】
# 月でまとめると、変化の大きい月で上旬と下旬の差が埋もれる。
# 旬（10日区切り）で実測を見ると、月内で気温が ±1.5℃、
# 湿度が ±0.09 も動く月があった（4〜6月）。
#
# 【なぜ旬の表をそのまま使わないか】
# 旬ごとの平均には天気のノイズが乗っている。たとえば4月中旬は
# 平均日射 15.38 MJ（平年比 +3.4）と明るい日に偏っており、
# 湿度 0.536 という値はその天気の反映であって季節の姿ではない。
# これを季節の値として固定すると、下の日射補正と二重に効いてしまう。
#
# 【やったこと】
# 季節成分と天気成分を同時に推定した。
#
#     値 = 季節の滑らかな曲線(日付) + 傾き × (その日の日射 − 平年の日射(日付))
#
# 季節曲線は1年で1周する調和関数（sin/cos を2次まで）。
# 10日区切りの標本ノイズを自然にならしてくれる。
#
# 【当てはまり（実測257日）】
#            季節だけ   季節+日射
#   気温      0.448  →  0.805
#   湿度      0.375  →  0.848
#   CO2       0.559  →  0.655
#
# 日射の傾きは、単純に月内偏差から求めた値とほぼ一致した
# （気温 +0.396 vs +0.382、湿度 -0.0171 vs -0.0168）。
#
# 【当てはめた範囲】
# 2025年9月〜2026年6月の実測。7・8月は栽培していないため、
# その時期の値は曲線の外挿になる。

#: ★★★ 平均を取る時間帯は「日の出〜日没」★★★（2026-09-24 に直した）
#:
#: 以前は 6〜17時 の平均に当てはめていた。ところが advise() はこの値を
#: 日の出〜日没に展開して使う。5月なら日の出4.8時〜日没19.0時で14.2時間あり、
#: 11時間ぶんの平均を14.2時間ぶんとして使っていたことになる。
#: 朝夕の湿って涼しい時間帯が平均に入らないので、全時間が乾きすぎていた。
#:
#: 窓をそろえると、実測257日に対する蒸散量の成績がこうなった。
#:
#:                          比の中央値  ±20%以内  平均絶対誤差
#:   6〜17時の平均（旧）        1.011      94%     0.147 L/m²
#:   日の出〜日没の平均（新）    1.021      95%     0.071 L/m²  ★
#:
#: 平均絶対誤差が半分になる。日射の強さによる偏りも平らになった。
#:
#:   日射[MJ]    0〜4   4〜8  8〜12 12〜16 16〜30   幅
#:   旧         0.936  0.966  1.007  1.087  1.132  0.196
#:   新         0.938  1.025  1.041  1.013  1.028  0.103
#:
#: 【入力側の R² は下がるが、それでよい】
#: 当てはまりは 気温 0.805→0.719、湿度 0.848→0.808 と下がる。
#: 窓が広がって朝夕のばらつきが入るぶん、日ごとの平均が季節の曲線から
#: 外れやすくなるため。だが目的は蒸散量を当てることであって、
#: 入力値そのものの当てはまりではない。下流の成績で判断すること。

#: 調和関数の係数。並びは [定数, sin(θ), cos(θ), sin(2θ), cos(2θ)]。
#: θ = 2π × 通日 / 365.25
#:
#: 【night_ 付きの項目について】（2026-09-24 に追加）
#: 日の出〜日没を「日中」、それ以外を「夜」として別々に当てはめた。
#: 夜の値は日射とほとんど関係しない（当てはめた傾きは気温 +0.007 ℃/MJ、
#: ea −0.005 kPa/MJ）ので、下の RADIATION_SLOPES には入れていない。
#:
#: 【ea_kpa を足して rh を計算から外した理由】
#: 相対湿度は気温の従属変数なので、気温を揺らすと二重に効く。
#: 状態量として持つべきは水蒸気圧。config.py 第9-2節に実測の裏づけ。
#: rh は画面表示用に残してある。
#:
#: 【当てはまり（実測257日・残差RMS）】
#:   日中の気温  1.23 ℃      夜の気温  1.03 ℃
#:   日中の ea   0.188 kPa   夜の ea   0.143 kPa
#:   日中のCO2   82.1 ppm    夜のCO2   74.9 ppm
HARMONIC_COEFFICIENTS: dict[str, list[float]] = {
    # 日射だけは24時間の積算なので、窓を変えても値は同じ。
    "radiation_mj": [9.00289, 2.70096, -2.42777, -0.82946, -0.00510],
    "temp_c": [21.77877, -0.97899, -2.29576, 0.63724, 0.52860],
    "ea_kpa": [1.94372, -0.27169, -0.24689, 0.23352, 0.11321],
    # ★CO2 だけは日射に対して直線ではない。下の CO2_RADIATION_RESPONSE と
    #   セットで使う（この係数は指数項と同時に解いたもの）。
    "co2_ppm": [464.99124, 33.75404, 111.83040, 6.10927, 25.79224],
    # 表示用。蒸散の計算には ea_kpa を使う。
    "rh": [0.74320, -0.05341, 0.00896, 0.05446, 0.01835],
    # 夜（日没〜日の出）。暖房の設定温度がそのまま出ている。
    "night_temp_c": [16.84740, -2.24999, -3.81500, 2.03915, 0.77151],
    "night_ea_kpa": [1.69018, -0.26598, -0.39187, 0.27362, 0.07805],
    "night_co2_ppm": [610.77696, 25.45644, 174.47879, 22.68490, 18.08976],
}

#: 平年の日射からのずれ1MJあたり、気温・湿度を動かす量。
#: 季節成分と同時に推定した値。
#: 夜の項目は日射でほとんど動かないので入れていない（上のコメント参照）。
#: CO2 は直線では表せないので、下の CO2_RADIATION_RESPONSE で別扱い。
RADIATION_SLOPES: dict[str, float] = {
    "temp_c": 0.36107,
    "ea_kpa": -0.00187,
    "rh": -0.01606,
}

#: CO2 が日射に対してどう動くか。★直線ではなく指数で減衰する。
#:
#:     CO2 = 季節の水準 + AMPLITUDE_PPM × exp( −日射[MJ] / SCALE_MJ )
#:
#: 【なぜ直線ではいけないか】
#: 暗い日は締め切るうえに光合成でも消費されないので急に溜まる。
#: 明るくなると換気と光合成で下がるが、外気の水準で下げ止まる。
#: つまり「暗い側に張り出し、明るい側で頭打ち」という飽和の形になる。
#: 直線で当てはめると、暗い日に上げ足りず（−105 ppm）、
#: 明るい日に下げすぎる（20MJ で −58 ppm）。
#:
#: 【当てはまり（実測257日）】
#:   形                残差RMS   R²     4MJ未満のずれ
#:   直線（前）          82.1   0.593    −40 ppm
#:   exp(−日射/4)       72.2   0.685     +5 ppm  ★
#:   1/(1+日射/4)       73.4   0.675     +1 ppm
#:   直線 + 2乗         72.5   0.683    +14 ppm
#:
#: 【効き方】
#: CO2 が低いとモデルは気孔を開けるので蒸散が増える。暗い日の CO2 を
#: 115 ppm も低く見積もっていたのが、暗い日（0〜2MJ）の蒸散が
#: 1.67倍に膨らんでいた主因だった。
#:
#: 夜の CO2（night_co2_ppm）にも同じ形が見られるが、そちらは日射項を
#: 入れていない。改善がわずか（残差RMS 74.9 → 73.7）で、夜は気孔が
#: 閉じていて蒸散にほとんど効かないため。
CO2_RADIATION_RESPONSE = {
    "AMPLITUDE_PPM": 324.098,   # 日射ゼロのときの上乗せ [ppm]
    "SCALE_MJ": 4.0,            # この日射量ごとに 1/e に減る [MJ/m²]
}


def _harmonic_value(day_of_year: int, coefficients: list[float]) -> float:
    """調和関数で、その日の季節の値を求める。"""
    angle = 2.0 * math.pi * day_of_year / 365.25
    value = coefficients[0]
    for order in range(1, (len(coefficients) - 1) // 2 + 1):
        value += coefficients[2 * order - 1] * math.sin(order * angle)
        value += coefficients[2 * order] * math.cos(order * angle)
    return value


def normal_radiation_mj(day_of_year: int) -> float:
    """その日の平年の日射 [MJ/m²]（ハウスセンサー基準）。"""
    return max(_harmonic_value(day_of_year, HARMONIC_COEFFICIENTS["radiation_mj"]), 0.0)

#: その日の日射が平年からずれた分だけ、気温・湿度・CO2 を動かす傾き。
#:
#: 【なぜ要るか】
#: 月平均だけで気温・湿度を決めると、曇天日の蒸散を大きく見積もりすぎる。
#: 実際には暗い日はハウスを閉めるので湿度が上がり、気温も上がらない。
#: 日射4MJ未満の日で、予測が実測の2.8倍にもなっていた。
#:
#: 【裏づけ】
#: 中央ハウス242日の日中（6〜17時）平均を、日射で層別すると
#:
#:   日射帯     日中RH  日中気温  日中CO2
#:   0〜4 MJ     0.845    17.7℃    722 ppm
#:   16〜30 MJ   0.539    24.9℃    430 ppm
#:
#: 相関は RH −0.89、気温 +0.80、CO2 −0.65。
#: 月内の偏差だけで見ても −0.88 / +0.79 / −0.47 と残るので、
#: 季節の効果ではなく日射そのものの効果。
RADIATION_RESPONSE = {
    "RH_PER_MJ": -0.0168,       # 日射1MJ明るいごとに相対湿度が下がる量
    "TEMP_C_PER_MJ": 0.382,     # 日射1MJ明るいごとに日中気温が上がる量 [℃]
    "CO2_PPM_PER_MJ": -10.62,   # 日射1MJ明るいごとにCO2が下がる量 [ppm]
}

#: 上の傾きで外挿しすぎないための範囲。
#: 水蒸気圧の下限 0.3 kPa は 0℃ の飽和水蒸気圧の約半分、
#: 上限 4.0 kPa は 29℃ の飽和水蒸気圧にあたる。ハウス内でこの外へ出ない。
NORM_BOUNDS = {
    "RH": (0.35, 0.95),
    "TEMP_C": (10.0, 35.0),
    "CO2_PPM": (350.0, 1200.0),
    "EA_KPA": (0.3, 4.0),
    "NIGHT_TEMP_C": (5.0, 30.0),
    "NIGHT_EA_KPA": (0.3, 4.0),
    "NIGHT_CO2_PPM": (350.0, 1500.0),
}


def forecast_from_day_of_year(
    sensor_radiation_mj: float,
    day_of_year: int,
    temp_c: float | None = None,
    relative_humidity: float | None = None,
    co2_ppm: float | None = None,
) -> DayForecast:
    """日射と日付から見込みを組み立てる。省いた項目は実測から推定する。

    推定は「季節の曲線（日付）＋ 日射が平年からずれた分」で行う。
    手で値を渡せば、そちらが優先される。

    Args:
        sensor_radiation_mj: ハウスセンサー基準の日積算日射 [MJ/m²]
        day_of_year: 通日（1〜366）
        temp_c: 日中の平均気温 [℃]。省くと推定する。
        relative_humidity: 日中の平均相対湿度 [0〜1]。省くと推定する。
        co2_ppm: 日中の平均CO2 [ppm]。省くと推定する。
    """
    if not 1 <= day_of_year <= 366:
        raise ValueError(
            f"通日は 1〜366 で指定することだ。渡された値: {day_of_year}"
        )

    # その日が平年よりどれだけ明るいか。
    # 明るい日は乾いて暑く CO2 が下がる、という実測の傾きで補正する。
    brightness = sensor_radiation_mj - normal_radiation_mj(day_of_year)

    def estimate(key: str, bound_key: str) -> float:
        """季節の曲線に、その日の明るさによる補正を足す。"""
        seasonal = _harmonic_value(day_of_year, HARMONIC_COEFFICIENTS[key])
        value = seasonal + RADIATION_SLOPES[key] * brightness
        low, high = NORM_BOUNDS[bound_key]
        return min(max(value, low), high)

    def estimate_night(key: str, bound_key: str) -> float:
        """夜の値。日射の項を持たない（夜は日射と無関係なため）。"""
        value = _harmonic_value(day_of_year, HARMONIC_COEFFICIENTS[key])
        low, high = NORM_BOUNDS[bound_key]
        return min(max(value, low), high)

    def estimate_co2() -> float:
        """日中のCO2。日射に対して指数で減衰する（上のコメント参照）。"""
        seasonal = _harmonic_value(day_of_year, HARMONIC_COEFFICIENTS["co2_ppm"])
        buildup = CO2_RADIATION_RESPONSE["AMPLITUDE_PPM"] * math.exp(
            -sensor_radiation_mj / CO2_RADIATION_RESPONSE["SCALE_MJ"]
        )
        low, high = NORM_BOUNDS["CO2_PPM"]
        return min(max(seasonal + buildup, low), high)

    # 日射を配る時間帯は、その日の実際の日の出〜日没にする。
    # 推定した気温・湿度も同じ時間帯の平均として当てはめてあるので、
    # ここをそろえないと平均の意味がずれる（DayForecast のコメント参照）。
    day = solar_day(day_of_year)

    # 気温と相対湿度を手で渡された場合は、水蒸気圧もその2つから作る。
    # そうしないと「渡した湿度」と「使われる湿り具合」が食い違う。
    given_air = temp_c is not None and relative_humidity is not None
    mean_temp = estimate("temp_c", "TEMP_C") if temp_c is None else temp_c
    mean_rh = (
        estimate("rh", "RH") if relative_humidity is None else relative_humidity
    )
    vapor_pressure = (
        saturation_vapor_pressure_kpa(mean_temp) * mean_rh
        if given_air
        else estimate("ea_kpa", "EA_KPA")
    )

    return DayForecast(
        sensor_radiation_mj=sensor_radiation_mj,
        mean_temp_c=mean_temp,
        mean_relative_humidity=mean_rh,
        mean_co2_ppm=estimate_co2() if co2_ppm is None else co2_ppm,
        daylight_hours=(day.sunrise_h, day.sunset_h),
        mean_vapor_pressure_kpa=vapor_pressure,
        night_temp_c=estimate_night("night_temp_c", "NIGHT_TEMP_C"),
        night_vapor_pressure_kpa=estimate_night("night_ea_kpa", "NIGHT_EA_KPA"),
        night_co2_ppm=estimate_night("night_co2_ppm", "NIGHT_CO2_PPM"),
    )


def forecast_from_date(
    sensor_radiation_mj: float,
    date: dt.date,
    temp_c: float | None = None,
    relative_humidity: float | None = None,
    co2_ppm: float | None = None,
) -> DayForecast:
    """日付（datetime.date）から見込みを組み立てる。ふだんはこれを使う。"""
    return forecast_from_day_of_year(
        sensor_radiation_mj, date.timetuple().tm_yday,
        temp_c=temp_c, relative_humidity=relative_humidity, co2_ppm=co2_ppm,
    )


def forecast_from_month(
    sensor_radiation_mj: float,
    month: int,
    temp_c: float | None = None,
    relative_humidity: float | None = None,
    co2_ppm: float | None = None,
) -> DayForecast:
    """月だけから見込みを組み立てる（その月の15日として扱う）。

    日付が分かるなら forecast_from_date を使ったほうが当たる。
    月内でも上旬と下旬で気温が1.5℃ほど違う時期がある。
    """
    if not 1 <= month <= 12:
        raise ValueError(f"月は 1〜12 で指定することだ。渡された値: {month}")
    return forecast_from_day_of_year(
        sensor_radiation_mj, dt.date(2026, month, 15).timetuple().tm_yday,
        temp_c=temp_c, relative_humidity=relative_humidity, co2_ppm=co2_ppm,
    )


def _diurnal_shape(forecast: DayForecast) -> list[float]:
    """1日を5分刻みにしたときの、日射の配分（合計1）。

    日の出から日没まで正弦の半周期。夜はゼロ。
    """
    sunrise, sunset = forecast.daylight_hours
    if not 0.0 <= sunrise < sunset <= 24.0:
        raise ValueError(
            f"日の出・日没は 0 ≦ 日の出 < 日没 ≦ 24 で指定する。"
            f"渡された値: {forecast.daylight_hours}"
        )

    hours_per_step = 24.0 / STEPS_PER_DAY
    weights = []
    for index in range(STEPS_PER_DAY):
        hour = (index + 0.5) * hours_per_step
        if sunrise <= hour < sunset:
            phase = (hour - sunrise) / (sunset - sunrise)
            weights.append(math.sin(math.pi * phase))
        else:
            weights.append(0.0)

    total = sum(weights)
    if total <= 0.0:
        raise ValueError("日中の時間帯に1ステップも入らなかった。")
    return [w / total for w in weights]


# =============================================================================
# 2. 計算する
# =============================================================================

@dataclass(frozen=True)
class IrrigationAdvice:
    """その日の潅水のめやす。"""

    # --- 入力の読みかえ ---
    sensor_radiation_mj: float       # ハウスセンサー基準 [MJ/m²]
    outside_radiation_mj: float      # 外の日射（近似）[MJ/m²]
    effective_radiation_mj: float    # 実質日射（既存アプリが使う値）[MJ/m²]

    # --- 蒸散 ---
    transpiration_l_per_m2: float    # 予測蒸散量 [L/m²/日]
    transpiration_l_per_house: float # ハウス1棟あたり [L/日]
    peak_transpiration_mm_per_h: float   # 正午前後の蒸散速度 [mm/h]

    # --- 光合成（★蒸散より信頼度が低い。較正係数を見ること）---
    co2_fixed_g_per_m2: float        # CO2固定量 [g-CO2/m²/日]
    sugar_g_per_m2: float            # 糖換算 [g-糖/m²/日]
    sugar_kg_per_house: float        # ハウス1棟あたり [kg-糖/日]

    # --- 潅水のめやす ---
    water_per_10mj_l_per_m2: float   # ★液肥混入機レシピに入れる数字
    recommended_irrigation_l_per_m2: float
    recommended_irrigation_l_per_house: float

    # --- 土壌 ---
    expected_drainage_l_per_m2: float
    drainage_fraction: float
    end_pf: float
    min_air_filled_porosity: float
    hours_above_field_capacity: float

    # --- 添えもの ---
    warnings: list[str] = field(default_factory=list)

    @property
    def is_over_irrigating(self) -> bool:
        """空気率が根に苦しい水準まで下がる見込みか。"""
        return self.min_air_filled_porosity < 0.08


def advise(
    forecast: DayForecast,
    house: str = "中央",
    lai: float = 2.16,
    start_water_content: float | None = None,
    leaching_fraction: float = 0.20,
    wind_speed_m_per_s: float = WIND_SPEED_M_PER_S,
    water_balance_settings: WaterBalanceSettings | None = None,
    recipe_radiation_coef: float = RECIPE_RADIATION_COEF,
    film_degradation_factor: float = FILM_DEGRADATION_FACTOR,
) -> IrrigationAdvice:
    """その日の潅水量のめやすを計算する。

    Args:
        forecast: その日の見込み
        house: "中央" または "東"
        lai: 葉面積指数
        start_water_content: 朝の土壌体積含水率。省くと圃場容水量から始める。
        leaching_fraction: 塩を流すために上乗せする割合。
            0.20 なら蒸散量の1.2倍を潅水する。土耕で EC を抑えるには
            ある程度の余剰が要るため、0 にはしない。
        wind_speed_m_per_s: 群落内の風速
        water_balance_settings: 土壌側の設定
        recipe_radiation_coef: 液肥混入機レシピと合わせるための換算係数。
            ★物理量ではない。「実質日射」＝ センサー値 ÷ この係数 で、
            レシピ側の coef（6.5 = 0.65）と同じ値にしておけば、
            分子と分母で打ち消し合って狙った量の水が出る。
            ハウス被覆の実測透過率（季節で 0.62〜0.72）とは別物。
        film_degradation_factor: フィルム劣化による減衰係数。
            これもレシピ側の film と揃えるための値。

    Returns:
        IrrigationAdvice
    """
    if house not in HOUSE_SPECS:
        raise ValueError(
            f"ハウス名は {list(HOUSE_SPECS.keys())} のどれか。"
            f"渡された値: '{house}'"
        )
    if forecast.sensor_radiation_mj < 0.0:
        raise ValueError(
            f"日射予測が負: {forecast.sensor_radiation_mj} MJ/m²"
        )
    if leaching_fraction < 0.0:
        raise ValueError(
            f"塩を流すための上乗せ割合が負: {leaching_fraction}"
        )
    if not 0.0 < recipe_radiation_coef <= 1.0:
        raise ValueError(
            f"レシピ換算係数は 0 より大きく 1 以下で指定する。"
            f"渡された値: {recipe_radiation_coef}"
        )
    if film_degradation_factor <= 0.0:
        raise ValueError(
            f"フィルム劣化係数は正の値で指定する。"
            f"渡された値: {film_degradation_factor}"
        )

    floor_area = HOUSE_SPECS[house]["floor_area_m2"]
    shape = _diurnal_shape(forecast)

    # --- 5分ごとの日射に展開する ---
    # 日積算 [MJ/m²] → 各ステップの瞬時日射 [W/m²]
    #   そのステップのエネルギー [J/m²] = 日積算 × 配分
    #   瞬時日射 [W/m²] = エネルギー ÷ 300秒
    total_j = forecast.sensor_radiation_mj / MJ_PER_J
    radiations = [total_j * w / SECONDS_PER_STEP for w in shape]

    # --- 気温も日射の形に合わせて揺らす ---
    #
    # 渡された気温・水蒸気圧は「日中の平均」なので、揺らしたあとの日中平均が
    # その値どおりになるように中心をそろえる。単純に (w/peak − 0.5) で
    # 振ると、正弦の日中平均が 2/π ≒ 0.64 なので 0.8℃ ほど上振れする。
    peak_shape = max(shape) if max(shape) > 0 else 1.0
    normalized = [w / peak_shape for w in shape]
    daytime = [n for n in normalized if n > 0.0]
    daytime_center = sum(daytime) / len(daytime) if daytime else 0.5

    # --- 日中と夜、それぞれの空気の状態を決める ---
    #
    # 湿り具合は水蒸気圧 ea で持つ。相対湿度は気温の従属変数なので、
    # 気温を揺らすと水蒸気が出入りしなくても動いてしまう（config.py 第9-2節）。
    day_vapor_kpa = (
        saturation_vapor_pressure_kpa(forecast.mean_temp_c)
        * forecast.mean_relative_humidity
        if forecast.mean_vapor_pressure_kpa is None
        else forecast.mean_vapor_pressure_kpa
    )
    if day_vapor_kpa <= 0.0:
        raise ValueError(
            f"日中の水蒸気圧が 0 以下: {day_vapor_kpa} kPa。"
            f"平均気温 {forecast.mean_temp_c} ℃・"
            f"相対湿度 {forecast.mean_relative_humidity} から作った値。"
        )

    # 夜は日中の振幅から外挿せず、暖房の効いた夜の水準をそのまま使う。
    # 外挿していたころは夜間蒸散が1日の15%を占めていた（実測は3〜10%）。
    night_temp = (
        forecast.mean_temp_c
        if forecast.night_temp_c is None
        else forecast.night_temp_c
    )
    night_vapor_kpa = (
        day_vapor_kpa
        if forecast.night_vapor_pressure_kpa is None
        else forecast.night_vapor_pressure_kpa
    )
    night_co2 = (
        forecast.mean_co2_ppm
        if forecast.night_co2_ppm is None
        else forecast.night_co2_ppm
    )

    # --- 5分ごとに群落のガス交換を解く ---
    transpiration_g_per_m2 = 0.0
    peak_transpiration_g_per_m2_s = 0.0
    net_photosynthesis_mmol_per_m2 = 0.0   # 1日の積算 [mmol/m²]
    for radiation, n in zip(radiations, normalized):
        if n <= 0.0:
            # 夜（日射がゼロの時刻）
            temp_c = night_temp
            vapor_kpa = night_vapor_kpa
            co2_ppm = night_co2
        else:
            temp_c = forecast.mean_temp_c + forecast.temp_amplitude_c * (
                n - daytime_center
            )
            vapor_kpa = day_vapor_kpa + DIURNAL_VAPOR_PRESSURE_AMPLITUDE_KPA * (
                n - daytime_center
            )
            co2_ppm = forecast.mean_co2_ppm

        # 飽和を超えた水蒸気圧は物理的にありえない（結露する）ので頭を押さえる。
        # 冬の朝など、気温が下がりきった時刻に起こりうる。
        saturation_kpa = saturation_vapor_pressure_kpa(temp_c)
        vapor_kpa = min(max(vapor_kpa, 0.01), saturation_kpa * 0.995)

        canopy = calculate_canopy_gas_exchange(
            air_temp_c=temp_c,
            vapor_pressure_kpa=vapor_kpa,
            solar_above_w_per_m2=radiation,
            air_co2_ppm=co2_ppm,
            wind_speed_m_per_s=wind_speed_m_per_s,
            pressure_kpa=LOCAL_PRESSURE_KPA,
            characteristic_length_m=LEAF_CHARACTERISTIC_LENGTH_M,
            lai=lai,
            extinction_coefficient=EXTINCTION_COEFFICIENT_K,
            solve_energy_balance=True,
        )
        rate = canopy.transpiration_mg_per_m2_s / 1000.0   # g/(m²·s)
        transpiration_g_per_m2 += rate * SECONDS_PER_STEP
        peak_transpiration_g_per_m2_s = max(peak_transpiration_g_per_m2_s, rate)
        net_photosynthesis_mmol_per_m2 += (
            canopy.net_photosynthesis_mmol_per_m2_s * SECONDS_PER_STEP
        )

    # 1日を1本のなめらかな曲線で代表させたぶんの目減りを戻す。
    # 雲のギザギザと日変化のとがりで稼いでいた水が、平滑化で消えている
    # （イェンセンの不等式。大きさの測り方は config.py 第9-3節）。
    # ★積算にだけ掛ける。上の peak は瞬時値なので掛けない。
    transpiration_g_per_m2 *= DIURNAL_SMOOTHING_CORRECTION

    transpiration_l = transpiration_g_per_m2 / GRAMS_PER_LITER

    # --- 光合成を目に見える単位へ ---
    # CO2 1 mol = 44 g。糖はCO2 6分子からグルコース(180)1分子ができる。
    co2_fixed_g = net_photosynthesis_mmol_per_m2 / 1000.0 * MOLAR_MASS_CO2_G_PER_MOL
    sugar_g = net_photosynthesis_mmol_per_m2 / 1000.0 * (GLUCOSE_G_PER_MOL / 6.0)

    # --- 日射の読みかえ（既存アプリと同じ約束）---
    # ★これは物理の換算ではなく、レシピと数字を合わせるための取り決め。
    #   「実質日射」は本当のハウス外日射ではない（config.py 第4-2節）。
    outside_mj = forecast.sensor_radiation_mj / recipe_radiation_coef
    effective_mj = outside_mj * film_degradation_factor

    # --- 潅水のめやす ---
    recommended = transpiration_l * (1.0 + leaching_fraction)
    water_per_10mj = (
        recommended / effective_mj * 10.0 if effective_mj > 0 else 0.0
    )

    # --- 土壌の水収支を1日ぶん回す ---
    settings = water_balance_settings or WaterBalanceSettings()
    theta_start = (
        FIELD_CAPACITY_WATER_CONTENT
        if start_water_content is None
        else start_water_content
    )
    balance = step_one_day(
        theta_start,
        DailyWaterInput(
            date="今日",
            irrigation_mm=recommended,
            potential_transpiration_mm=transpiration_l,
        ),
        settings,
    )

    # --- 気をつけることを言葉にする ---
    warnings: list[str] = []
    if balance.min_air_filled_porosity < 0.08:
        warnings.append(
            f"空気率が {balance.min_air_filled_porosity * 100:.1f}% まで下がる見込み。"
            f"根の呼吸には 10% 以上ほしいので、回数を分けて1回量を減らすとよい。"
        )
    if balance.drainage_fraction > 0.40:
        warnings.append(
            f"潅水の {balance.drainage_fraction * 100:.0f}% が根の下へ抜ける見込み。"
            f"塩を流す分を差し引いても多いので、減らす余地がある。"
        )
    if balance.hours_above_field_capacity >= 20.0:
        warnings.append(
            "ほぼ1日中、圃場容水量を超えたままになる見込み。"
            "土が乾く時間がないと根が更新されにくい。"
        )
    if forecast.sensor_radiation_mj < 3.0:
        warnings.append(
            "日射が少ない日。蒸散も少ないので、潅水を控えめにしないと"
            "土が過湿になりやすい。"
        )
    end_state = describe_soil_state(balance.end_water_content)
    if end_state.available_water < 0.3:
        warnings.append(
            f"夕方の有効水分度が {end_state.available_water:.2f} まで下がる見込み。"
            f"翌朝までに乾きすぎないか見ておくとよい。"
        )

    return IrrigationAdvice(
        sensor_radiation_mj=forecast.sensor_radiation_mj,
        outside_radiation_mj=outside_mj,
        effective_radiation_mj=effective_mj,
        transpiration_l_per_m2=transpiration_l,
        transpiration_l_per_house=transpiration_l * floor_area,
        peak_transpiration_mm_per_h=peak_transpiration_g_per_m2_s * 3600.0 / GRAMS_PER_LITER,
        co2_fixed_g_per_m2=co2_fixed_g,
        sugar_g_per_m2=sugar_g,
        sugar_kg_per_house=sugar_g * floor_area / 1000.0,
        water_per_10mj_l_per_m2=water_per_10mj,
        recommended_irrigation_l_per_m2=recommended,
        recommended_irrigation_l_per_house=recommended * floor_area,
        expected_drainage_l_per_m2=balance.drainage_mm,
        drainage_fraction=balance.drainage_fraction,
        end_pf=balance.end_pf,
        min_air_filled_porosity=balance.min_air_filled_porosity,
        hours_above_field_capacity=balance.hours_above_field_capacity,
        warnings=warnings,
    )
