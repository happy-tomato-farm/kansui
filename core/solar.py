"""
天文計算から、その日の日射量を求める。

【このモジュールの立ち位置】
core/ の他のモジュールが「測った値を使って計算する」のに対し、
ここは**測らなくても分かること**を計算する。太陽の位置は暦と緯度経度だけで
決まるので、センサーもデータも要らない。

    日付 + 緯度経度 → 太陽の高さ → 大気を通る距離 → 地表に届く日射

【出どころ】
Dropbox/作業記録/アプリ/日照予測.html に入っている式をそのまま移植した。
元は『生物環境物理学の基礎』第11章の式11.1〜11.13。
同じ入力に対して同じ数字が出るようにしてあるので、
あちらのアプリと突き合わせて確かめられる。

【大気透過率 τ（タウ）とは】
大気の外側に届いた日射のうち、大気を垂直に1回通り抜けたときに
残る割合。空気分子・水蒸気・塵による散乱と吸収で減る。

    澄んだ冬の晴天  τ ≒ 0.75〜0.80
    夏の晴天        τ ≒ 0.60〜0.70（水蒸気が多い）
    黄砂・花粉      τ が 0.05〜0.15 ほど下がる
    曇天            τ ≒ 0.45
    本降りの雨      τ 0.2 以下

日射量への効き方は「べき乗」なので、τ の差は太陽が低いときほど大きく響く。
太陽高度30°なら大気を2倍の距離通るので、τ の差が2乗で効く。

【★重要: ハウス内センサーとの関係】
この計算はハウスの**外**の日射を出す。ハウス内のセンサー基準に直すには
その日の透過率を掛ける。

    センサー基準 = 快晴計算 × 透過率(その日)

透過率は季節で動く。宇都宮気象台と突き合わせた実測で

    5月中旬 0.72  ／  11月中旬 0.64  （最大/最小 = 1.13 倍）

冬は太陽高度が低く、フィルムへの入射角が浅くなるぶん反射で失われるため。
年ごとの水準はフィルムの齢で下がる（2019年 0.76 → 2026年 0.67）。

【★もうひとつの「透過率」と混同しないこと】
液肥混入機レシピに渡す「実質日射」は センサー値 ÷ 0.65 で計算する。
この 0.65 は物理量ではなく、レシピと数字を合わせるためだけの
単位の取り決め（config.py の RECIPE_RADIATION_COEF）。
分子と分母で打ち消し合うので、実測に合っていなくてよい。
このモジュールが扱うのは物理のほうだけ。

【廃止した較正係数について】
以前は一律 1.20 を掛けていたが、ずれの向きが季節で逆だったため
一律では直らなかった。原因（τ の振幅不足と透過率の過小評価）を
直したので廃止した。経緯は config.py の第4-3節に残してある。
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass

from config import (
    CLEAR_SKY_TAU_COS_COEF,
    CLEAR_SKY_TAU_MEAN,
    CLEAR_SKY_TAU_SIN_COEF,
    CLOUDY_TAU,
    COVER_TRANSMITTANCE,
    DIFFUSE_FRACTION,
    SOLAR_CONSTANT_W_PER_M2,
    SOLAR_SITE,
    STANDARD_MERIDIAN_DEG,
    STANDARD_PRESSURE_KPA,
    TRANSMITTANCE_SEASON_COS,
    TRANSMITTANCE_SEASON_SIN,
)

#: 1日を何分割して積算するか。日照予測アプリに合わせて1分刻み。
MINUTES_PER_DAY = 1440

#: 国際標準大気で気圧を高度から求めるときの係数。
#:   Pa = 101.3 × (1 − 2.25577e-5 × h)^5.25588
PRESSURE_LAPSE_COEF = 2.25577e-5
PRESSURE_EXPONENT = 5.25588


# =============================================================================
# 1. 気圧と大気透過率
# =============================================================================

def pressure_from_altitude_kpa(altitude_m: float) -> float:
    """標高から大気圧を求める [kPa]（国際標準大気）。

    大気路程（太陽光が通る空気の量）を求めるのに使う。
    標高が高いほど空気が薄く、日射は強くなる。
    """
    base = 1.0 - PRESSURE_LAPSE_COEF * altitude_m
    if base <= 0.0:
        raise ValueError(
            f"標高が高すぎて大気圧が求まらない。渡された値: {altitude_m} m"
        )
    return STANDARD_PRESSURE_KPA * base ** PRESSURE_EXPONENT


def clear_sky_tau(day_of_year: int) -> float:
    """快晴時の大気透過率 τ を日付から求める。

    夏は水蒸気が多くて透過率が落ち、冬は乾いて上がる。
    その季節変化を1年周期の波で表したもの。

        τ = 0.7185 + 0.1429 × cos(2πJ/365) + 0.0400 × sin(2πJ/365)
        → 1月中旬 0.87 ／ 7月中旬 0.57 ／ 春秋 0.64〜0.80

    sin の項は「山が元日ちょうどではなく1月16日ごろにある」というずれを
    表している。cos だけでは山を元日からずらせないため。

    【宇都宮気象台の実測から当てはめた値】
    晴れた日341日について、その日の実測を再現する τ を1日ずつ解いて
    当てはめた。解けた τ は 0.57〜0.87 で、すべて物理的にありうる範囲。

    【この式は目安であって、その日の実際の τ ではない】
    黄砂・花粉・煙霧でこれより大きく下がる日がある。
    逆に、寒気が入って大気が澄んだ日はこれより上がる。
    日ごとのばらつきは標準偏差で 0.04 ほどある。
    アプリでは τ を手で動かせるようにしてあるので、
    この関数の値は「ふつうの快晴日の出発点」として使う。
    """
    _validate_day_of_year(day_of_year)
    angle = 2.0 * math.pi * day_of_year / 365.0
    return (
        CLEAR_SKY_TAU_MEAN
        + CLEAR_SKY_TAU_COS_COEF * math.cos(angle)
        + CLEAR_SKY_TAU_SIN_COEF * math.sin(angle)
    )


def cover_transmittance_for_date(
    day_of_year: int,
    base_transmittance: float = COVER_TRANSMITTANCE,
) -> float:
    """その日のハウス被覆の日射透過率を求める（外の日射 → センサー値）。

    ★これは実測にもとづく物理量である。
      「実質日射」を出すときに使う RECIPE_RADIATION_COEF とは別物。
      あちらは液肥混入機レシピと数字を合わせるためだけの取り決めで、
      物理的に正しい必要がない。混同しないこと。

        透過率(J) = base × exp(−0.0390 × cos(2πJ/365) + 0.0475 × sin(2πJ/365))
        → 5月中旬が最大、11月中旬が最小。最大/最小 = 1.13 倍

    【なぜ季節で動くのか】
    冬は太陽高度が低く、フィルムへの入射角が浅くなるぶん反射で失われる。
    骨材やハウス間の影も長くなる。快晴日だけで見ても同じ形が出るので、
    天気のせいではなく光学的なもの。

    Args:
        day_of_year: 通日（1〜366）
        base_transmittance: 年ごとの水準。フィルムの齢で動く。
            既定は直近3年の実測 0.678。張り替えたらここを測り直す。

    Returns:
        その日の透過率。
    """
    _validate_day_of_year(day_of_year)
    if base_transmittance <= 0.0:
        raise ValueError(
            f"透過率の水準は正の値で指定する。渡された値: {base_transmittance}"
        )
    angle = 2.0 * math.pi * day_of_year / 365.0
    season = math.exp(
        TRANSMITTANCE_SEASON_COS * math.cos(angle)
        + TRANSMITTANCE_SEASON_SIN * math.sin(angle)
    )
    return base_transmittance * season


def _validate_day_of_year(day_of_year: int) -> None:
    if not 1 <= day_of_year <= 366:
        raise ValueError(
            f"通日は 1〜366 で指定する。渡された値: {day_of_year}"
        )


def _validate_tau(tau: float) -> None:
    if not 0.0 < tau < 1.0:
        raise ValueError(
            f"大気透過率 τ は 0 と 1 のあいだで指定する。渡された値: {tau}"
        )


# =============================================================================
# 2. その日の太陽の動き（時刻に依らない値）
# =============================================================================

@dataclass(frozen=True)
class SolarDay:
    """その日の太陽の動き。時刻によらず、日付と場所だけで決まる。"""

    day_of_year: int                 #: 通日（1〜366）
    declination_rad: float           #: 赤緯 δ。太陽が赤道からどれだけ北にあるか
    equation_of_time_h: float        #: 均時差 ET [時]。時計と太陽のずれ
    solar_noon_h: float              #: 南中時刻 [時・日本標準時]
    half_daylength_h: float          #: 可照時間の半分 [時]

    @property
    def sunrise_h(self) -> float:
        """日の出の時刻 [時]。大気による屈折は考えていない。"""
        return self.solar_noon_h - self.half_daylength_h

    @property
    def sunset_h(self) -> float:
        """日の入りの時刻 [時]。"""
        return self.solar_noon_h + self.half_daylength_h

    @property
    def daylength_h(self) -> float:
        """可照時間 [時]。"""
        return 2.0 * self.half_daylength_h

    @property
    def declination_deg(self) -> float:
        """赤緯を度で返す。夏至 +23.4°、冬至 −23.4°。"""
        return math.degrees(self.declination_rad)


def solar_day(
    day_of_year: int,
    latitude_deg: float = SOLAR_SITE["LATITUDE_DEG"],
    longitude_deg: float = SOLAR_SITE["LONGITUDE_DEG"],
) -> SolarDay:
    """その日の太陽の動きを求める。

    Args:
        day_of_year: 通日（1〜366）
        latitude_deg: 緯度 [°N]
        longitude_deg: 経度 [°E]

    【式の中身】
    式11.4 均時差 ET
        地球の公転軌道が楕円で、自転軸が傾いているために、
        太陽が南中する時刻は日によって最大±16分ずれる。
    式11.3 南中時刻
        t0 = 12 − 経度補正 − 均時差
        日本標準時は東経135°基準なので、それより東なら早く南中する。
    式11.2 赤緯 δ
        太陽が天の赤道からどれだけ北（南）にあるか。季節そのもの。
    """
    _validate_day_of_year(day_of_year)
    if not -90.0 < latitude_deg < 90.0:
        raise ValueError(
            f"緯度は −90 〜 90 の範囲で指定する。渡された値: {latitude_deg}"
        )

    # 経度補正 [時]。標準子午線から1°離れるごとに4分ずれる。
    longitude_correction_h = (longitude_deg - STANDARD_MERIDIAN_DEG) / 15.0

    # 式11.4 均時差
    angle_rad = math.radians(279.575 + 0.9856 * day_of_year)
    equation_of_time_h = (
        -104.7 * math.sin(angle_rad)
        + 596.2 * math.sin(2.0 * angle_rad)
        + 4.3 * math.sin(3.0 * angle_rad)
        - 12.7 * math.sin(4.0 * angle_rad)
        - 429.3 * math.cos(angle_rad)
        - 2.0 * math.cos(2.0 * angle_rad)
        + 19.3 * math.cos(3.0 * angle_rad)
    ) / 3600.0

    # 式11.3 南中時刻
    solar_noon_h = 12.0 - longitude_correction_h - equation_of_time_h

    # 式11.2 赤緯
    sin_declination = 0.39785 * math.sin(math.radians(
        278.97 + 0.9856 * day_of_year
        + 1.9165 * math.sin(math.radians(356.6 + 0.9856 * day_of_year))
    ))
    sin_declination = min(max(sin_declination, -1.0), 1.0)
    declination_rad = math.asin(sin_declination)

    # 可照時間の半分。cosψ = 0 になる時刻（太陽が地平線に来る）から求める。
    hour_angle_cos = -math.tan(math.radians(latitude_deg)) * math.tan(declination_rad)
    if hour_angle_cos >= 1.0:
        half_daylength_h = 0.0        # 極夜。日本では起きない。
    elif hour_angle_cos <= -1.0:
        half_daylength_h = 12.0       # 白夜。日本では起きない。
    else:
        half_daylength_h = math.degrees(math.acos(hour_angle_cos)) / 15.0

    return SolarDay(
        day_of_year=day_of_year,
        declination_rad=declination_rad,
        equation_of_time_h=equation_of_time_h,
        solar_noon_h=solar_noon_h,
        half_daylength_h=half_daylength_h,
    )


def solar_day_from_date(
    date: dt.date,
    latitude_deg: float = SOLAR_SITE["LATITUDE_DEG"],
    longitude_deg: float = SOLAR_SITE["LONGITUDE_DEG"],
) -> SolarDay:
    """日付（datetime.date）から太陽の動きを求める。ふだんはこちらを使う。"""
    return solar_day(date.timetuple().tm_yday, latitude_deg, longitude_deg)


# =============================================================================
# 3. 指定時刻の日射
# =============================================================================

@dataclass(frozen=True)
class RadiationComponents:
    """ある瞬間の日射の内訳 [W/m²]。"""

    cos_zenith: float            #: cosψ。太陽が真上なら1、地平線で0
    altitude_deg: float          #: 太陽高度 [°]
    air_mass: float              #: 大気路程 m。真上を1とした空気の厚み
    direct_normal_w: float       #: 法線面直達日射 Sp（太陽に正対した面）
    direct_w: float              #: 水平面直達日射 Sb
    diffuse_w: float             #: 水平面散乱日射 Sd（空全体からの青空光）
    global_w: float              #: 全天日射 St = Sb + Sd


def instantaneous_radiation(
    hour: float,
    day: SolarDay,
    tau: float,
    latitude_deg: float = SOLAR_SITE["LATITUDE_DEG"],
    pressure_kpa: float | None = None,
) -> RadiationComponents:
    """指定時刻の全天日射 [W/m²] を求める。

    Args:
        hour: 時刻 [時・日本標準時]。12.5 なら 12時30分。
        day: solar_day() の結果
        tau: 大気透過率
        latitude_deg: 緯度 [°N]
        pressure_kpa: 現地の大気圧。省くと SOLAR_SITE の標高から求める。

    【式の中身】
    式11.1 cosψ = sinΦ sinδ + cosΦ cosδ cos(15(t − t0))
        太陽の天頂角。緯度Φ・赤緯δ・時角から幾何学的に決まる。
    式11.12 m = Pa / (101.3 cosψ)
        大気路程。太陽が低いほど空気を長く通る。
        真上なら1、高度30°なら約2、地平線近くでは10を超える。
    式11.11 Sp = Spo τ^m
        直達日射。τ を m 乗するのがポイントで、
        空気を2倍通れば残るのは τ² になる。
    式11.13 Sd = 0.3 (1 − τ^m) Spo cosψ
        散乱日射。大気で失われた分（1 − τ^m）の3割が、
        青空からの光として地表に届くという近似。
    """
    _validate_tau(tau)
    if pressure_kpa is None:
        pressure_kpa = pressure_from_altitude_kpa(SOLAR_SITE["ALTITUDE_M"])

    # 式11.1 太陽の天頂角
    hour_angle_deg = 15.0 * (hour - day.solar_noon_h)
    cos_zenith = (
        math.sin(math.radians(latitude_deg)) * math.sin(day.declination_rad)
        + math.cos(math.radians(latitude_deg)) * math.cos(day.declination_rad)
        * math.cos(math.radians(hour_angle_deg))
    )
    cos_zenith = min(max(cos_zenith, -1.0), 1.0)
    altitude_deg = math.degrees(math.asin(cos_zenith))

    # 太陽が地平線の下なら日射はゼロ
    if cos_zenith <= 0.0:
        return RadiationComponents(
            cos_zenith=cos_zenith, altitude_deg=altitude_deg, air_mass=float("inf"),
            direct_normal_w=0.0, direct_w=0.0, diffuse_w=0.0, global_w=0.0,
        )

    # 式11.12 大気路程
    air_mass = pressure_kpa / (STANDARD_PRESSURE_KPA * cos_zenith)
    transmitted = tau ** air_mass

    # 式11.11 / 11.8 直達日射
    direct_normal_w = SOLAR_CONSTANT_W_PER_M2 * transmitted
    direct_w = direct_normal_w * cos_zenith

    # 式11.13 散乱日射
    diffuse_w = (
        DIFFUSE_FRACTION * (1.0 - transmitted)
        * SOLAR_CONSTANT_W_PER_M2 * cos_zenith
    )

    return RadiationComponents(
        cos_zenith=cos_zenith,
        altitude_deg=altitude_deg,
        air_mass=air_mass,
        direct_normal_w=direct_normal_w,
        direct_w=direct_w,
        diffuse_w=diffuse_w,
        global_w=direct_w + diffuse_w,
    )


def extraterrestrial_radiation_w(hour: float, day: SolarDay,
                                 latitude_deg: float = SOLAR_SITE["LATITUDE_DEG"]) -> float:
    """大気の外側での水平面日射 [W/m²]。**物理の天井**。

    大気も屋根もない宇宙空間での値なので、地表で測った値がこれを
    超えることはありえない。センサーの点検に使う。
    """
    hour_angle_deg = 15.0 * (hour - day.solar_noon_h)
    cos_zenith = (
        math.sin(math.radians(latitude_deg)) * math.sin(day.declination_rad)
        + math.cos(math.radians(latitude_deg)) * math.cos(day.declination_rad)
        * math.cos(math.radians(hour_angle_deg))
    )
    return SOLAR_CONSTANT_W_PER_M2 * max(cos_zenith, 0.0)


# =============================================================================
# 4. 1日の積算
# =============================================================================

def daily_radiation_mj(
    date: dt.date,
    tau: float | None = None,
    latitude_deg: float = SOLAR_SITE["LATITUDE_DEG"],
    longitude_deg: float = SOLAR_SITE["LONGITUDE_DEG"],
    altitude_m: float = SOLAR_SITE["ALTITUDE_M"],
) -> float:
    """ハウスの**外**の1日の積算日射量 [MJ/m²]。

    Args:
        date: 日付
        tau: 大気透過率。省くとその日の快晴値を使う。
        latitude_deg / longitude_deg / altitude_m: 場所

    1分刻みで全天日射を求めて足し上げる（日照予測アプリと同じ刻み）。
    """
    day = solar_day_from_date(date, latitude_deg, longitude_deg)
    if tau is None:
        tau = clear_sky_tau(day.day_of_year)
    _validate_tau(tau)
    pressure_kpa = pressure_from_altitude_kpa(altitude_m)

    total_j = 0.0
    for minute in range(1, MINUTES_PER_DAY + 1):
        radiation = instantaneous_radiation(
            minute / 60.0, day, tau, latitude_deg, pressure_kpa
        )
        total_j += radiation.global_w * 60.0      # 各分の代表値 × 60秒
    return total_j * 1.0e-6


def sensor_basis_radiation_mj(
    date: dt.date,
    tau: float | None = None,
    cover_transmittance: float = COVER_TRANSMITTANCE,
) -> float:
    """ハウス**内**センサー基準に直した1日の積算日射量 [MJ/m²]。

        センサー基準 = 外の日射 × その日の透過率

    透過率は季節で動く（5月 0.72 ／ 11月 0.64）。
    cover_transmittance にはその年の水準を渡し、季節の形は中で掛ける。

    【以前あった較正係数 1.20 は廃止した】
    快晴計算が低く出るのを一律の係数で埋めていたが、ずれの向きが
    季節で逆だったため、冬は埋まっても秋は3割も過大になっていた。
    τ と透過率を実測から当て直したので、この係数は要らなくなった。
    経緯は config.py の第4-3節に残してある。
    """
    if cover_transmittance <= 0.0:
        raise ValueError(
            f"被覆の透過率は正の値で指定する。渡された値: {cover_transmittance}"
        )
    day_of_year = date.timetuple().tm_yday
    transmittance = cover_transmittance_for_date(day_of_year, cover_transmittance)
    outside_mj = daily_radiation_mj(date, tau)
    return outside_mj * transmittance


def clear_sky_sensor_mj(date: dt.date, **kwargs) -> float:
    """快晴の日のセンサー基準日射 [MJ/m²]。その日の上限のめやす。"""
    return sensor_basis_radiation_mj(date, tau=None, **kwargs)


def cloudy_sensor_mj(date: dt.date, **kwargs) -> float:
    """曇天の日のセンサー基準日射 [MJ/m²]。

    τ = 0.45 は「よく曇った日」であって最低値ではない。
    実測では半分の日がこれを下回る（本降りの雨なら 1/5 まで落ちる）。
    """
    return sensor_basis_radiation_mj(date, tau=CLOUDY_TAU, **kwargs)
