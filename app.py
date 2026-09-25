"""
今日の潅水量のめやす — ブラウザで動く画面。

【起動のしかた】
このフォルダで、コマンドプロンプト（または PowerShell）から

    py -m streamlit run app.py

と打つ。ブラウザが自動で開く。閉じるときは黒い画面で Ctrl+C。

【この画面の役割】
朝、日射予測と葉の枚数を入れると、潅水量のめやすが出る。
いちばん大事なのは **「10MJあたり潅水量」** で、これを
液肥混入機レシピのアプリにそのまま入れれば、潅水時間が決まる。

【画面の作り】
スマホでも使えるよう、サイドバーは使わず本体に縦に並べている。
上から順に「よく触るもの → ときどき触るもの → 結果 → 参考」。
よく触る6つ（日付・葉の枚数・葉の面積・上乗せ・気温・湿度）を最上部に置き、
残りはタブの中にしまってある。

【計算の立ち位置】
潅水量は物理式から出している（エネルギー収支・気孔コンダクタンス・
Beer則・Campbell 式9.20/9.22）。過去の潅水記録はモデルには使っておらず、
答え合わせのために画面に並べているだけ。
"""

from __future__ import annotations

import datetime as dt

import streamlit as st

from config import (
    CLOUDY_TAU,
    COVER_TRANSMITTANCE,
    DRIP_WETTED_FRACTION,
    HOUSE_SPECS,
    LEAF_AREA_PER_LEAF_M2,
    RECIPE_RADIATION_COEF,
    ROOT_SYSTEM_MAX_UPTAKE_MM_PER_DAY,
    ROOT_ZONE_DEPTH_M,
    SOLAR_SITE,
    STOMATA_SCALE_FACTOR,
    WIND_SPEED_M_PER_S,
)
from core.advisor import (
    MONTHS_WITH_WEAK_DATA,
    advise,
    forecast_from_date,
    normal_radiation_mj,
)
from core.canopy import leaf_count_per_m2
from core.irrigation_history import describe_source, lookup
from core.psychrometry import saturation_vapor_pressure_kpa
from core.soil import (
    FIELD_CAPACITY_WATER_CONTENT,
    WILTING_POINT_WATER_CONTENT,
    describe_soil_state,
    pf_from_potential,
    potential_from_pf,
    water_content_from_potential,
)
from core.solar import (
    clear_sky_tau,
    cover_transmittance_for_date,
    sensor_basis_radiation_mj,
    solar_day_from_date,
)
from core.water_balance import WaterBalanceSettings, steady_state

st.set_page_config(page_title="今日の潅水量のめやす", page_icon="💧", layout="wide")

#: 1アールは100m²。10アール＝1000m²。
M2_PER_10A = 1000.0


def panel(title: str, body: str, accent: str, background: str) -> None:
    """色のついた四角で囲んで表示する。

    st.info などより目立たせたい、結果の要点だけに使う。
    Streamlit の素の部品では枠の色を変えられないので HTML を書いている。
    """
    st.markdown(
        f"""
        <div style="
            border-left: 6px solid {accent};
            background: {background};
            border-radius: 6px;
            padding: 0.9rem 1.1rem;
            margin-bottom: 0.8rem;
        ">
          <div style="
              font-weight: 700; font-size: 1.0rem;
              color: {accent}; margin-bottom: 0.45rem;
          ">{title}</div>
          <div style="font-size: 0.95rem; line-height: 1.7;">{body}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


# =============================================================================
# 1. よく触る入力（最上部）
# =============================================================================

st.title("💧 今日の潅水量のめやす")

today = dt.date.today()

row = st.columns([1.2, 1, 1.4])
with row[0]:
    target_date = st.date_input("① 日付", today)
with row[1]:
    house = st.selectbox("ハウス", list(HOUSE_SPECS.keys()))
with row[2]:
    preset = st.radio(
        "きょうの空",
        ["快晴", "平年なみ", "曇天", "自分で入れる"],
        index=0,
        horizontal=True,
        help=(
            "快晴・曇天は太陽の位置から計算した値。"
            "平年なみは実測の同じ日付の平均。"
            "気象予報の日射量があるなら「自分で入れる」が一番当たる。"
        ),
    )

month = target_date.month
day_of_year = target_date.timetuple().tm_yday
planned_leaves = leaf_count_per_m2(target_date)

# --- 詳しい設定はたたんでおく ---
#
# 【なぜスライダーをやめたか】（2026-09-25）
# スマホで画面を上下にスワイプすると、指が触れたスライダーが動いてしまう。
# 気づかないうちに値が変わって、出てくる水の量がずれる。
# 数値入力欄なら、タップして数字を打つまで動かない。
#
# 【なぜタブをやめて expander にしたか】
# タブは「常にどれか1つが開いている」ので、触るつもりのない設定が
# いつも画面に出ている。expander なら閉じておける＝誤って触らない。
#
# ★ここで決めた値が下の日射計算に要るので、位置は画面の上のまま。
#   閉じていれば場所は取らない。

with st.expander("☀ 日射の細かい設定", expanded=False):
    st.caption("ふだんは触らなくてよい。空がふつうでない日だけ τ を動かす。")

    st.markdown("**物理の設定（外の日射 → ハウス内センサー値）**")
    cols = st.columns(2)
    with cols[0]:
        tau = st.number_input(
            "大気透過率 τ", min_value=0.20, max_value=0.95,
            value=float(clear_sky_tau(day_of_year)), step=0.01, format="%.2f",
            key=f"tau_{target_date}",
            help=(
                "空の澄み具合。0.20〜0.95 の範囲。"
                "初期値は日付から決まる快晴値で、"
                "宇都宮気象台の実測341日から当てはめた。"
                "黄砂・花粉・煙霧の日は 0.05〜0.15 下げる。"
            ),
        )
    with cols[1]:
        transmittance_base = st.number_input(
            "ハウス被覆の透過率（年の水準）", min_value=0.50, max_value=0.90,
            value=float(COVER_TRANSMITTANCE), step=0.005, format="%.3f",
            help=(
                "0.50〜0.90 の範囲。フィルムの齢で決まる水準。"
                "実測では 2019年 0.76 → 2026年 0.67。"
                "★2026年6月に張り替えた後のデータはまだ無いので、"
                "快晴の日が10日ほどたまったら測り直すこと。"
                "季節の変化は下に書いたとおり自動で掛かる。"
            ),
        )
    transmittance_today = cover_transmittance_for_date(
        day_of_year, transmittance_base
    )
    st.caption(
        f"この日の快晴 τ は **{clear_sky_tau(day_of_year):.3f}**"
        f"（1月 0.87 ／ 7月 0.57）。曇天なら 0.45 あたり。　／　"
        f"この日の透過率は **{transmittance_today:.3f}**"
        f"（水準 {transmittance_base:.3f} × 季節 "
        f"{transmittance_today / transmittance_base:.3f}）。"
        f"5月中旬が最大、11月中旬が最小で 1.13 倍の幅がある。"
    )

    st.divider()
    st.markdown("**液肥混入機レシピと合わせる係数（物理量ではない）**")
    cols2 = st.columns(2)
    with cols2[0]:
        recipe_coef = st.number_input(
            "レシピの日射センサー係数", min_value=0.50, max_value=0.90,
            value=float(RECIPE_RADIATION_COEF), step=0.01, format="%.2f",
            help=(
                "0.50〜0.90 の範囲。「実質日射」＝ センサー値 ÷ この係数。"
                "作業日誌の「センサー補正係数 6.5」に対応する。"
                "★実測の透過率に合わせる必要はない。"
                "レシピ側と同じ値であることだけが大事。"
            ),
        )
    with cols2[1]:
        film_factor = st.number_input(
            "フィルム劣化係数", min_value=0.70, max_value=1.00,
            value=float(HOUSE_SPECS[house]["film_degradation_factor"]),
            step=0.01, format="%.2f",
            help=(
                "0.70〜1.00 の範囲。これもレシピ側の film と同じ値にすること。"
                "ずれると、実際に出る水の量がその比率でずれる。"
            ),
        )
    st.caption(
        f"液肥混入機レシピの「{house}」の設定は "
        f"film = {HOUSE_SPECS[house]['film_degradation_factor']:.2f}、"
        f"日射センサーの係数 = {RECIPE_RADIATION_COEF * 10:.1f}。"
        f"　この2つは単位の取り決めなので、分子（潅水量）と分母（実質日射）で"
        f"打ち消し合う。実測に合っていなくてよい。"
    )

    st.divider()
    st.markdown("**日射をどちらの基準で入れるか**")
    radiation_basis = st.radio(
        "日射をどちらの基準で入れるか",
        ["ハウスの中（センサー値）", "気象予報の全天日射量"],
        label_visibility="collapsed",
        help=(
            "ハウスのセンサーは屋根の下にあるので、外より小さい値が出る。"
            "気象庁やウェザーニュースの予報を使うときは下を選ぶ。"
            "選ぶと、下の入力欄の単位が切り替わる。"
        ),
    )
    inside_basis = radiation_basis.startswith("ハウスの中")

with st.expander("🪴 土と群落の設定", expanded=False):
    st.caption("土と群落の前提。実測できていない値が多いので、感度を見るのに使う。")
    cols = st.columns(2)
    with cols[0]:
        wind = st.number_input(
            "群落内の風速 [m/s]", min_value=0.1, max_value=3.0,
            value=float(WIND_SPEED_M_PER_S), step=0.1, format="%.1f",
            help=(
                "0.1〜3.0 の範囲。ダクト送風を24時間しているので 0.2 を既定にしている。"
                "サイドカーテンを開けると上がるが、0.2→2.0 でも蒸散は1.21倍止まり。"
            ),
        )
        root_depth = st.number_input(
            "根群域の深さ [m]", min_value=0.20, max_value=0.80,
            value=float(ROOT_ZONE_DEPTH_M), step=0.05, format="%.2f",
            help="0.20〜0.80 の範囲。実測できていないが、水の配分にはほとんど効かない。",
        )
        start_pf = st.number_input(
            "朝の土壌 pF", min_value=1.2, max_value=3.5,
            value=1.8, step=0.1, format="%.1f",
            help="1.2〜3.5 の範囲。1.8 が圃場容水量。"
                 "前日の潅水が効いていればこのあたり。",
        )
    with cols[1]:
        wetted = st.number_input(
            "点滴で濡れる床面積の割合", min_value=0.3, max_value=1.0,
            value=float(DRIP_WETTED_FRACTION), step=0.05, format="%.2f",
            help="0.3〜1.0 の範囲。",
        )
        max_uptake = st.number_input(
            "根系の最大吸水速度 [mm/日]", min_value=3.0, max_value=12.0,
            value=float(ROOT_SYSTEM_MAX_UPTAKE_MM_PER_DAY), step=0.5, format="%.1f",
            help="3.0〜12.0 の範囲。",
        )
    st.caption(
        f"気孔コンダクタンス倍率は {STOMATA_SCALE_FACTOR}（文献に合わせた値）で固定。"
    )

# --- 日射の3つの目安 ---
# 外の日射 → センサー値 は、その日の透過率（季節変化つき）だけで換算する。
# 以前あった較正係数 1.20 は廃止した（config.py 第4-3節）。
outside_to_sensor = transmittance_today
normal_mj = normal_radiation_mj(day_of_year)
clear_mj = sensor_basis_radiation_mj(
    target_date, tau=tau, cover_transmittance=transmittance_base,
)
cloudy_mj = sensor_basis_radiation_mj(
    target_date, tau=CLOUDY_TAU, cover_transmittance=transmittance_base,
)
preset_values = {"快晴": clear_mj, "平年なみ": normal_mj, "曇天": cloudy_mj}
default_inside = preset_values.get(preset, clear_mj)

# =============================================================================
# 1-2. 毎朝触るところ（たたまずに出しておく）
# =============================================================================

st.divider()

# --- 日射の入力 ---
# 基準の選び分けは「日射の細かい設定」の中に移してある。
# ふだんはハウス内センサー基準のままなので、画面に出しておく必要がない。
bottom = st.columns([1.4, 1])
with bottom[0]:
    if inside_basis:
        entered = st.number_input(
            "ハウス内の日射 [MJ/m²]",
            min_value=0.0, max_value=30.0,
            value=float(round(default_inside, 1)), step=0.5,
            key=f"rad_in_{target_date}_{preset}",
            help="液肥混入機レシピのアプリに入れるのと同じ値。",
        )
        radiation = entered
    else:
        entered = st.number_input(
            "気象予報の全天日射量 [MJ/m²]",
            min_value=0.0, max_value=45.0,
            value=float(round(default_inside / outside_to_sensor, 1)), step=0.5,
            key=f"rad_out_{target_date}_{preset}",
            help="気象庁やウェザーニュースなどの全天日射量の予報値。",
        )
        radiation = entered * outside_to_sensor
with bottom[1]:
    st.caption(
        f"この日の目安（ハウス内）\n\n"
        f"快晴 **{clear_mj:.1f}** ／ 平年 **{normal_mj:.1f}** ／ "
        f"曇天 **{cloudy_mj:.1f}** MJ/m²"
    )

# --- 葉と潅水（ふだんは日付から引いた値のままでよいので、たたんでおく）---
with st.expander(
    f"🌿 葉と潅水の設定"
    f"（いまは {planned_leaves:.1f} 枚/m² ／ "
    f"LAI {planned_leaves * LEAF_AREA_PER_LEAF_M2:.2f} ／ 上乗せ 0%）",
    expanded=False,
):
    st.caption(
        "葉の枚数は作業計画シートの週ごとの目標を日付から引いている。"
        "実際の枚数が違うときや、塩を流すために水を増やすときだけ開く。"
    )
    top = st.columns(3)
    with top[0]:
        st.markdown("**② 葉の枚数 [枚/m²]**")
        leaves_per_m2 = st.number_input(
            "葉の枚数", min_value=8.0, max_value=48.0,
            value=float(round(planned_leaves, 1)), step=0.5, format="%.1f",
            label_visibility="collapsed",
            key=f"leaves_{target_date}",
            help=(
                "8〜48 の範囲。作業計画シート「葉枚数管理／必要枚数/㎡」の"
                "週ごとの目標値。初期値は日付から引いている。"
            ),
        )
        st.caption(
            f"作業計画（週{target_date.isocalendar()[1]}）の目標は "
            f"**{planned_leaves:.1f} 枚/m²**。冬16枚→春40枚と倍以上動く。"
        )
    with top[1]:
        st.markdown("**③ 葉1枚の面積 [m²]**")
        leaf_area = st.number_input(
            "葉1枚の面積", min_value=0.06, max_value=0.22,
            value=float(LEAF_AREA_PER_LEAF_M2), step=0.01, format="%.2f",
            label_visibility="collapsed",
            help="0.06〜0.22 の範囲。★実測してほしい値。結果にいちばん効く。",
        )
        lai = leaves_per_m2 * leaf_area
        st.caption(f"**LAI = {lai:.2f}**（葉の枚数 × 葉1枚の面積）")
        if lai > 5.0:
            st.caption("⚠ LAI 5 超。葉1枚の面積が過大でないか確かめること。")
    with top[2]:
        st.markdown("**④ 塩を流すための上乗せ [%]**")
        leaching_percent = st.number_input(
            "上乗せ", min_value=0, max_value=150, value=0, step=5,
            label_visibility="collapsed",
            help=(
                "0〜150 の範囲。蒸散量に対して何％多く入れるか。既定は 0%。"
                "排液の EC が上がってきたら増やす。"
            ),
        )
        if leaching_percent == 0:
            st.caption(
                "0% ＝ 蒸散量とちょうど同じ量。"
                "実際の潅水記録はこれより3〜4割多い（塩を流すぶん）。"
            )
        else:
            st.caption(f"蒸散量の {1 + leaching_percent / 100:.2f} 倍を入れる。")

# 気温・湿度は日射が決まってから推定するので、ここで計算
auto = forecast_from_date(radiation, target_date)

# --- 気温・湿度・CO2 は推定値のままでよいので、たたんでおく ---
with st.expander(
    f"🌡 気温・湿度・CO₂ を打ち替える"
    f"（いまは推定値 {auto.mean_temp_c:.1f} ℃ ／ "
    f"{auto.mean_relative_humidity * 100:.0f} % ／ "
    f"{auto.mean_co2_ppm:.0f} ppm）",
    expanded=False,
):
    st.caption(
        "日付と日射から推定した値が入っている。"
        "ハウスの実測があるときだけ打ち替える。"
        "日付か「きょうの空」を変えると推定し直される。"
    )
    weather = st.columns(3)
    with weather[0]:
        st.markdown("**⑤ 日中の平均気温 [℃]**")
        manual_temp = st.number_input(
            "気温", min_value=5.0, max_value=40.0,
            value=float(round(auto.mean_temp_c, 1)), step=0.5, format="%.1f",
            label_visibility="collapsed",
            key=f"temp_{target_date}_{preset}",
            help="5〜40 の範囲。",
        )
    with weather[1]:
        st.markdown("**⑥ 日中の平均湿度 [%]**")
        manual_rh = st.number_input(
            "湿度", min_value=30, max_value=95,
            value=int(round(auto.mean_relative_humidity * 100)), step=1,
            label_visibility="collapsed",
            key=f"rh_{target_date}_{preset}",
            help="30〜95 の範囲。",
        )
    with weather[2]:
        st.markdown("**日中の平均CO₂ [ppm]**")
        manual_co2 = st.number_input(
            "CO2", min_value=300.0, max_value=1500.0,
            value=float(round(auto.mean_co2_ppm)), step=10.0, format="%.0f",
            label_visibility="collapsed",
            key=f"co2_{target_date}_{preset}",
            help="300〜1500 の範囲。",
        )

# 日射の読みかえ結果
# 外の日射は物理の透過率で、実質日射はレシピの取り決めで割る。別の仕事なので別の値。
forecast_outside_mj = radiation / transmittance_today
effective_mj = radiation / recipe_coef * film_factor

if radiation > clear_mj * 1.05:
    st.warning(
        f"入れた日射 {radiation:.1f} MJ/m² は、この日の快晴の計算値 "
        f"{clear_mj:.1f} MJ/m² を {(radiation / clear_mj - 1) * 100:.0f}% "
        f"上回っている。τ を上げるか、値を見直すこと。"
    )


# =============================================================================
# 2. 計算
# =============================================================================

forecast = forecast_from_date(
    radiation, target_date,
    temp_c=manual_temp,
    relative_humidity=manual_rh / 100.0,
    co2_ppm=manual_co2,
)
water_settings = WaterBalanceSettings(
    root_zone_depth_m=root_depth,
    wetted_fraction=wetted,
    max_uptake_mm_per_day=max_uptake,
)
advice = advise(
    forecast,
    house=house,
    lai=lai,
    start_water_content=water_content_from_potential(potential_from_pf(start_pf)),
    leaching_fraction=leaching_percent / 100.0,
    wind_speed_m_per_s=wind,
    water_balance_settings=water_settings,
    recipe_radiation_coef=recipe_coef,
    film_degradation_factor=film_factor,
)
# ★入れた日射を渡して、同じくらいの明るさだった日だけを選んでもらう。
# 渡さないと天気を問わない平均になり、快晴の計算値と曇りこみの平均を
# 並べることになってしまう（10月なら 1.7 倍も違う）。
history = lookup(target_date, house, sensor_radiation_mj=radiation)

# 飽差（VPD）＝ 空気があとどれだけ水蒸気を受け取れるかの余力。蒸散の駆動力。
# 日中と夜で別々に出す。夜は暖房の設定温度の水準を使っている（config 第9-4節）。
day_vpd_kpa = (
    saturation_vapor_pressure_kpa(forecast.mean_temp_c)
    - forecast.mean_vapor_pressure_kpa
)
night_saturation_kpa = saturation_vapor_pressure_kpa(forecast.night_temp_c)
night_vpd_kpa = night_saturation_kpa - forecast.night_vapor_pressure_kpa
night_rh = forecast.night_vapor_pressure_kpa / night_saturation_kpa


# =============================================================================
# 3. 結果の要点（色のついた囲み）
# =============================================================================

st.divider()

panel(
    "液肥混入機レシピに入れる数字",
    f"""
    <div style="font-size:2.2rem; font-weight:800; line-height:1.2;">
      {advice.water_per_10mj_l_per_m2:.2f} <span style="font-size:1.1rem">L/m²</span>
    </div>
    <div style="margin:0.3rem 0 0.7rem 0;">「10MJあたり潅水量」の欄に入れる</div>
    <table style="width:100%; border-collapse:collapse;">
      <tr><td>日射予測の欄に入れる値</td>
          <td style="text-align:right"><b>{advice.sensor_radiation_mj:.1f}</b> MJ/m²</td></tr>
    </table>
    <div style="margin-top:0.6rem; font-size:0.85rem; opacity:0.85;">
      日射予測は<b>ハウス内センサー基準</b>のまま入れる。
      ÷透過率 ×フィルム係数 はレシピ側がやる。
    </div>
    """,
    accent="#1b7a3d",
    background="#eaf6ee",
)

# 前提条件はふだん見なくてよいので、たたんでおく。
# 数字が思ったのと違うときに、ここを開いて元をたどる。
with st.expander("📋 この日の前提条件", expanded=False):
    panel(
        "この日の前提条件",
        f"""
        <table style="width:100%; border-collapse:collapse;">
          <tr><td>実質日射（10MJあたりの分母）</td>
              <td style="text-align:right"><b>{advice.effective_radiation_mj:.1f}</b> MJ/m²</td></tr>
          <tr><td>日中の平均気温</td>
              <td style="text-align:right"><b>{forecast.mean_temp_c:.1f}</b> ℃</td></tr>
          <tr><td>日中の平均湿度</td>
              <td style="text-align:right"><b>{forecast.mean_relative_humidity * 100:.0f}</b> %
              （飽差 <b>{day_vpd_kpa:.2f}</b> kPa）</td></tr>
          <tr><td>日中の平均CO₂</td>
              <td style="text-align:right"><b>{forecast.mean_co2_ppm:.0f}</b> ppm</td></tr>
          <tr><td>夜（日没〜日の出）</td>
              <td style="text-align:right"><b>{forecast.night_temp_c:.1f}</b> ℃ ／
              <b>{night_rh * 100:.0f}</b> %
              （飽差 <b>{night_vpd_kpa:.2f}</b> kPa）</td></tr>
          <tr><td>葉の枚数 ／ LAI</td>
              <td style="text-align:right"><b>{leaves_per_m2:.1f}</b> 枚/m² ／ <b>{lai:.2f}</b></td></tr>
        </table>
        """,
        accent="#1f6f8b",
        background="#eaf4f8",
    )


# =============================================================================
# 4. その日の見通し
# =============================================================================

st.subheader("その日の見通し")

# 日射の2つは「この日の前提条件」から移してきた。
# どの明るさを前提にした数字なのかが、開かなくても分かるようにするため。
radiation_columns = st.columns(2)
radiation_columns[0].metric(
    "気象予報の全天日射量", f"{forecast_outside_mj:.1f} MJ/m²",
    help=(
        "ハウスの外の日射。気象庁やウェザーニュースの予報値と"
        "直接くらべられる。ハウス内の値 ÷ その日の透過率。"
    ),
)
radiation_columns[1].metric(
    "ハウスの中の日射", f"{advice.sensor_radiation_mj:.1f} MJ/m²",
    help=(
        "作物が実際に浴びる光。蒸散の計算に使っているのはこの値。"
        "液肥混入機レシピの「日射予測」の欄にもこの値を入れる。"
    ),
)

columns = st.columns(4)
columns[0].metric("予測蒸散量", f"{advice.transpiration_l_per_m2:.2f} L/m²")
columns[1].metric(
    "合計潅水量", f"{advice.recommended_irrigation_l_per_m2:.2f} L/m²",
    help="蒸散量に上乗せ分を足したもの。1日ぶんの合計。",
)
columns[2].metric(
    "10アールあたり",
    f"{advice.recommended_irrigation_l_per_m2 * M2_PER_10A / 1000:.2f} t",
    help="10アール＝1000m²。L/m² の値がそのまま t/10a になる。",
)
columns[3].metric("正午前後の蒸散", f"{advice.peak_transpiration_mm_per_h:.2f} mm/h")

# --- 過去の実績と並べる（★同じくらいの明るさの日どうしで比べる）---
if history is not None:
    if history.matched_brightness:
        headline = f"同じ明るさの日（{history.median_radiation_mj:.1f} MJ 前後）"
        note = (
            f"入れた日射 {radiation:.1f} MJ/m² に近い日だけを選んである。"
            f"前後{history.window_days}日・{history.days}日ぶん。"
        )
    else:
        headline = "同じ時期の全部の天気"
        note = (
            f"日射 {radiation:.1f} MJ/m² に近い日が足りなかったので、"
            f"天気を問わない値を出している。明るさがそろっていないので"
            f"読むときは気をつけること。"
        )

    compare = st.columns(4)
    compare[0].metric(
        "実績・ふだん", f"{history.median_l_per_m2:.2f} L/m²",
        delta=f"{advice.recommended_irrigation_l_per_m2 - history.median_l_per_m2:+.2f} "
              f"L/m² との差",
        delta_color="off",
        help=f"{headline}の中央値。{note}",
    )
    compare[1].metric(
        "実績・最大", f"{history.max_l_per_m2:.2f} L/m²",
        help=f"{headline}でいちばん多く出した日。",
    )
    compare[2].metric("実績・平均", f"{history.mean_l_per_m2:.2f} L/m²")
    compare[3].metric(
        "集計した日数", f"{history.days} 日",
        help=f"{history.season_range} 年産の記録から。前後{history.window_days}日。",
    )
    ratio = (
        advice.recommended_irrigation_l_per_m2 / history.median_l_per_m2
        if history.median_l_per_m2 > 0 else 0.0
    )
    st.caption(
        f"**比べたのは「{history.basis}」。**"
        f"モデルの合計潅水量は、そのふだんの量の **{ratio:.2f} 倍**。"
        f"　出どころ: {describe_source()}。"
        f"**この実績はモデルの計算にはいっさい使っていない**（答え合わせ用）。"
    )
    if not history.matched_brightness:
        st.info(
            "同じ明るさの日が足りず、天気を問わない値で比べている。"
            "晴れの日の計算値を曇りこみの平均と並べることになるので、"
            "モデルが大きめに見えるのはふつう。"
        )
else:
    st.caption(
        f"この時期（{target_date.month}月{target_date.day}日前後）は"
        f"潅水の記録が足りないので、過去との比較は出せない。"
    )

columns = st.columns(4)
columns[0].metric("地下へ抜ける量", f"{advice.expected_drainage_l_per_m2:.2f} L/m²")
columns[1].metric("流亡率", f"{advice.drainage_fraction * 100:.0f}%")
columns[2].metric("夕方の土壌 pF", f"{advice.end_pf:.2f}")
columns[3].metric(
    "最小空気率",
    f"{advice.min_air_filled_porosity * 100:.1f}%",
    delta=f"{(advice.min_air_filled_porosity - 0.10) * 100:+.1f} ポイント",
    delta_color="normal",
    help="根の呼吸には 10% 以上ほしい。",
)

for warning in advice.warnings:
    st.warning(warning)

if month in MONTHS_WITH_WEAK_DATA:
    st.warning(
        f"**{month}月は「平年なみ」が当てにならない。**"
        f"平年の曲線は 2025年9月〜2026年6月の実測に当てはめたもので、"
        f"7・8月は栽培しておらず外挿になっている。"
        f"あとから入ったデータと比べると、この時期は平年が実測より"
        f"2〜4割低く出る（7月 0.73倍・8月 0.62倍・9月 0.77倍）。"
        f"「快晴」からの割合で見るか、日射・気温・湿度を手で入れたほうがよい。"
    )

_day = solar_day_from_date(target_date)
_share = radiation / clear_mj if clear_mj > 0 else 0.0
st.caption(
    f"平年（{target_date.month}月{target_date.day}日）のハウス内日射は "
    f"{normal_mj:.1f} MJ/m²。今日はそれより {radiation - normal_mj:+.1f} MJ/m² "
    f"{'明るい' if radiation >= normal_mj else '暗い'}見込み。"
    f"入れた日射は快晴の **{_share * 100:.0f}%**。　／　"
    f"可照時間 **{_day.daylength_h:.1f} 時間**"
    f"（日の出 {_day.sunrise_h:.2f} ／ 南中 {_day.solar_noon_h:.2f} ／ "
    f"日の入り {_day.sunset_h:.2f} 時）"
)


# =============================================================================
# 5. 潅水基準を変えたらどうなるか
# =============================================================================

st.divider()
st.subheader("潅水基準を変えたらどうなるか")
st.caption(
    "実務で使う **L/10MJ**（実質日射10MJあたりの潅水量）で振ったときに、"
    "流亡・pF・空気率がどう落ち着くか。"
    "**同じ天気が毎日つづいたときの落ち着き先**を出している。"
)

rows = []
unsettled = []
for basis in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
    irrigation = basis * advice.effective_radiation_mj / 10.0
    steady = steady_state(
        irrigation_mm=irrigation,
        potential_transpiration_mm=advice.transpiration_l_per_m2,
        settings=water_settings,
        initial_water_content=water_content_from_potential(potential_from_pf(start_pf)),
    )
    balance = steady.balance
    if not steady.settled:
        unsettled.append(basis)
    rows.append({
        "潅水基準 [L/10MJ]": basis,
        "1日の潅水 [L/m²]": round(irrigation, 2),
        "蒸散 [L/m²]": round(balance.actual_transpiration_mm, 2),
        "潅水/蒸散": round(
            irrigation / balance.actual_transpiration_mm, 2
        ) if balance.actual_transpiration_mm > 0 else None,
        "流亡率 [%]": round(balance.drainage_fraction * 100),
        "夕方 pF": round(balance.end_pf, 2),
        "最小空気率 [%]": round(balance.min_air_filled_porosity * 100, 1),
        # 落ち着かない行は「—」を入れるので、列全体を文字列にそろえる
        # （数値と文字列が混ざると表の変換で警告が出る）
        "落ち着くまで [日]": f"{steady.days}" if steady.settled else "—",
    })
st.dataframe(rows, width="stretch", hide_index=True)

if unsettled:
    st.warning(
        f"**{'・'.join(f'{b:.1f}' for b in unsettled)} L/10MJ は落ち着かない。**"
        f"潅水が蒸散に足りず、土が乾きつづける。"
        f"表の pF はその時点の値で、日を追うごとにさらに上がる（＝もっと乾く）。"
    )

# 落ち着かない行はひとつ上で別に知らせているので、ここでは除く
dry_rows = [
    r for r in rows
    if r["夕方 pF"] >= 2.7 and r["潅水基準 [L/10MJ]"] not in unsettled
]
if dry_rows:
    st.warning(
        f"**{'・'.join(f'{r['潅水基準 [L/10MJ]']:.1f}' for r in dry_rows)} L/10MJ は"
        f"乾かしすぎ。**"
        f"落ち着いてはいるが、それは土が乾いて根が吸えなくなり、"
        f"蒸散が潅水量まで抑えられた結果（教科書 式9.22）。"
        f"**落ち着いた＝健全、ではない。**pF と蒸散量を必ず見ること。"
    )

st.caption(
    f"いまの見立て（上乗せ {leaching_percent}%）は "
    f"**{advice.water_per_10mj_l_per_m2:.2f} L/10MJ** にあたる。"
    f"デルフィのコーチ値は 1.5〜3.0 L/10MJ。"
    f"蒸散量が潅水を増やしても頭打ちになるのは、圃場容水量を超えた水が"
    f"吸水能力を増やさないため（式9.20 で有効水分度が1で頭打ちになる）。"
    f"余った水は地下へ抜けるだけ。"
)
st.caption(
    "**空気率は鵜呑みにしないこと。**目安の10%を下回る値が出るが、"
    "実際に栽培できている以上、実測していない前提"
    "（飽和透水係数・根群域の深さ・点滴で濡れる面積の割合）の"
    "どれかが実態と合っていない可能性がある。"
)


# =============================================================================
# 6. 光合成（参考）
# =============================================================================

st.divider()
st.subheader("その日の光合成（参考）")
st.caption(
    "★蒸散より信頼度が低い。実収量から逆算すると約3倍足りない。"
    "日ごとの上下や季節の傾向を見るぶんには使えるが、絶対量は当てにしない。"
)

columns = st.columns(3)
columns[0].metric("CO₂ 固定量", f"{advice.co2_fixed_g_per_m2:.1f} g/m²")
columns[1].metric("糖の生産量", f"{advice.sugar_g_per_m2:.1f} g/m²")
columns[2].metric(
    "10アールあたりの糖",
    f"{advice.sugar_g_per_m2 * M2_PER_10A / 1000:.1f} kg",
)


# =============================================================================
# 7. 参考（たたむ）
# =============================================================================

st.divider()

with st.expander("なぜ光合成は当てにならないのか"):
    st.write(
        "**光ではなく CO₂ の通り道が頭打ちを作っている。**\n\n"
        "光応答の式は直角双曲線から非直角双曲線に直したので、"
        "光飽和の形はまともになった（PPFD 2000 での飽和度が 0.74 → 0.91）。\n\n"
        "ただし飽和光・25℃・400ppm での葉の純光合成は、"
        "Pmax をいくら上げても 26.7 μmol/(m²·s) で止まる。"
        "CO₂ の総コンダクタンス gc = 0.164 mol/(m²·s) が律速している。\n\n"
        "比較先の文献値（20〜28）はキュベット内の測定で、強制送風により"
        "境界層コンダクタンスが野外の5〜10倍ある。土俵が違う。\n\n"
        "直すには Michaelis 型の KC を Farquhar 型"
        "（Rubisco律速とRuBP再生律速の小さいほう、葉肉コンダクタンスを含む）に"
        "書き換える必要がある。**蒸散側には一切影響しない。**"
    )

with st.expander("この土のこと（シルト質壌土・実測のpF曲線から）"):
    fc = describe_soil_state(FIELD_CAPACITY_WATER_CONTENT)
    pwp_theta = WILTING_POINT_WATER_CONTENT
    st.write(
        f"- 圃場容水量 pF **{pf_from_potential(fc.potential_j_per_kg):.2f}**"
        f"（θ = {FIELD_CAPACITY_WATER_CONTENT:.3f}、空気率 "
        f"{fc.air_filled_porosity * 100:.1f}%）\n"
        f"- 永久しおれ点 pF **4.18**（θ = {pwp_theta:.3f}）\n"
        f"- 有効水分量 **{FIELD_CAPACITY_WATER_CONTENT - pwp_theta:.3f}**"
        f"　→ 根群域 {root_depth * 100:.0f} cm で "
        f"**{(FIELD_CAPACITY_WATER_CONTENT - pwp_theta) * root_depth * 1000:.0f} mm**\n\n"
        f"この土は圃場容水量でも空気率が {fc.air_filled_porosity * 100:.1f}% しかなく、"
        f"根の目安（10%以上）をわずかに下回る。水をやりすぎる余地は小さい。"
    )

with st.expander("快晴・曇天の日射はどう計算しているか"):
    st.write(
        f"太陽の位置は暦と緯度経度だけで決まるので、"
        f"センサーがなくても「その日、快晴ならどれだけ日射があるか」を"
        f"計算できる。式は日照予測アプリ（{SOLAR_SITE['LATITUDE_DEG']}°N, "
        f"{SOLAR_SITE['LONGITUDE_DEG']}°E, 標高 {SOLAR_SITE['ALTITUDE_M']:.0f} m）"
        f"と同じもので、生物環境物理学の基礎 第11章の式11.1〜11.13。\n\n"
        f"**大気透過率 τ** は、大気を垂直に1回通り抜けたときに残る割合。"
        f"既定値は日付から決まる快晴値（7月 0.57 ／ 1月 0.87）で、"
        f"夏に低いのは水蒸気が多いため。黄砂や花粉の日は 0.05〜0.15 下がる。"
        f"「日射の細かい設定」タブから動かせる。"
    )
    st.write(
        f"τ の式は**宇都宮気象台の実測から当てはめた**もの。"
        f"作業日誌の「日射量(宇都宮)」列 2251 日から晴れた上位15%（341日）を選び、"
        f"その日の実測を再現する τ を1日ずつ解いて調和関数を当てた。\n\n"
        f"　τ = 0.7185 + 0.1429·cos(2πJ/365) + 0.0400·sin(2πJ/365)\n\n"
        f"解けた τ は 0.57〜0.87 で、すべて物理的にありうる範囲に収まった。"
        f"日ごとのばらつきは標準偏差で 0.04 ほどある。"
    )
    st.write(
        f"**ハウス内センサーに直すときは、その日の透過率を掛ける。**"
        f"透過率も季節で動く（実測: 5月中旬 0.72 ／ 11月中旬 0.64）。"
        f"冬は太陽高度が低く、フィルムへの入射角が浅くなるぶん反射で失われる。"
        f"快晴の日だけで見ても同じ形が出るので、天気のせいではなく光学的なもの。\n\n"
        f"この日は **τ {tau:.3f} × 透過率 {transmittance_today:.3f}** で計算している。"
    )

with st.expander("「気象予報の全天日射量」と「実質日射」は何が違うのか"):
    st.write(
        f"**仕事が違う。片方は物理、片方は単位の取り決め。**\n\n"
        f"| | 計算 | この日の値 |\n|---|---|---|\n"
        f"| 気象予報の全天日射量 | センサー値 ÷ その日の透過率 "
        f"({transmittance_today:.3f}) | {forecast_outside_mj:.1f} MJ/m² |\n"
        f"| 実質日射 | センサー値 ÷ レシピ係数 ({recipe_coef:.2f}) × フィルム "
        f"| {advice.effective_radiation_mj:.1f} MJ/m² |\n\n"
        f"上は**本当のハウス外の日射**。気象予報と比べるならこちら。\n\n"
        f"下は**液肥混入機レシピと数字を合わせるためだけの値**で、"
        f"物理的に正しい必要がない。レシピが出す水の量は\n\n"
        f"　潅水量 ＝ 10MJあたり潅水量 × 面積 × 実質日射 ÷ 10\n\n"
        f"で、こちらも同じ係数 {recipe_coef:.2f} で割って"
        f"「10MJあたり潅水量」を作っているので、**分子と分母で打ち消し合う**。"
        f"だから係数が実測とずれていても、出る水は狙ったとおりになる。"
    )
    st.write(
        f"この日の2つの比は {advice.effective_radiation_mj / forecast_outside_mj:.2f} 倍。"
        f"実測の透過率は季節で 0.62〜0.72 と動くのに対し、"
        f"レシピ係数は年間固定なので、春は実質日射が本物より大きめ、"
        f"秋は小さめに出る。\n\n"
        f"**デルフィ指標（1.5〜3 L/10MJ・ハウス外基準）と比べるときだけ注意。**"
        f"春の 3.0 L/10MJ は、本当の外の日射に対しては 3.3 相当にあたる。"
    )

with st.expander("以前あった「較正係数 1.20」をなぜやめたのか"):
    st.write(
        "**一律の係数では直らない誤りだったから。**"
        "快晴の計算値が実測より低く出るのを 1.20 倍で埋めていたが、"
        "気象台の実測（快晴日341日）と比べると、"
        "ずれの向きが季節で逆だった。\n\n"
        "| 月 | 1月 | 3月 | 5月 | 9月 | 10月 |\n|---|---|---|---|---|---|\n"
        "| τ そのまま | 0.905 | 0.934 | 0.933 | 1.084 | 0.970 |\n"
        "| ×1.20 | 1.086 | 1.121 | 1.120 | **1.300** | 1.165 |\n\n"
        "冬は1割足りないが、8〜9月はすでに足りている。"
        "そこへ一律1.20を掛けたので、秋が3割も過大になっていた。"
    )
    st.write(
        "**1.20 の正体は、2つの誤りが同じ向きに重なったもの。**\n\n"
        "| | 倍率 |\n|---|---|\n"
        "| τ の振幅が足りない（0.100 → 実測 0.148） | 1.077 |\n"
        "| 透過率が実測より低い（0.65 → 実測 0.678） | 1.046 |\n"
        "| **積** | **1.127** |\n\n"
        "実際に使っていた 1.20 は、これより 6% 行きすぎていた。"
        "原因のほうを直したので係数は要らなくなり、"
        "月ごとのばらつきは幅 0.21 → **0.086** と 2.5 倍そろった。"
    )
    st.write(
        "**教訓: 合わない量を一律の係数で埋めようとしたら、まず疑う。**"
        "ずれが季節や条件で向きを変えるなら原因は複数あり、"
        "一律の係数はそれを隠すだけで直さない。"
    )

with st.expander("この見積りはどれくらい当たるか"):
    st.write(
        "実測242日（中央ハウス 2025年11月〜2026年6月）で、"
        "5分ごとの実測値を流した結果と突き合わせた成績。\n\n"
        "- 比の中央値 **0.991**\n"
        "- 相関係数 **0.984**\n"
        "- ±20%以内に収まる日 **78%**、±30%以内 **91%**\n"
        "- 平均絶対誤差 **0.22 L/m²/日**\n\n"
        "日射が 4 MJ/m² を下回る暗い日は当たりが悪く、1.2倍ほど多めに出る。"
    )
    st.write(
        "**いちばん効く前提は LAI。**"
        "葉1枚の面積 0.12 m² は暫定値なので、"
        "代表的な葉を5〜10枚測れば幅がぐっと縮む"
        "（複葉の面積 ≒ 葉長 × 最大幅 × 0.5〜0.6 で概算できる）。"
    )
    st.write(
        "**過去の潅水記録はモデルには使っていない。**"
        "蒸散量は物理式（エネルギー収支・気孔コンダクタンス・Beer則・"
        "Campbell 式9.20/9.22）だけで出している。"
        "実績は画面に並べて答え合わせに使うだけ。"
    )
