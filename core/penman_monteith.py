"""
Penman-Monteith 式による群落蒸散の計算（教科書 第14章 式14.12）。

【なぜ2つ目の計算方法を作るのか】
core/canopy.py は群落を20層に分け、層ごとに葉のエネルギー収支を解いている。
一方こちらは、群落全体を「1枚の大きな葉」とみなして解析的に解く。

    多層モデル      ： 層ごとに葉温を数値的に探す
    Penman-Monteith ： 葉温を式の上で消去してあり、一発で解ける

同じ現象を別の道筋で計算しているので、両者が一致すれば、
どちらの実装も正しいという強い裏づけになる。食い違えば、
どちらかに誤りがあるということなので、原因を追える。

【式14.12】

                s(R_abs − εs σ Ta⁴ − G) + γ* λ gv D/pa
    λE_canopy = ────────────────────────────────────────
                              s + γ*

  λE_canopy : 潜熱フラックス（蒸散のエネルギー）[W/m²]
  s         : 飽和水蒸気圧曲線の傾きを大気圧で割ったもの [1/K]（式3.10）
  R_abs     : 群落が吸収した放射 [W/m²]
  εs σ Ta⁴  : 群落から出ていく長波放射 [W/m²]
  G         : 土壌熱フラックス [W/m²]
  γ*        : 見かけの乾湿計定数 [1/K]
  λ         : 水のモル気化潜熱 [J/mol]
  gv        : 群落の水蒸気コンダクタンス [mol/(m²·s)]
  D         : 飽差 [kPa]
  pa        : 大気圧 [kPa]

【式の読み方】
分子の第1項は「放射で与えられたエネルギーのうち蒸散に回る分」、
第2項は「乾いた空気が水を引き出す分」を表す。
前者を放射項、後者を移流項と呼ぶ。
ハウス内は風が弱いので放射項が支配的になる。

【γ* とは】
見かけの乾湿計定数。本来の乾湿計定数 γ に、熱と水蒸気の
通りやすさの比を掛けたもの。

    γ* = γ · gHa / gv

教科書P251の基準蒸発散の例（gv = 0.6·0.2u/(0.6+0.2u)、gHa = 0.2u）に
この定義を当てはめると γ* = 6.67e-4(1 + u/3) となり、本文の式と一致する。

【前提条件】
1. 群落を1枚の葉とみなすので、群落内の日射・温度・湿度の勾配は
   平均化される。多層モデルより粗い扱いになる。
2. 土壌熱フラックス G は、白マルチ被覆率100%のため小さい。
   マルチが日射を反射し、地面への熱の流入を抑えるため。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from config import (
    LEAF_EMISSIVITY,
    LEAF_SIDES_FOR_HEAT,
    MOLAR_LATENT_HEAT_WATER_J_PER_MOL,
    MOLAR_MASS_WATER_G_PER_MOL,
    PSYCHROMETRIC_GAMMA_PER_K,
    STEFAN_BOLTZMANN,
    TETENS_A_KPA,
    TETENS_B,
    TETENS_C_C,
)

# 土壌熱フラックスを、群落を通り抜けた日射の何割とみなすか。
# 白マルチが敷いてあり被覆率100%なので、地面に届いた日射の多くは
# 反射され、土壌に入る熱は小さい。裸地なら0.3程度だが、
# マルチありでは0.05〜0.1が目安。
SOIL_HEAT_FLUX_FRACTION = 0.07


@dataclass(frozen=True)
class PenmanMonteithResult:
    """Penman-Monteith 式による計算結果。"""

    latent_heat_w_per_m2: float        # 潜熱フラックス λE [W/m²]
    transpiration_mg_per_m2_s: float   # 蒸散速度 [mg/(m²·s)]
    transpiration_mol_per_m2_s: float  # 蒸散速度 [mol/(m²·s)]

    # 内訳（どちらが効いているかを見るため）
    radiation_term_w_per_m2: float     # 放射項
    advection_term_w_per_m2: float     # 移流項
    net_radiation_w_per_m2: float      # 正味放射 R_abs − εσTa⁴ − G
    soil_heat_flux_w_per_m2: float     # 土壌熱フラックス G
    slope_s: float                     # 飽和水蒸気圧曲線の傾き s
    gamma_star: float                  # 見かけの乾湿計定数 γ*

    @property
    def radiation_fraction(self) -> float:
        """蒸散のうち放射項が占める割合。"""
        total = self.radiation_term_w_per_m2 + self.advection_term_w_per_m2
        if total == 0:
            return 0.0
        return self.radiation_term_w_per_m2 / total


def saturation_slope_per_k(temp_c: float, pressure_kpa: float) -> float:
    """
    飽和水蒸気圧曲線の傾き s を求める（教科書 式3.10）。

        s = (dCs/dT) = Δ / pa

    Cs は飽和水蒸気のモル分率なので、s の単位は 1/K になる。

    【物理的な意味】
    気温が1℃上がったとき、空気が抱えられる水蒸気がどれだけ増えるか。
    気温が高いほど大きくなる（指数関数の微分なので）。
    この値が大きいほど、同じ放射エネルギーがより多く蒸散に回る。

    Args:
        temp_c: 気温 [℃]
        pressure_kpa: 大気圧 [kPa]

    Returns:
        傾き s [1/K]
    """
    denominator = temp_c + TETENS_C_C
    if denominator <= 0:
        raise ValueError(
            f"飽和水蒸気圧曲線の傾きを計算できません: 気温 {temp_c} ℃ は "
            f"Tetens式の適用範囲外です"
        )

    # es(T) = a·exp(bT/(T+c)) を T で微分すると
    #   des/dT = es(T) · b·c / (T+c)²
    es = TETENS_A_KPA * math.exp(TETENS_B * temp_c / denominator)
    des_dt = es * TETENS_B * TETENS_C_C / denominator ** 2

    return des_dt / pressure_kpa


def apparent_psychrometric_constant(
    heat_conductance: float,
    vapor_conductance: float,
) -> float:
    """
    見かけの乾湿計定数 γ* を求める。

        γ* = γ · gHa / gv

    【物理的な意味】
    熱の逃げやすさと水蒸気の逃げやすさの比。
    気孔が閉じて gv が小さくなると γ* が大きくなり、
    式14.12 の分母が大きくなって蒸散が減る。

    Args:
        heat_conductance: 熱のコンダクタンス gHa [mol/(m²·s)]
        vapor_conductance: 水蒸気のコンダクタンス gv [mol/(m²·s)]

    Returns:
        見かけの乾湿計定数 γ* [1/K]
    """
    if vapor_conductance <= 0:
        # 気孔が完全に閉じている。蒸散は起きないので、
        # γ* を無限大として扱えば式14.12 の結果は 0 に近づく。
        return float("inf")
    if heat_conductance < 0:
        raise ValueError(f"熱コンダクタンスが負です: {heat_conductance}")

    return PSYCHROMETRIC_GAMMA_PER_K * heat_conductance / vapor_conductance


def calculate_penman_monteith(
    air_temp_c: float,
    vapor_pressure_kpa: float,
    absorbed_radiation_w_per_m2: float,
    below_canopy_solar_w_per_m2: float,
    vapor_conductance: float,
    heat_conductance: float,
    pressure_kpa: float,
    soil_heat_flux_fraction: float = SOIL_HEAT_FLUX_FRACTION,
) -> PenmanMonteithResult:
    """
    Penman-Monteith 式（式14.12）で群落の蒸散を求める。

    Args:
        air_temp_c: 気温 Ta [℃]
        vapor_pressure_kpa: 大気の水蒸気圧 ea [kPa]
        absorbed_radiation_w_per_m2: 群落が吸収した日射 R_abs [W/m²]
        below_canopy_solar_w_per_m2: 群落を通り抜けた日射 [W/m²]
            土壌熱フラックス G の推定に使う
        vapor_conductance: 群落の水蒸気コンダクタンス gv [mol/(m²·s)]
        heat_conductance: 群落の熱コンダクタンス gHa [mol/(m²·s)]
        pressure_kpa: 大気圧 pa [kPa]
        soil_heat_flux_fraction: 群落を抜けた日射のうち土壌に入る割合

    Returns:
        PenmanMonteithResult
    """
    air_temp_k = air_temp_c + 273.15

    # --- 飽差 D ---------------------------------------------------------
    es = TETENS_A_KPA * math.exp(TETENS_B * air_temp_c / (air_temp_c + TETENS_C_C))
    vapor_deficit_kpa = max(es - vapor_pressure_kpa, 0.0)

    # --- 傾き s ---------------------------------------------------------
    s = saturation_slope_per_k(air_temp_c, pressure_kpa)

    # --- 見かけの乾湿計定数 γ* -------------------------------------------
    gamma_star = apparent_psychrometric_constant(heat_conductance, vapor_conductance)

    if math.isinf(gamma_star):
        # 気孔が完全に閉じている
        return PenmanMonteithResult(
            latent_heat_w_per_m2=0.0,
            transpiration_mg_per_m2_s=0.0,
            transpiration_mol_per_m2_s=0.0,
            radiation_term_w_per_m2=0.0,
            advection_term_w_per_m2=0.0,
            net_radiation_w_per_m2=0.0,
            soil_heat_flux_w_per_m2=0.0,
            slope_s=s,
            gamma_star=gamma_star,
        )

    # --- 放射の収支 ------------------------------------------------------
    # 群落から出ていく長波放射。群落を1枚の葉とみなすので、
    # 気温と同じ温度の面から上下に放射されるとして扱う。
    outgoing_longwave = (
        LEAF_SIDES_FOR_HEAT * LEAF_EMISSIVITY * STEFAN_BOLTZMANN * air_temp_k ** 4
    )
    # 周囲（被覆・地面）から受ける長波放射も同じだけあるので、
    # 気温と群落温度が等しければ正味はゼロになる。
    # 式14.12 では R_abs に長波の入射が含まれている前提なので、
    # ここでは吸収した日射に長波の入射を加える。
    incoming_longwave = outgoing_longwave

    # --- 土壌熱フラックス G ----------------------------------------------
    # 白マルチ被覆率100%なので、地面に届いた日射の多くは反射される。
    soil_heat_flux = below_canopy_solar_w_per_m2 * soil_heat_flux_fraction

    # 正味放射
    net_radiation = (
        absorbed_radiation_w_per_m2 + incoming_longwave
        - outgoing_longwave - soil_heat_flux
    )

    # --- 式14.12 --------------------------------------------------------
    # 放射項: s(R_abs − εσTa⁴ − G)
    radiation_term = s * net_radiation

    # 移流項: γ* λ gv D/pa
    advection_term = (
        gamma_star
        * MOLAR_LATENT_HEAT_WATER_J_PER_MOL
        * vapor_conductance
        * vapor_deficit_kpa
        / pressure_kpa
    )

    latent_heat = (radiation_term + advection_term) / (s + gamma_star)

    # 潜熱 [W/m²] を蒸散速度に直す
    transpiration_mol = latent_heat / MOLAR_LATENT_HEAT_WATER_J_PER_MOL
    transpiration_mg = transpiration_mol * MOLAR_MASS_WATER_G_PER_MOL * 1000.0

    return PenmanMonteithResult(
        latent_heat_w_per_m2=latent_heat,
        transpiration_mg_per_m2_s=transpiration_mg,
        transpiration_mol_per_m2_s=transpiration_mol,
        radiation_term_w_per_m2=radiation_term / (s + gamma_star),
        advection_term_w_per_m2=advection_term / (s + gamma_star),
        net_radiation_w_per_m2=net_radiation,
        soil_heat_flux_w_per_m2=soil_heat_flux,
        slope_s=s,
        gamma_star=gamma_star,
    )


def calculate_from_canopy(
    canopy,
    air_temp_c: float,
    vapor_pressure_kpa: float,
    pressure_kpa: float,
    soil_heat_flux_fraction: float = SOIL_HEAT_FLUX_FRACTION,
) -> PenmanMonteithResult:
    """
    多層モデルの結果（CanopyGasExchange）から Penman-Monteith 式を計算する。

    多層モデルが求めた群落コンダクタンスと吸収放射をそのまま使うので、
    2つの方法の差は「葉温をどう扱うか」だけになる。
    この条件で結果が一致すれば、両方の実装が正しいと確かめられる。

    Args:
        canopy: canopy.calculate_canopy_gas_exchange() の戻り値
        air_temp_c: 気温 [℃]
        vapor_pressure_kpa: 水蒸気圧 [kPa]
        pressure_kpa: 大気圧 [kPa]
        soil_heat_flux_fraction: 群落を抜けた日射のうち土壌に入る割合

    Returns:
        PenmanMonteithResult
    """
    return calculate_penman_monteith(
        air_temp_c=air_temp_c,
        vapor_pressure_kpa=vapor_pressure_kpa,
        absorbed_radiation_w_per_m2=canopy.absorbed_radiation_w_per_m2,
        below_canopy_solar_w_per_m2=canopy.below_canopy_solar_w_per_m2,
        vapor_conductance=canopy.canopy_conductance_water,
        heat_conductance=canopy.canopy_conductance_heat,
        pressure_kpa=pressure_kpa,
        soil_heat_flux_fraction=soil_heat_flux_fraction,
    )
