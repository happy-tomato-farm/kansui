"""過去に、同じ時期・同じくらいの明るさの日にどれだけ潅水したかを引く。

【何のためにあるか】
モデルが出した潅水量が妥当かを、自分の過去の実績と並べて確かめるため。

★モデルの計算にはいっさい使わない。答え合わせのための参考値。

【★明るさをそろえて比べること】
以前は「その時期の全部の日」の中央値と比べていたが、これは比べ方が悪かった。
アプリの既定は「快晴」なので、晴れた日の計算値を、曇りも含めた平均と
並べていたことになる。10月なら快晴 10.7 MJ に対し平年なみ 6.2 MJ と
1.7倍も違うので、モデルが過大に見えて当然だった。

そこで、引くときに「その日の日射」を渡してもらい、
**同じくらいの明るさだった日だけ**を選んで返すようにしてある。

【データの範囲】
2023年8月以降の3作期のみ（中央 759 日・東 764 日）。
それより前は記録の取り方も潅水のやり方も違うので使わない。

【作り直し方】
作業日誌が更新されたら、このフォルダで

    py tools/潅水実績を作り直す.py
"""

from __future__ import annotations

import datetime as dt
import json
import statistics
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

#: JSON の置き場所。このファイルから見た相対位置。
HISTORY_PATH = Path(__file__).resolve().parent.parent / "data" / "潅水実績.json"

#: 日付の窓を広げていく順番 [日]。足りなければ次へ。
WINDOW_STEPS = (10, 15, 21, 30)

#: 明るさをそろえるときの許容幅。割合と最小幅の大きいほう。
#: 割合だけだと暗い日に窓が狭くなりすぎるので、最小幅を持たせる。
BRIGHTNESS_TOLERANCE_RATIO = 0.25
BRIGHTNESS_TOLERANCE_MIN_MJ = 1.5

#: これだけ集まれば窓を広げるのをやめる。
ENOUGH_DAYS = 12

#: これを下回ったら、統計として意味がないので返さない。
MINIMUM_DAYS = 3

#: 1年の日数。通日の距離を求めるのに使う（うるう年を含めて366）。
DAYS_IN_YEAR = 366


@dataclass(frozen=True)
class IrrigationHistory:
    """その条件で実際どれだけ潅水していたか [L/m²/日]。"""

    house: str
    days: int                   #: 集計に使った日数
    median_l_per_m2: float      #: 中央値（ふだんの量）
    max_l_per_m2: float         #: 同じ条件で出した最大
    mean_l_per_m2: float        #: 平均
    median_radiation_mj: float  #: 選ばれた日の日射の中央値
    window_days: int            #: 前後何日まで広げたか
    matched_brightness: bool    #: 明るさをそろえられたか
    seasons: tuple[int, ...]    #: 集計に含まれる作期（8月始まりの年）

    @property
    def season_range(self) -> str:
        """作期の範囲を「2023〜2025」の形で返す。"""
        if not self.seasons:
            return "—"
        if len(self.seasons) == 1:
            return str(self.seasons[0])
        return f"{min(self.seasons)}〜{max(self.seasons)}"

    @property
    def basis(self) -> str:
        """何と比べた値なのかを一言で返す（画面の注記用）。"""
        if self.matched_brightness:
            return (
                f"前後{self.window_days}日・日射 {self.median_radiation_mj:.1f} MJ 前後の"
                f"{self.days}日"
            )
        return f"前後{self.window_days}日の全部の天気・{self.days}日"


@lru_cache(maxsize=1)
def _load() -> dict:
    """JSON を読む。何度呼ばれても1回しか読まない（lru_cache）。"""
    if not HISTORY_PATH.exists():
        raise FileNotFoundError(
            f"潅水実績のファイルが見つからない: {HISTORY_PATH}\n"
            f"`py tools/潅水実績を作り直す.py` を実行して作ること。"
        )
    try:
        return json.loads(HISTORY_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"潅水実績のファイルが壊れている: {HISTORY_PATH}（{error}）"
        ) from error


def _records(house: str) -> list[list[float]]:
    """そのハウスの生データ [[通日, 年, 日射MJ, 潅水L], ...] を返す。"""
    houses = _load().get("ハウス", {})
    if house not in houses:
        raise KeyError(
            f"潅水実績に '{house}' の記録がない。あるのは {list(houses.keys())}。"
        )
    return houses[house]


def _day_distance(day_a: int, day_b: int) -> int:
    """通日どうしの距離 [日]。年をまたぐ側も見る（12月31日と1月1日は1日）。"""
    gap = abs(day_a - day_b)
    return min(gap, DAYS_IN_YEAR - gap)


def _season_of(day_of_year: int, year: int) -> int:
    """その日が属する作期（8月始まり）を返す。"""
    # 8月1日は平年で通日213、うるう年で214。境目は214で足りる。
    return year if day_of_year >= 213 else year - 1


def _summarise(
    house: str,
    selected: list[list[float]],
    window_days: int,
    matched_brightness: bool,
) -> IrrigationHistory:
    """選ばれた日から統計を作る。"""
    waters = [record[3] for record in selected]
    radiations = [record[2] for record in selected]
    seasons = {_season_of(int(record[0]), int(record[1])) for record in selected}
    return IrrigationHistory(
        house=house,
        days=len(selected),
        median_l_per_m2=statistics.median(waters),
        max_l_per_m2=max(waters),
        mean_l_per_m2=statistics.fmean(waters),
        median_radiation_mj=statistics.median(radiations),
        window_days=window_days,
        matched_brightness=matched_brightness,
        seasons=tuple(sorted(seasons)),
    )


def lookup(
    date: dt.date,
    house: str = "中央",
    sensor_radiation_mj: float | None = None,
) -> IrrigationHistory | None:
    """その時期・その明るさの日に、実際どれだけ潅水していたかを返す。

    Args:
        date: 知りたい日付（年は無視して月日だけ使う）
        house: "中央" または "東"
        sensor_radiation_mj: その日のハウス内日射 [MJ/m²]。
            渡すと、**同じくらいの明るさだった日だけ**を選ぶ。
            省くと天気を問わず、その時期の全部の日を使う。

    Returns:
        IrrigationHistory。記録が足りない時期（作期外など）は None。
        明るさをそろえると日が足りない場合は、窓を広げ、
        それでも足りなければ天気を問わない値に落として返す
        （`matched_brightness` が False になる）。
    """
    if sensor_radiation_mj is not None and sensor_radiation_mj < 0.0:
        raise ValueError(
            f"日射は 0 以上で指定する。渡された値: {sensor_radiation_mj}"
        )

    records = _records(house)
    target_day = date.timetuple().tm_yday

    # ★まず「そもそも作期の中か」を見る。
    # いちばん狭い窓に記録がないなら作期外なので、窓を広げてはいけない。
    # 広げると、8月（栽培していない）に9月の記録を持ってきてしまう。
    in_season = [
        record for record in records
        if _day_distance(int(record[0]), target_day) <= WINDOW_STEPS[0]
    ]
    if len(in_season) < MINIMUM_DAYS:
        return None

    if sensor_radiation_mj is not None:
        tolerance = max(
            sensor_radiation_mj * BRIGHTNESS_TOLERANCE_RATIO,
            BRIGHTNESS_TOLERANCE_MIN_MJ,
        )
        # 窓を広げながら、十分に集まったところで止める。
        best: list[list[float]] = []
        best_window = WINDOW_STEPS[0]
        for window in WINDOW_STEPS:
            selected = [
                record for record in records
                if _day_distance(int(record[0]), target_day) <= window
                and abs(record[2] - sensor_radiation_mj) <= tolerance
            ]
            best, best_window = selected, window
            if len(selected) >= ENOUGH_DAYS:
                break
        if len(best) >= MINIMUM_DAYS:
            return _summarise(house, best, best_window, matched_brightness=True)

    # 明るさをそろえられなかった（か、指定がなかった）ので天気は問わない。
    selected, used_window = in_season, WINDOW_STEPS[0]
    for window in WINDOW_STEPS:
        if len(selected) >= ENOUGH_DAYS:
            break
        selected = [
            record for record in records
            if _day_distance(int(record[0]), target_day) <= window
        ]
        used_window = window
    return _summarise(house, selected, used_window, matched_brightness=False)


def describe_source() -> str:
    """データの出どころを一文で返す（画面の注記用）。"""
    payload = _load()
    seasons = payload.get("作期", [])
    span = (f"{min(seasons)}〜{max(seasons)}年産" if seasons else "—")
    return f"{span}（{payload.get('期間', '—')}）の作業日誌より"
