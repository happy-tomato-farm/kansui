@echo off
chcp 65001 > nul
cd /d "%~dp0"
echo.
echo   今日の潅水量のめやす を起動します。
echo   ブラウザが自動で開きます。閉じるときはこの黒い画面で Ctrl+C。
echo.
py -m streamlit run app.py
pause
