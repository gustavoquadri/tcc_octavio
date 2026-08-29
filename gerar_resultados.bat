@echo off
chcp 65001 >nul
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
  echo Ambiente .venv nao encontrado.
  echo Crie com:  py -3.12 -m venv .venv
  echo Depois:    .\.venv\Scripts\python.exe -m pip install -r requirements.txt
  pause
  exit /b 1
)

echo Gerando resultados do Capitulo 5...
".venv\Scripts\python.exe" simulador_fermentacao.py
if errorlevel 1 (
  echo Falha na simulacao.
  pause
  exit /b 1
)

echo.
echo Pronto. Saidas em: saidas\5_1_referencia ... saidas\5_6_indicadores
pause
