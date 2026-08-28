# Simulador de fermentação alcoólica (TCC)

Batelada com *S. cerevisiae*: biomassa (X), substrato (S) e etanol (P); Monod ou Andrews; parâmetros de Zentou et al. (2019).

**Para estudar e apresentar o código na banca:** [`GUIA_APRESENTACAO.md`](GUIA_APRESENTACAO.md) (explicações ficam no guia; o código permanece enxuto).

## Rodar

Python 3.11–3.13 (evitar 3.15 por falta de wheel de SciPy/Matplotlib).

```bash
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
gerar_resultados.bat
abrir_interface.bat
```

## Saídas (`saidas/`)

Pastas `5_1_referencia` … `5_6_indicadores` alinhadas ao Capítulo 5; `tabelas_completas` tem as séries brutas.

## Arquivos

| Arquivo | Função |
|---|---|
| `simulador_fermentacao.py` | núcleo + exportação Cap. 5 |
| `interface.py` | Tkinter |
| `GUIA_APRESENTACAO.md` | roteiro de defesa do código |
| `gerar_resultados.bat` / `abrir_interface.bat` | atalhos |
