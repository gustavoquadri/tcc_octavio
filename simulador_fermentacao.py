"""
Simulador de fermentação alcoólica em batelada (S. cerevisiae).

Variáveis: X (biomassa), S (substrato), P (etanol), em g/L; tempo em horas.
Cinéticas: Monod e Andrews. Parâmetros: Zentou et al. (2019).

    python simulador_fermentacao.py   # tabelas e figuras do Cap. 5
    python interface.py               # interface Tkinter

Saídas em saidas/ (pastas 5_1 … 5_6).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


PASTA_RAIZ = Path(__file__).resolve().parent
PASTA_SAIDAS = PASTA_RAIZ / "saidas"
PASTA_5_1 = PASTA_SAIDAS / "5_1_referencia"
PASTA_5_2 = PASTA_SAIDAS / "5_2_consistencia"
PASTA_5_3 = PASTA_SAIDAS / "5_3_substrato_S0"
PASTA_5_4 = PASTA_SAIDAS / "5_4_biomassa_X0"
PASTA_5_5 = PASTA_SAIDAS / "5_5_monod_andrews"
PASTA_5_6 = PASTA_SAIDAS / "5_6_indicadores"
PASTA_COMPLETAS = PASTA_SAIDAS / "tabelas_completas"

# 0–72 h, 1441 pontos (passo de 0,05 h).
TF_H = 72.0
N_PONTOS = 1441

# Tolerâncias explícitas; overshoot de S na ordem do ATOL não é mascarado.
METODO_SOLVER = "RK45"
RTOL = 1e-6
ATOL = 1e-9


@dataclass(frozen=True)
class Parametros:
    """Parâmetros de um cenário. Ki só entra em Andrews."""

    nome: str
    modelo: str  # "monod" | "andrews"
    mu_max: float  # h⁻¹
    Ks: float  # g/L
    Ki: float | None  # g/L
    Yxs: float  # g/g
    Yps: float  # g/g
    X0: float  # g/L
    S0: float  # g/L
    P0: float  # g/L
    tf: float  # h


# Zentou et al. (2019) — conjuntos ajustados em separado para cada modelo.
MONOD = dict(
    modelo="monod",
    mu_max=0.179,
    Ks=11.37,
    Ki=None,
    Yxs=0.280,
    Yps=0.420,
    P0=0.00,
    tf=TF_H,
)

ANDREWS = dict(
    modelo="andrews",
    mu_max=0.508,
    Ks=47.53,
    Ki=181.02,
    Yxs=0.286,
    Yps=0.431,
    P0=0.00,
    tf=TF_H,
)

ROTULO_CENARIO = {
    "referencia": "Referência",
    "S1": "S1",
    "S2": "S2",
    "X1": "X1",
    "X2": "X2",
    "M1": "M1",
    "M2": "M2",
}


def definir_cenarios() -> dict[str, Parametros]:
    """Cenários do TCC. M1 replica a referência; M2 usa Andrews."""

    return {
        "referencia": Parametros(nome="referencia", X0=1.00, S0=150.0, **MONOD),
        "S1": Parametros(nome="S1", X0=1.00, S0=75.0, **MONOD),
        "S2": Parametros(nome="S2", X0=1.00, S0=225.0, **MONOD),
        "X1": Parametros(nome="X1", X0=0.50, S0=150.0, **MONOD),
        "X2": Parametros(nome="X2", X0=2.00, S0=150.0, **MONOD),
        "M1": Parametros(nome="M1", X0=1.00, S0=150.0, **MONOD),
        "M2": Parametros(nome="M2", X0=1.00, S0=150.0, **ANDREWS),
    }


def calcular_monod(S: float, mu_max: float, Ks: float) -> float:
    """μ = μmáx · S / (Ks + S)  [h⁻¹]."""

    return mu_max * S / (Ks + S)


def calcular_andrews(S: float, mu_max: float, Ks: float, Ki: float) -> float:
    """μ = μmáx · S / [Ks + S + (S²/Ki)]  [h⁻¹]."""

    return mu_max * S / (Ks + S + (S**2) / Ki)


def sistema_edo(t: float, y: np.ndarray, parametros: Parametros) -> list[float]:
    """
    dX/dt = μ X
    dS/dt = −(1/Yx/s) · dX/dt
    dP/dt = −(Yp/s) · dS/dt
    """

    X, S, _P = y
    _ = t  # modelo autônomo; t exigido pelo solve_ivp

    if parametros.modelo == "monod":
        mu = calcular_monod(S, parametros.mu_max, parametros.Ks)
    elif parametros.modelo == "andrews":
        mu = calcular_andrews(S, parametros.mu_max, parametros.Ks, parametros.Ki)
    else:
        raise ValueError(f"Modelo cinético desconhecido: {parametros.modelo}")

    dX_dt = mu * X
    dS_dt = -(1.0 / parametros.Yxs) * dX_dt
    dP_dt = -parametros.Yps * dS_dt
    return [dX_dt, dS_dt, dP_dt]


def simular_cenario(
    parametros: Parametros,
    n_pontos: int = N_PONTOS,
    metodo: str = METODO_SOLVER,
    rtol: float = RTOL,
    atol: float = ATOL,
) -> tuple[pd.DataFrame, dict]:
    """Integra de 0 a tf com solve_ivp e devolve tabela + diagnóstico."""

    t_eval = np.linspace(0.0, parametros.tf, n_pontos)
    y0 = [parametros.X0, parametros.S0, parametros.P0]

    solucao = solve_ivp(
        fun=sistema_edo,
        t_span=(0.0, parametros.tf),
        y0=y0,
        method=metodo,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        args=(parametros,),
    )

    tabela = pd.DataFrame(
        {
            "tempo_h": solucao.t,
            "X_g_L": solucao.y[0],
            "S_g_L": solucao.y[1],
            "P_g_L": solucao.y[2],
        }
    )

    diagnostico = {
        "cenario": parametros.nome,
        "modelo": parametros.modelo,
        "solver_sucesso": bool(solucao.success),
        "mensagem_solver": solucao.message,
        "nfev": int(solucao.nfev),
        "metodo": metodo,
        "rtol": rtol,
        "atol": atol,
        "X_min": float(tabela["X_g_L"].min()),
        "S_min": float(tabela["S_g_L"].min()),
        "P_min": float(tabela["P_g_L"].min()),
        "houve_valor_negativo": bool(
            (tabela[["X_g_L", "S_g_L", "P_g_L"]] < 0).any().any()
        ),
    }

    if not solucao.success:
        print(
            f"[AVISO] Solver não convergiu no cenário {parametros.nome}: "
            f"{solucao.message}"
        )

    return tabela, diagnostico


def calcular_indicadores(tabela: pd.DataFrame, parametros: Parametros) -> dict:
    """
    Sconsumido = S0 − Sfinal
    Pproduzido = Pfinal − P0
    Yp/s aparente = Pproduzido / Sconsumido
    Qp = Pproduzido / tf
    """

    S_final = float(tabela["S_g_L"].iloc[-1])
    P_final = float(tabela["P_g_L"].iloc[-1])
    X_final = float(tabela["X_g_L"].iloc[-1])

    s_consumido = parametros.S0 - S_final
    p_produzido = P_final - parametros.P0
    yps_aparente = float("nan") if s_consumido == 0 else p_produzido / s_consumido
    qp = p_produzido / parametros.tf

    return {
        "cenario": parametros.nome,
        "modelo": parametros.modelo,
        "X0_g_L": parametros.X0,
        "S0_g_L": parametros.S0,
        "P0_g_L": parametros.P0,
        "tf_h": parametros.tf,
        "X_final_g_L": X_final,
        "S_final_g_L": S_final,
        "P_final_g_L": P_final,
        "Sconsumido_g_L": s_consumido,
        "Pproduzido_g_L": p_produzido,
        "Yps_aparente_g_g": yps_aparente,
        "Qp_g_L_h": qp,
    }


def calcular_residuos_balanco(tabela: pd.DataFrame, parametros: Parametros) -> pd.DataFrame:
    """
    R_X = (X − X0) − Yx/s · (S0 − S)
    R_P = (P − P0) − Yp/s · (S0 − S)
    """

    s_consumido = parametros.S0 - tabela["S_g_L"]
    r_x = (tabela["X_g_L"] - parametros.X0) - parametros.Yxs * s_consumido
    r_p = (tabela["P_g_L"] - parametros.P0) - parametros.Yps * s_consumido

    return pd.DataFrame(
        {
            "tempo_h": tabela["tempo_h"],
            "R_X_g_L": r_x,
            "R_P_g_L": r_p,
            "R_abs_max_g_L": np.maximum(np.abs(r_x), np.abs(r_p)),
        }
    )


def _configurar_estilo() -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.30,
            "savefig.dpi": 300,
            "figure.autolayout": True,
        }
    )


def gerar_grafico_conjunto(
    tabela: pd.DataFrame,
    parametros: Parametros,
    arquivo: Path,
    titulo: str | None = None,
) -> None:
    """X, S e P no mesmo eixo."""

    arquivo.parent.mkdir(parents=True, exist_ok=True)
    rotulo = ROTULO_CENARIO.get(parametros.nome, parametros.nome)
    if titulo is None:
        titulo = f"Cenário {rotulo} — X, S e P"

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.plot(tabela["tempo_h"], tabela["X_g_L"], label="Biomassa X", linewidth=2.0, linestyle="-")
    ax.plot(tabela["tempo_h"], tabela["S_g_L"], label="Substrato S", linewidth=2.0, linestyle="--")
    ax.plot(tabela["tempo_h"], tabela["P_g_L"], label="Etanol P", linewidth=2.0, linestyle="-.")
    ax.set_xlabel("Tempo (h)")
    ax.set_ylabel("Concentração (g/L)")
    ax.set_title(titulo)
    ax.set_xlim(0.0, parametros.tf)
    ax.legend()
    fig.savefig(arquivo)
    plt.close(fig)


def gerar_graficos(tabela: pd.DataFrame, parametros: Parametros, pasta: Path) -> None:
    """Figuras individuais (X, S, P) e o gráfico conjunto."""

    pasta.mkdir(parents=True, exist_ok=True)
    t = tabela["tempo_h"]
    prefixo = parametros.nome
    rotulo = ROTULO_CENARIO.get(parametros.nome, parametros.nome)

    series = [
        ("biomassa", tabela["X_g_L"], "Concentração de biomassa (g/L)", "C0", "-"),
        ("substrato", tabela["S_g_L"], "Concentração de substrato (g/L)", "C1", "-"),
        ("etanol", tabela["P_g_L"], "Concentração de etanol (g/L)", "C2", "-"),
    ]

    for sufixo, valores, ylabel, cor, estilo in series:
        fig, ax = plt.subplots(figsize=(8.5, 5.0))
        ax.plot(t, valores, color=cor, linestyle=estilo, linewidth=2.0)
        ax.set_xlabel("Tempo (h)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Cenário {rotulo} — {ylabel.split(' (')[0]}")
        ax.set_xlim(0.0, parametros.tf)
        fig.savefig(pasta / f"{prefixo}_{sufixo}.png")
        plt.close(fig)

    gerar_grafico_conjunto(tabela, parametros, pasta / f"{prefixo}_conjunto.png")


def gerar_grafico_residuos(residuos: pd.DataFrame, parametros: Parametros, arquivo: Path) -> None:
    """Resíduos R_X e R_P ao longo do tempo."""

    arquivo.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.plot(residuos["tempo_h"], residuos["R_X_g_L"], label="Resíduo de biomassa R_X", linewidth=2.0)
    ax.plot(
        residuos["tempo_h"],
        residuos["R_P_g_L"],
        label="Resíduo de etanol R_P",
        linewidth=2.0,
        linestyle="--",
    )
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("Tempo (h)")
    ax.set_ylabel("Resíduo do balanço (g/L)")
    rotulo = ROTULO_CENARIO.get(parametros.nome, parametros.nome)
    ax.set_title(f"Cenário {rotulo} — resíduo do balanço de massa")
    ax.set_xlim(0.0, parametros.tf)
    ax.legend()
    fig.savefig(arquivo)
    plt.close(fig)


def comparar_cenarios(
    resultados: dict[str, pd.DataFrame],
    nomes: list[str],
    rotulos: dict[str, str],
    titulo: str,
    arquivo: Path,
    tf: float = TF_H,
) -> None:
    """Três eixos (X, S, P) com estilos de linha distintos (melhor em P&B)."""

    estilos = ["-", "--", "-.", ":"]
    variaveis = [
        ("X_g_L", "Biomassa X (g/L)"),
        ("S_g_L", "Substrato S (g/L)"),
        ("P_g_L", "Etanol P (g/L)"),
    ]

    arquivo.parent.mkdir(parents=True, exist_ok=True)
    fig, eixos = plt.subplots(nrows=3, ncols=1, figsize=(8.5, 10.5), sharex=True)

    for ax, (coluna, ylabel) in zip(eixos, variaveis):
        for i, nome in enumerate(nomes):
            tabela = resultados[nome]
            ax.plot(
                tabela["tempo_h"],
                tabela[coluna],
                label=rotulos[nome],
                linewidth=2.0,
                linestyle=estilos[i % len(estilos)],
            )
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.0, tf)
        ax.legend()

    eixos[-1].set_xlabel("Tempo (h)")
    fig.suptitle(titulo, fontsize=14)
    fig.savefig(arquivo)
    plt.close(fig)


def comparar_variavel(
    resultados: dict[str, pd.DataFrame],
    nomes: list[str],
    rotulos: dict[str, str],
    coluna: str,
    ylabel: str,
    titulo: str,
    arquivo: Path,
    tf: float = TF_H,
) -> None:
    """Uma figura por variável."""

    estilos = ["-", "--", "-.", ":"]
    arquivo.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))

    for i, nome in enumerate(nomes):
        tabela = resultados[nome]
        ax.plot(
            tabela["tempo_h"],
            tabela[coluna],
            label=rotulos[nome],
            linewidth=2.0,
            linestyle=estilos[i % len(estilos)],
        )

    ax.set_xlabel("Tempo (h)")
    ax.set_ylabel(ylabel)
    ax.set_title(titulo)
    ax.set_xlim(0.0, tf)
    ax.legend()
    fig.savefig(arquivo)
    plt.close(fig)


def salvar_tabela(
    tabela: pd.DataFrame,
    caminho: Path,
    float_format: str = "%.6f",
) -> None:
    """CSV com ';' e vírgula decimal (Excel pt-BR)."""

    caminho.parent.mkdir(parents=True, exist_ok=True)
    tabela.to_csv(
        caminho,
        sep=";",
        decimal=",",
        index=False,
        encoding="utf-8-sig",
        float_format=float_format,
    )


COLUNAS_INDICADORES_TCC = [
    "cenario",
    "modelo",
    "X0_g_L",
    "S0_g_L",
    "P0_g_L",
    "tf_h",
    "X_final_g_L",
    "S_final_g_L",
    "P_final_g_L",
    "Sconsumido_g_L",
    "Pproduzido_g_L",
    "Yps_aparente_g_g",
    "Qp_g_L_h",
]


def _ordenar_cenarios(df: pd.DataFrame, ordem: list[str]) -> pd.DataFrame:
    saida = df.copy()
    saida["cenario"] = saida["cenario"].astype(str)
    saida["cenario"] = pd.Categorical(saida["cenario"], categories=ordem, ordered=True)
    return saida.sort_values("cenario").reset_index(drop=True)


def _tabela_indicadores_tcc(df: pd.DataFrame, ordem: list[str]) -> pd.DataFrame:
    return _ordenar_cenarios(df[COLUNAS_INDICADORES_TCC], ordem)


def _imprimir_resumo(indicadores: pd.DataFrame, diagnosticos: list[dict]) -> None:
    print("\n=== Indicadores finais ===")
    colunas = [
        "cenario",
        "modelo",
        "S0_g_L",
        "X0_g_L",
        "X_final_g_L",
        "S_final_g_L",
        "P_final_g_L",
        "Sconsumido_g_L",
        "Pproduzido_g_L",
        "Yps_aparente_g_g",
        "Qp_g_L_h",
    ]
    print(indicadores[colunas].to_string(index=False, float_format=lambda v: f"{v:10.4f}"))

    print("\n=== Diagnóstico numérico ===")
    for item in diagnosticos:
        print(
            f"{item['cenario']:>12}  sucesso={item['solver_sucesso']}  "
            f"nfev={item['nfev']:<5}  "
            f"negativo={item['houve_valor_negativo']}  "
            f"S_min={item['S_min']:.4e}  "
            f"X_min={item['X_min']:.4e}  "
            f"P_min={item['P_min']:.4e}"
        )

    com_overshoot = [d for d in diagnosticos if d["houve_valor_negativo"]]
    if com_overshoot:
        pior = min(com_overshoot, key=lambda d: d["S_min"])
        print(
            "\n[NOTA] Overshoot numérico de S após esgotamento "
            f"(ordem do ATOL={ATOL:.0e}). "
            f"Ex.: cenário {pior['cenario']} com S_min={pior['S_min']:.4e} g/L. "
            "Valores não foram forçados a zero — ver diagnostico_solver.csv."
        )


def _exportar_cap5(
    resultados: dict[str, pd.DataFrame],
    cenarios: dict[str, Parametros],
    df_indicadores: pd.DataFrame,
    df_diagnostico: pd.DataFrame,
) -> None:
    # 5.1
    gerar_graficos(resultados["referencia"], cenarios["referencia"], PASTA_5_1)
    salvar_tabela(resultados["referencia"], PASTA_5_1 / "referencia_serie_temporal.csv")
    salvar_tabela(
        df_indicadores.loc[df_indicadores["cenario"] == "referencia", COLUNAS_INDICADORES_TCC],
        PASTA_5_1 / "referencia_valores_finais.csv",
    )

    # 5.2 — conjunto + resíduo
    residuos = calcular_residuos_balanco(resultados["referencia"], cenarios["referencia"])
    salvar_tabela(residuos, PASTA_5_2 / "residuos_balanco_referencia.csv", float_format="%.6e")
    gerar_grafico_residuos(
        residuos,
        cenarios["referencia"],
        PASTA_5_2 / "referencia_residuo_balanco.png",
    )
    gerar_grafico_conjunto(
        resultados["referencia"],
        cenarios["referencia"],
        PASTA_5_2 / "referencia_conjunto_verificacao.png",
        titulo="Verificação da consistência — X, S e P (referência)",
    )

    resumo_consistencia = pd.DataFrame(
        [
            {
                "cenario": "referencia",
                "R_X_max_abs_g_L": float(np.abs(residuos["R_X_g_L"]).max()),
                "R_P_max_abs_g_L": float(np.abs(residuos["R_P_g_L"]).max()),
                "R_abs_max_g_L": float(residuos["R_abs_max_g_L"].max()),
                "S_min_g_L": float(resultados["referencia"]["S_g_L"].min()),
                "X_min_g_L": float(resultados["referencia"]["X_g_L"].min()),
                "P_min_g_L": float(resultados["referencia"]["P_g_L"].min()),
            }
        ]
    )
    salvar_tabela(resumo_consistencia, PASTA_5_2 / "resumo_consistencia.csv", float_format="%.6e")

    # 5.3
    rotulos_s0 = {
        "S1": "S1 (S0 = 75 g/L)",
        "referencia": "Referência (S0 = 150 g/L)",
        "S2": "S2 (S0 = 225 g/L)",
    }
    nomes_s0 = ["S1", "referencia", "S2"]
    comparar_cenarios(
        resultados,
        nomes=nomes_s0,
        rotulos=rotulos_s0,
        titulo="Comparação entre S1, referência e S2",
        arquivo=PASTA_5_3 / "comparacao_S0_conjunto.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_s0,
        rotulos=rotulos_s0,
        coluna="X_g_L",
        ylabel="Biomassa X (g/L)",
        titulo="Comparação de biomassa — S0 = 75, 150 e 225 g/L",
        arquivo=PASTA_5_3 / "comparacao_S0_biomassa.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_s0,
        rotulos=rotulos_s0,
        coluna="S_g_L",
        ylabel="Substrato S (g/L)",
        titulo="Comparação de substrato — S0 = 75, 150 e 225 g/L",
        arquivo=PASTA_5_3 / "comparacao_S0_substrato.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_s0,
        rotulos=rotulos_s0,
        coluna="P_g_L",
        ylabel="Etanol P (g/L)",
        titulo="Comparação de etanol — S0 = 75, 150 e 225 g/L",
        arquivo=PASTA_5_3 / "comparacao_S0_etanol.png",
    )
    salvar_tabela(
        _tabela_indicadores_tcc(df_indicadores[df_indicadores["cenario"].isin(nomes_s0)], nomes_s0),
        PASTA_5_3 / "indicadores_S0.csv",
    )

    # 5.4
    rotulos_x0 = {
        "X1": "X1 (X0 = 0,5 g/L)",
        "referencia": "Referência (X0 = 1,0 g/L)",
        "X2": "X2 (X0 = 2,0 g/L)",
    }
    nomes_x0 = ["X1", "referencia", "X2"]
    comparar_cenarios(
        resultados,
        nomes=nomes_x0,
        rotulos=rotulos_x0,
        titulo="Comparação entre X1, referência e X2",
        arquivo=PASTA_5_4 / "comparacao_X0_conjunto.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_x0,
        rotulos=rotulos_x0,
        coluna="X_g_L",
        ylabel="Biomassa X (g/L)",
        titulo="Comparação de biomassa — X0 = 0,5, 1,0 e 2,0 g/L",
        arquivo=PASTA_5_4 / "comparacao_X0_biomassa.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_x0,
        rotulos=rotulos_x0,
        coluna="S_g_L",
        ylabel="Substrato S (g/L)",
        titulo="Comparação de substrato — X0 = 0,5, 1,0 e 2,0 g/L",
        arquivo=PASTA_5_4 / "comparacao_X0_substrato.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_x0,
        rotulos=rotulos_x0,
        coluna="P_g_L",
        ylabel="Etanol P (g/L)",
        titulo="Comparação de etanol — X0 = 0,5, 1,0 e 2,0 g/L",
        arquivo=PASTA_5_4 / "comparacao_X0_etanol.png",
    )
    salvar_tabela(
        _tabela_indicadores_tcc(df_indicadores[df_indicadores["cenario"].isin(nomes_x0)], nomes_x0),
        PASTA_5_4 / "indicadores_X0.csv",
    )

    # 5.5
    rotulos_m = {"M1": "M1 — Monod", "M2": "M2 — Andrews"}
    nomes_m = ["M1", "M2"]
    comparar_cenarios(
        resultados,
        nomes=nomes_m,
        rotulos=rotulos_m,
        titulo="Comparação entre Monod (M1) e Andrews (M2)",
        arquivo=PASTA_5_5 / "comparacao_monod_andrews_conjunto.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_m,
        rotulos=rotulos_m,
        coluna="X_g_L",
        ylabel="Biomassa X (g/L)",
        titulo="Comparação de biomassa — Monod e Andrews",
        arquivo=PASTA_5_5 / "comparacao_monod_andrews_biomassa.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_m,
        rotulos=rotulos_m,
        coluna="S_g_L",
        ylabel="Substrato S (g/L)",
        titulo="Comparação de substrato — Monod e Andrews",
        arquivo=PASTA_5_5 / "comparacao_monod_andrews_substrato.png",
    )
    comparar_variavel(
        resultados,
        nomes=nomes_m,
        rotulos=rotulos_m,
        coluna="P_g_L",
        ylabel="Etanol P (g/L)",
        titulo="Comparação de etanol — Monod e Andrews",
        arquivo=PASTA_5_5 / "comparacao_monod_andrews_etanol.png",
    )
    salvar_tabela(
        _tabela_indicadores_tcc(df_indicadores[df_indicadores["cenario"].isin(nomes_m)], nomes_m),
        PASTA_5_5 / "indicadores_monod_andrews.csv",
    )

    # 5.6
    ordem_todos = ["referencia", "S1", "S2", "X1", "X2", "M1", "M2"]
    salvar_tabela(
        _tabela_indicadores_tcc(df_indicadores, ordem_todos),
        PASTA_5_6 / "indicadores_todos_cenarios.csv",
    )
    salvar_tabela(df_diagnostico, PASTA_5_6 / "diagnostico_solver.csv", float_format="%.6e")


def executar() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    _configurar_estilo()
    for pasta in (
        PASTA_5_1,
        PASTA_5_2,
        PASTA_5_3,
        PASTA_5_4,
        PASTA_5_5,
        PASTA_5_6,
        PASTA_COMPLETAS,
    ):
        pasta.mkdir(parents=True, exist_ok=True)

    cenarios = definir_cenarios()
    resultados: dict[str, pd.DataFrame] = {}
    lista_indicadores: list[dict] = []
    lista_diagnosticos: list[dict] = []

    print("Iniciando simulações (0 a 72 h, 1441 pontos)...\n")

    for nome, parametros in cenarios.items():
        print(f"  - {nome} ({parametros.modelo}, S0={parametros.S0:g} g/L, X0={parametros.X0:g} g/L)")
        tabela, diagnostico = simular_cenario(parametros)
        resultados[nome] = tabela
        lista_diagnosticos.append(diagnostico)

        indicadores = calcular_indicadores(tabela, parametros)
        indicadores.update(
            {
                "solver_sucesso": diagnostico["solver_sucesso"],
                "nfev": diagnostico["nfev"],
                "houve_valor_negativo": diagnostico["houve_valor_negativo"],
                "X_min_g_L": diagnostico["X_min"],
                "S_min_g_L": diagnostico["S_min"],
                "P_min_g_L": diagnostico["P_min"],
            }
        )
        lista_indicadores.append(indicadores)
        salvar_tabela(tabela, PASTA_COMPLETAS / f"{nome}.csv")

    df_indicadores = pd.DataFrame(lista_indicadores)
    df_diagnostico = pd.DataFrame(lista_diagnosticos)
    ordem_todos = ["referencia", "S1", "S2", "X1", "X2", "M1", "M2"]
    df_indicadores = _ordenar_cenarios(df_indicadores, ordem_todos)
    salvar_tabela(df_indicadores, PASTA_COMPLETAS / "indicadores.csv")
    salvar_tabela(df_diagnostico, PASTA_COMPLETAS / "diagnostico_solver.csv", float_format="%.6e")

    print("\nExportando figuras e tabelas do Capítulo 5...")
    _exportar_cap5(resultados, cenarios, df_indicadores, df_diagnostico)

    _imprimir_resumo(df_indicadores, lista_diagnosticos)
    print(f"\nSaídas do Cap. 5 em: {PASTA_SAIDAS}")
    print("  5_1_referencia / 5_2_consistencia / 5_3_substrato_S0 /")
    print("  5_4_biomassa_X0 / 5_5_monod_andrews / 5_6_indicadores")
    print(f"Tabelas completas em: {PASTA_COMPLETAS}")


if __name__ == "__main__":
    executar()
