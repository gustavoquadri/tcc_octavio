"""
Interface Tkinter do simulador.

Reutiliza simular_cenario e calcular_indicadores de simulador_fermentacao.py.
Não altera o modelo: só escolhe cenário/modelo, S0 e X0.

    .\\.venv\\Scripts\\python.exe interface.py
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib

matplotlib.use("TkAgg")

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from simulador_fermentacao import (
    ANDREWS,
    MONOD,
    ROTULO_CENARIO,
    TF_H,
    Parametros,
    calcular_indicadores,
    definir_cenarios,
    simular_cenario,
)

ORDEM_CENARIOS = ["referencia", "S1", "S2", "X1", "X2", "M1", "M2"]


def _ler_numero(texto: str, nome: str) -> float:
    """Aceita vírgula ou ponto decimal."""

    texto = texto.strip().replace(",", ".")
    if not texto:
        raise ValueError(f"Informe um valor para {nome}.")
    valor = float(texto)
    if valor < 0:
        raise ValueError(f"{nome} não pode ser negativo.")
    return valor


def _montar_parametros(nome: str, modelo: str, x0: float, s0: float) -> Parametros:
    cinetica = dict(MONOD if modelo == "monod" else ANDREWS)
    return Parametros(nome=nome, X0=x0, S0=s0, **cinetica)


class JanelaSimulador(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Simulador de fermentação alcoólica — TCC")
        self.geometry("1180x720")
        self.minsize(960, 620)

        self.cenarios = definir_cenarios()
        self._preenchendo = False

        self._montar_layout()
        self._aplicar_cenario("referencia")
        self._simular()

    def _montar_layout(self) -> None:
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        painel = ttk.Frame(self, padding=12)
        painel.grid(row=0, column=0, sticky="ns")

        grafico = ttk.Frame(self, padding=(0, 12, 12, 12))
        grafico.grid(row=0, column=1, sticky="nsew")
        grafico.columnconfigure(0, weight=1)
        grafico.rowconfigure(0, weight=1)

        ttk.Label(painel, text="Fermentação alcoólica em batelada", font=("Segoe UI", 12, "bold")).pack(anchor="w")
        ttk.Label(
            painel,
            text="Saccharomyces cerevisiae  |  X, S e P  |  0 a 72 h",
            wraplength=280,
        ).pack(anchor="w", pady=(0, 10))

        ttk.Label(painel, text="Cenário do TCC").pack(anchor="w")
        self.var_cenario = tk.StringVar()
        combo = ttk.Combobox(
            painel,
            textvariable=self.var_cenario,
            state="readonly",
            width=36,
            values=[self._texto_cenario(chave) for chave in ORDEM_CENARIOS] + ["Personalizado"],
        )
        combo.pack(anchor="w", pady=(2, 8))
        combo.bind("<<ComboboxSelected>>", self._ao_escolher_cenario)

        ttk.Label(painel, text="Modelo cinético").pack(anchor="w")
        self.var_modelo = tk.StringVar(value="monod")
        modelos = ttk.Frame(painel)
        modelos.pack(anchor="w", pady=(2, 8))
        ttk.Radiobutton(
            modelos, text="Monod", value="monod", variable=self.var_modelo, command=self._marcar_personalizado
        ).pack(side="left")
        ttk.Radiobutton(
            modelos, text="Andrews", value="andrews", variable=self.var_modelo, command=self._marcar_personalizado
        ).pack(side="left", padx=(12, 0))

        self.entrada_s0 = self._campo(painel, "S0 — substrato inicial (g/L)")
        self.entrada_x0 = self._campo(painel, "X0 — biomassa inicial (g/L)")

        ttk.Label(painel, text="P0 = 0 g/L    tf = 72 h", foreground="#444").pack(anchor="w", pady=(4, 8))
        ttk.Label(
            painel,
            text="μmáx, Ks, Ki, Yx/s e Yp/s vêm de Zentou et al. (2019) e não são editados aqui.",
            wraplength=280,
            foreground="#444",
        ).pack(anchor="w", pady=(0, 10))

        ttk.Button(painel, text="Simular", command=self._simular).pack(fill="x", pady=2)
        ttk.Button(painel, text="Comparar S1, referência e S2", command=self._comparar_s0).pack(fill="x", pady=2)
        ttk.Button(painel, text="Comparar X1, referência e X2", command=self._comparar_x0).pack(fill="x", pady=2)
        ttk.Button(painel, text="Comparar Monod e Andrews", command=self._comparar_modelos).pack(fill="x", pady=2)

        ttk.Separator(painel, orient="horizontal").pack(fill="x", pady=12)
        ttk.Label(painel, text="Indicadores", font=("Segoe UI", 10, "bold")).pack(anchor="w")

        self.texto_indicadores = tk.Text(painel, width=36, height=14, wrap="word", state="disabled")
        self.texto_indicadores.pack(fill="both", expand=True, pady=(4, 0))

        self.fig = Figure(figsize=(7.2, 6.0), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=grafico)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        barra = ttk.Frame(grafico)
        barra.grid(row=1, column=0, sticky="ew")
        NavigationToolbar2Tk(self.canvas, barra)

    def _campo(self, pai: ttk.Frame, rotulo: str) -> ttk.Entry:
        ttk.Label(pai, text=rotulo).pack(anchor="w")
        entrada = ttk.Entry(pai, width=18)
        entrada.pack(anchor="w", pady=(2, 8))
        entrada.bind("<KeyRelease>", lambda _evento: self._marcar_personalizado())
        return entrada

    def _texto_cenario(self, chave: str) -> str:
        p = self.cenarios[chave]
        modelo = "Monod" if p.modelo == "monod" else "Andrews"
        return f"{ROTULO_CENARIO[chave]}  |  {modelo}  |  S0={p.S0:g}  X0={p.X0:g}"

    def _chave_cenario(self) -> str | None:
        texto = self.var_cenario.get()
        if texto == "Personalizado":
            return None
        for chave in ORDEM_CENARIOS:
            if self._texto_cenario(chave) == texto:
                return chave
        return None

    def _ao_escolher_cenario(self, _evento=None) -> None:
        chave = self._chave_cenario()
        if chave:
            self._aplicar_cenario(chave)
            self._simular()

    def _aplicar_cenario(self, chave: str) -> None:
        parametros = self.cenarios[chave]
        self._preenchendo = True
        self.var_cenario.set(self._texto_cenario(chave))
        self.var_modelo.set(parametros.modelo)
        self.entrada_s0.delete(0, tk.END)
        self.entrada_s0.insert(0, f"{parametros.S0:g}".replace(".", ","))
        self.entrada_x0.delete(0, tk.END)
        self.entrada_x0.insert(0, f"{parametros.X0:g}".replace(".", ","))
        self._preenchendo = False

    def _marcar_personalizado(self) -> None:
        if not self._preenchendo:
            self.var_cenario.set("Personalizado")

    def _parametros_da_janela(self) -> Parametros:
        s0 = _ler_numero(self.entrada_s0.get(), "S0")
        x0 = _ler_numero(self.entrada_x0.get(), "X0")
        if s0 == 0 and x0 == 0:
            raise ValueError("S0 e X0 não podem ser ambos zero.")
        chave = self._chave_cenario() or "personalizado"
        return _montar_parametros(chave, self.var_modelo.get(), x0, s0)

    def _simular(self) -> None:
        try:
            parametros = self._parametros_da_janela()
        except ValueError as erro:
            messagebox.showerror("Valor inválido", str(erro))
            return

        tabela, diagnostico = simular_cenario(parametros)
        indicadores = calcular_indicadores(tabela, parametros)
        self._desenhar_conjunto(tabela, parametros)
        self._mostrar_indicadores(parametros, indicadores, diagnostico)

    def _desenhar_conjunto(self, tabela, parametros: Parametros) -> None:
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.plot(tabela["tempo_h"], tabela["X_g_L"], label="Biomassa X", linewidth=2.0, linestyle="-")
        ax.plot(tabela["tempo_h"], tabela["S_g_L"], label="Substrato S", linewidth=2.0, linestyle="--")
        ax.plot(tabela["tempo_h"], tabela["P_g_L"], label="Etanol P", linewidth=2.0, linestyle="-.")
        ax.set_xlabel("Tempo (h)")
        ax.set_ylabel("Concentração (g/L)")
        rotulo = ROTULO_CENARIO.get(parametros.nome, "personalizado")
        modelo = "Monod" if parametros.modelo == "monod" else "Andrews"
        ax.set_title(f"{rotulo} — {modelo}  |  X, S e P")
        ax.set_xlim(0.0, parametros.tf)
        ax.grid(True, alpha=0.3)
        ax.legend()
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def _desenhar_comparacao(self, nomes: list[str], rotulos: dict[str, str], titulo: str) -> None:
        resultados = {}
        for nome in nomes:
            tabela, _diagnostico = simular_cenario(self.cenarios[nome])
            resultados[nome] = tabela

        self.fig.clear()
        eixos = self.fig.subplots(nrows=3, ncols=1, sharex=True)
        estilos = ["-", "--", "-."]
        variaveis = [
            ("X_g_L", "Biomassa X (g/L)"),
            ("S_g_L", "Substrato S (g/L)"),
            ("P_g_L", "Etanol P (g/L)"),
        ]

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
            ax.set_xlim(0.0, TF_H)
            ax.grid(True, alpha=0.3)
            ax.legend()

        eixos[-1].set_xlabel("Tempo (h)")
        self.fig.suptitle(titulo)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self._escrever_texto(
            f"{titulo}\n\n"
            "Comparação com os cenários fixos do TCC.\n"
            "Para um caso pontual, use Simular."
        )

    def _comparar_s0(self) -> None:
        self._desenhar_comparacao(
            nomes=["S1", "referencia", "S2"],
            rotulos={
                "S1": "S1 (S0 = 75 g/L)",
                "referencia": "Referência (S0 = 150 g/L)",
                "S2": "S2 (S0 = 225 g/L)",
            },
            titulo="Comparação entre S1, referência e S2",
        )

    def _comparar_x0(self) -> None:
        self._desenhar_comparacao(
            nomes=["X1", "referencia", "X2"],
            rotulos={
                "X1": "X1 (X0 = 0,5 g/L)",
                "referencia": "Referência (X0 = 1,0 g/L)",
                "X2": "X2 (X0 = 2,0 g/L)",
            },
            titulo="Comparação entre X1, referência e X2",
        )

    def _comparar_modelos(self) -> None:
        self._desenhar_comparacao(
            nomes=["M1", "M2"],
            rotulos={"M1": "M1 — Monod", "M2": "M2 — Andrews"},
            titulo="Comparação entre Monod (M1) e Andrews (M2)",
        )

    def _mostrar_indicadores(self, parametros: Parametros, indicadores: dict, diagnostico: dict) -> None:
        modelo = "Monod" if parametros.modelo == "monod" else "Andrews"
        s_min = diagnostico["S_min"]
        if diagnostico["houve_valor_negativo"]:
            nota_s = (
                f"S mínimo = {s_min:.4e} g/L\n"
                "(overshoot numérico ~ ATOL; não forçado a zero)"
            )
        else:
            nota_s = f"S mínimo = {s_min:.4e} g/L"

        texto = (
            f"Modelo: {modelo}\n"
            f"S0 = {parametros.S0:g} g/L    X0 = {parametros.X0:g} g/L\n"
            f"tf = {parametros.tf:g} h\n"
            "\n"
            f"X final = {indicadores['X_final_g_L']:.4f} g/L\n"
            f"S final = {indicadores['S_final_g_L']:.4f} g/L\n"
            f"P final = {indicadores['P_final_g_L']:.4f} g/L\n"
            "\n"
            f"Sconsumido = {indicadores['Sconsumido_g_L']:.4f} g/L\n"
            f"Pproduzido = {indicadores['Pproduzido_g_L']:.4f} g/L\n"
            f"Yp/s aparente = {indicadores['Yps_aparente_g_g']:.4f} g/g\n"
            f"Qp = {indicadores['Qp_g_L_h']:.4f} g/L/h\n"
            "\n"
            f"Solver: {'ok' if diagnostico['solver_sucesso'] else 'falhou'}"
            f"  |  nfev = {diagnostico['nfev']}\n"
            f"{nota_s}"
        )
        self._escrever_texto(texto)

    def _escrever_texto(self, texto: str) -> None:
        self.texto_indicadores.configure(state="normal")
        self.texto_indicadores.delete("1.0", tk.END)
        self.texto_indicadores.insert("1.0", texto)
        self.texto_indicadores.configure(state="disabled")


def abrir_interface() -> None:
    JanelaSimulador().mainloop()


if __name__ == "__main__":
    abrir_interface()
