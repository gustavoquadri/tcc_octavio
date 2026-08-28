# Guia de apresentação do código — TCC Octávio

Roteiro para estudar e defender o simulador. O código-fonte principal é `simulador_fermentacao.py`; a janela é `interface.py`. As explicações ficam **neste arquivo** — use o código só para localizar as funções citadas abaixo.

---

## 1. O que dizer em 60 segundos

> “Desenvolvi um simulador em Python da fermentação alcoólica em batelada por *Saccharomyces cerevisiae*. O modelo acompanha três variáveis: biomassa (X), substrato (S) e etanol (P). A velocidade específica de crescimento μ pode ser calculada por Monod ou por Andrews. Com μ, resolvo o sistema de equações diferenciais com `solve_ivp` (SciPy), gero gráficos (Matplotlib) e calculo indicadores de desempenho. Os parâmetros vêm de Zentou et al. (2019). O escopo é didático e simplificado: sem morte celular, sem inibição por etanol, sem oxigênio e sem CFD.”

---

## 2. Como o programa está organizado

```
definir_cenarios()          → monta S0, X0 e o modelo de cada caso
        ↓
calcular_monod / andrews    → calcula μ a partir de S
        ↓
sistema_edo                 → devolve dX/dt, dS/dt, dP/dt
        ↓
simular_cenario             → resolve no tempo com solve_ivp
        ↓
calcular_indicadores        → Sconsumido, Pproduzido, Yp/s, Qp
calcular_residuos_balanco   → checa coerência do balanço (Cap. 5.2)
        ↓
gerar / comparar gráficos   → figuras do Cap. 5
executar()                  → roda tudo e grava em saidas/
```

A interface (`interface.py`) **não refaz a matemática**. Ela só chama `simular_cenario` e `calcular_indicadores` e mostra o resultado na tela.

| Arquivo | Papel na banca |
|---|---|
| `simulador_fermentacao.py` | Núcleo do Apêndice A — explicar este |
| `interface.py` | Demonstração ao vivo |
| `gerar_resultados.bat` | Regenera figuras/tabelas |
| `saidas/5_1` … `5_6` | O que foi colado no Cap. 5 |

---

## 3. Variáveis e unidades (fixir de cabeça)

| Símbolo | Significado | Unidade |
|---|---|---|
| X | concentração de biomassa | g/L |
| S | concentração de substrato | g/L |
| P | concentração de etanol | g/L |
| t | tempo | h |
| μ | velocidade específica de crescimento | h⁻¹ |
| μmáx | μ máxima | h⁻¹ |
| Ks | constante de saturação | g/L |
| Ki | constante de inibição por substrato (Andrews) | g/L |
| Yx/s | rendimento biomassa/substrato | g/g |
| Yp/s | rendimento etanol/substrato | g/g |

Condições iniciais: `X(0)=X0`, `S(0)=S0`, `P(0)=P0` (no TCC, P0 = 0).

No código, esses valores entram na classe `Parametros` e nos dicionários `MONOD` / `ANDREWS` (topo de `simulador_fermentacao.py`).

---

## 4. Fórmulas — o que a banca pode pedir no quadro

### Monod  → função `calcular_monod`

\[
\mu = \mu_{máx}\,\frac{S}{K_S + S}
\]

Quando S é alto, μ → μmáx. Quando S é baixo, μ cai.

### Andrews  → função `calcular_andrews`

\[
\mu = \mu_{máx}\,\frac{S}{K_S + S + S^2/K_i}
\]

O termo \(S^2/K_i\) representa **inibição por substrato** (não por etanol).

### Sistema de EDOs  → função `sistema_edo`

\[
\frac{dX}{dt} = \mu X
\qquad
\frac{dS}{dt} = -\frac{1}{Y_{X/S}}\frac{dX}{dt}
\qquad
\frac{dP}{dt} = -Y_{P/S}\frac{dS}{dt}
\]

Em palavras: a biomassa cresce; o substrato é consumido na proporção do crescimento; o etanol sobe porque `dS/dt` é negativo e o sinal na equação de P inverte isso.

**Detalhes úteis ao abrir essa função:**

- O argumento `y` é o vetor `[X, S, P]`.
- O argumento `t` existe porque o SciPy exige; o modelo **não depende explicitamente do tempo** (autônomo).
- `P` **não entra** no cálculo de μ neste escopo (não há inibição por produto).

### Indicadores  → função `calcular_indicadores`

\[
S_{consumido}=S_0-S_f
\qquad
P_{produzido}=P_f-P_0
\qquad
Y_{P/S}^{ap}=\frac{P_{produzido}}{S_{consumido}}
\qquad
Q_P=\frac{P_{produzido}}{t_f}
\]

Usa o **último ponto** da série. Se o substrato não for consumido, Yp/s aparente vira NaN (evita divisão por zero).

### Resíduo do balanço  → função `calcular_residuos_balanco` (Cap. 5.2)

\[
R_X=(X-X_0)-Y_{X/S}(S_0-S)
\qquad
R_P=(P-P_0)-Y_{P/S}(S_0-S)
\]

Se a integração estiver coerente, \(R_X\) e \(R_P\) ficam perto de zero (~10⁻¹⁴ g/L).

---

## 5. Funções do código — o que cada uma faz

### `Parametros` + `MONOD` / `ANDREWS`

Empacotam os números de Zentou. **Importante:** Monod e Andrews **não** compartilham o mesmo μmáx/Ks/Y. Cada modelo foi ajustado à parte no artigo. Na comparação M1×M2 você compara dois modelos calibrados, não a mesma cinética com “só a fórmula diferente”.

Constantes no topo do arquivo: `TF_H = 72`, `N_PONTOS = 1441`, `METODO_SOLVER = "RK45"`, `RTOL`, `ATOL`.

### `definir_cenarios()`

Monta os casos do TCC:

| Nome | Modelo | S0 | X0 |
|---|---|---|---|
| referência / M1 | Monod | 150 | 1,0 |
| S1 / S2 | Monod | 75 / 225 | 1,0 |
| X1 / X2 | Monod | 150 | 0,5 / 2,0 |
| M2 | Andrews | 150 | 1,0 |

`tf` = 72 h em todos. M1 é cópia da referência (mesmo Monod) para a comparação com M2.

### `simular_cenario`

1. Monta a malha de tempo (`np.linspace` de 0 a 72 h, 1441 pontos)
2. Chama `solve_ivp` com método RK45, passando `sistema_edo`
3. Monta um DataFrame com `tempo_h`, `X_g_L`, `S_g_L`, `P_g_L`
4. Guarda diagnóstico (sucesso, `nfev`, mínimos de X/S/P)

**Frase útil:** “O solver escolhe o passo interno; `t_eval` só define onde eu gravo a solução para tabela e gráfico.”

Sobre S negativo ~10⁻⁹: quando o substrato zera, o RK45 pode “passar” levemente do zero (ordem do `ATOL`). O valor **não é forçado a zero** de propósito. Nos resultados do texto, trate S final como ~0 g/L.

### Gráficos

| Função | O que gera |
|---|---|
| `gerar_grafico_conjunto` | X, S e P no mesmo eixo |
| `gerar_graficos` | individuais + conjunto |
| `gerar_grafico_residuos` | R_X e R_P |
| `comparar_cenarios` | três eixos (X, S, P) |
| `comparar_variavel` | uma figura por variável |

Estilos de linha diferentes (`-`, `--`, `-.`) ajudam se a figura for impressa em P&B.

### `executar()`

Orquestra o lote completo: simula todos os cenários → salva CSVs em `tabelas_completas` → `_exportar_cap5` grava as pastas `5_1` … `5_6` → imprime o resumo no terminal.

### Interface (`interface.py`)

| Parte | O que faz |
|---|---|
| `_montar_parametros` | junta S0/X0 da tela com `MONOD` ou `ANDREWS` |
| `_simular` | chama o mesmo núcleo do simulador |
| `_desenhar_comparacao` | re-simula os cenários **fixos** do TCC (S1/ref/S2 etc.), não o que está digitado no painel |
| `_preenchendo` | evita que o preenchimento automático do combo marque “Personalizado” |

P0 e tf continuam fixos (0 g/L e 72 h). μmáx, Ks, Ki e rendimentos **não** são editáveis na janela.

---

## 6. Perguntas frequentes da banca

**“Por que Python e SciPy?”**  
Linguagem acessível, bibliotecas científicas maduras; `solve_ivp` resolve EDOs com método bem estabelecido (RK45).

**“O que é RK45?”**  
Runge–Kutta de ordem 4(5), método explícito adequado a sistemas não rígidos como este. Está em `METODO_SOLVER` no topo do simulador.

**“Por que S fica um pouco negativo (~10⁻⁹)?”**  
Overshoot numérico quando S → 0, na ordem do `ATOL`. Não forçamos a zero. Nos gráficos/tabelas do texto, S final ≈ 0.

**“O resíduo do balanço é zero?”**  
Praticamente: ordem de 10⁻¹⁴ g/L (`calcular_residuos_balanco` / pasta `5_2_consistencia`).

**“Por que X1 e X2 têm o mesmo etanol final?”**  
Em 72 h o substrato acaba nos dois. Com o mesmo S0 e o mesmo Yp/s, P final é o mesmo. Muda a **velocidade** das curvas.

**“Por que Andrews parece mais rápido que Monod na referência?”**  
Porque os **parâmetros ajustados** em `ANDREWS` são diferentes dos de `MONOD` (μmáx maior no Zentou). Não é só a forma da equação.

**“O simulador substitui experimento?”**  
Não. Ferramenta didática/analítica. Sem morte celular, inibição por etanol, T/pH dinâmicos, oxigênio etc.

**“Onde está a interface no texto?”**  
Já implementada (`interface.py`). Nas Perspectivas, citar como recurso já disponível.

**“Você usou código de terceiros / GitHub?”**  
Não. Implementação própria com base na metodologia do TCC e nas bibliotecas NumPy, SciPy, Matplotlib, Pandas e Tkinter.

---

## 7. Demonstração ao vivo (2–3 minutos)

1. Rodar `abrir_interface.bat` (ou `python interface.py`).
2. Mostrar o cenário **Referência** → curvas X↑, S↓, P↑.
3. Trocar para **S2** → etanol final maior.
4. Trocar para **X2** → processo mais rápido, P final parecido.
5. Clicar **Comparar Monod e Andrews**.
6. Se pedirem o lote completo: `gerar_resultados.bat` e abrir `saidas/`.

---

## 8. O que NÃO está no modelo

- Morte celular  
- Inibição por etanol  
- Manutenção metabólica  
- Oxigênio dissolvido / aeração  
- CFD / geometria de impelidor  
- Controle industrial / planta completa  

Isso não é falha: é o **recorte** definido no TCC.

---

## 9. Roteiro de estudo (antes da banca)

1. Decorar as três EDOs e as duas cinéticas.  
2. Abrir no editor, nesta ordem: `calcular_monod` → `calcular_andrews` → `sistema_edo` → `simular_cenario` → `calcular_indicadores`. Explicar cada uma em voz alta.  
3. Rodar a interface uma vez e ensaiar a demo.  
4. Revisar a tabela de indicadores da referência (X≈43, P≈63, Yp/s=0,420).  
5. Ensaiar as respostas do S negativo e do “mesmo P em X1/X2”.

Se a banca pedir “mostre no código”, use **Ctrl+F** pelos nomes das funções deste guia — não depende de comentários extras no arquivo.
