> **Global Solution 2026.2 | CitizenTech** · **Disciplina:** Differentiated Problem Solving

# 💊 Memora · Por que o horário importa

Modelagem matemática em Python que mostra como atrasar, esquecer ou repetir uma dose muda a concentração de um remédio no corpo ao longo do tempo.

## 📌 Descrição do projeto

Para quem tem Alzheimer, o remédio traz duas dúvidas todos os dias: *"já está na hora?"* e *"será que eu já tomei?"*. Esquecer a dose atrapalha o tratamento, e tomar duas vezes pode ser perigoso. Mas quanto um atraso de alguns minutos realmente muda? E uma dose esquecida, ou dobrada para "compensar"?

O **Memora** nasceu para ajudar a pessoa idosa e a família a cuidar dos horários. A solução completa reúne uma caixa de remédios inteligente, que lembra a pessoa na hora certa, e a Central do Cuidador, onde a família acompanha cada dose.

Este repositório traz a **parte matemática do Memora**. Com um modelo simples de concentração do remédio no corpo, usamos **derivadas** para estudar a queda da concentração entre as doses, **limites** para entender o comportamento depois de muitas doses e **integrais** para medir a exposição acumulada. Depois comparamos cinco cenários: dose no horário, atrasada, esquecida, dobrada e esquecida e compensada. O programa também lê os registros da Central do Cuidador e mostra como seria a curva com as retiradas registradas.

O objetivo é simples: mostrar com números por que registrar os horários importa e por que um único percentual não conta toda a história da rotina de remédios.

> **Importante:** todos os parâmetros são **fictícios e didáticos**, não são de um remédio real. O modelo não substitui orientação médica e os resultados não comprovam efeito clínico.

---

## 🧮 Modelo e hipóteses

A concentração é representada por `C(t)`, em unidades arbitrárias (u.a.), com `t` em horas. O modelo considera um único compartimento, entrada instantânea do remédio e eliminação proporcional à concentração:

```
dC/dt = −kC   ⇒   C(t) = C₀ · e^(−kt)
```

Com várias doses, a concentração é a soma do que sobrou de cada uma: `C(t) = Σ dᵢ · e^(−k(t − tᵢ))`, somando só as doses que já aconteceram.

| Parâmetro | Valor | Papel no modelo |
|---|---|---|
| Meia-vida | 6 h | Tempo para a concentração cair pela metade |
| Constante de eliminação `k` | ln(2)/6 ≈ 0,1155 h⁻¹ | Velocidade da queda |
| Incremento `C₀` | 10 u.a. | Quanto cada dose aumenta a concentração |
| Intervalo entre doses `τ` | 12 h | 2 doses por dia nos cenários |
| Fração que sobra `r = e^(−kτ)` | 0,25 | O que resta de uma dose depois de 12 h |
| Referências inferior e superior | 3 e 16 u.a. | Faixa arbitrária para comparar as curvas |
| Referência para a média | 5 u.a. | Exemplo algébrico, sem meta clínica |

---

## 📊 O que o código calcula

| Parte | Conceito | O que mostra |
|---|---|---|
| 1 | **Derivada** | Como a concentração cai entre as doses e quanto o vale diminui com atrasos de 10, 20, 30, 60 e 180 min (cálculo exato x aproximação linear) |
| 2 | **Limites** | Os picos e vales depois de muitas doses (série geométrica), quantas doses são precisas para chegar a 99% do patamar e o que acontece se o intervalo tender a zero |
| 3 | **Integrais** | A área sob a curva (AUC), a concentração média e a comparação dos 5 cenários em uma janela de 72 h |
| 4 | **Registros do Memora** | Lê o `dados.json` da Central do Cuidador e compara a curva das retiradas registradas com a curva ideal |

### Principais resultados

| Cenário | AUC (u.a.·h) | Δ AUC | Pico (u.a.) | Vale (u.a.) | Horas abaixo de 3 | Horas acima de 16 |
|---|---|---|---|---|---|---|
| No horário | 519,34 | 0,00% | 13,33 | 3,33 | 0,00 | 0,00 |
| Atraso de 3 h | 518,78 | −0,11% | 14,37 | 2,36 | 2,09 | 0,00 |
| Dose esquecida | 434,13 | −16,41% | 13,33 | 0,83 | 11,97 | 0,00 |
| Dose dobrada | 604,55 | +16,41% | 23,33 | 3,33 | 0,00 | 3,27 |
| Esquecida e compensada | 515,28 | −0,78% | 20,83 | 0,83 | 11,09 | 2,28 |

O cenário **esquecida e compensada** tem quase a mesma AUC do cenário no horário, mas passa mais de 11 horas abaixo da referência inferior e mais de 2 horas acima da superior. Ou seja: **uma área parecida não garante uma curva parecida**.

Outros resultados:

- A concentração tende a um ciclo entre **13,33 u.a.** (pico) e **3,33 u.a.** (vale), alcançado com 99% depois de **4 doses**.
- Com estes parâmetros, um atraso de **26,6 min** derruba o vale em 5%.
- No exemplo da Central do Cuidador, **83,3%** das retiradas foram no horário, e a curva hipotética fica com **91,3%** da AUC de referência.

---

## 📈 Gráficos gerados

Ao executar o programa, seis gráficos são salvos na pasta `graficos/`:

| Arquivo | Conteúdo |
|---|---|
| `fig1_derivada.png` | Uma dose caindo pela metade a cada 6 h e a taxa de variação (derivada) |
| `fig2_limite.png` | Picos e vales convergindo e o pico limite quando o intervalo diminui (limites) |
| `fig3_cenarios.png` | Os cinco cenários, um por um e todos juntos |
| `fig4_comparacao.png` | Exposição acumulada em 72 h (integral) e horas fora da faixa de referência |
| `fig5_tolerancia_adesao.png` | Quanto o vale cai com o atraso e a média teórica pela fração de doses tomadas |
| `fig6_dados_memora.png` | Curva com os registros da Central do Cuidador |

---

## ▶️ Como executar

### Pré-requisitos

- Python 3 instalado
  Verifique com: `python --version` ou `python3 --version`
- Bibliotecas **numpy** e **matplotlib**:

```bash
pip install numpy matplotlib
# ou: python -m pip install numpy matplotlib
# ou: py -m pip install numpy matplotlib
```

### 🟢 Pelo terminal (qualquer sistema)

```bash
# Clone o repositório
git clone https://github.com/CitizenTech-Memora/modelagem-matematica.git

# Entre na pasta
cd modelagem-matematica

# Execute
python memora.py    # Windows (ou: py memora.py)
python3 memora.py   # Linux / Mac

# Para também abrir as janelas dos gráficos
python memora.py --mostrar
```

### 🟡 Pelo PyCharm

1. Abra o PyCharm e clique em **File > Open**
2. Selecione a pasta do projeto (modelagem-matematica)
3. Instale **numpy** e **matplotlib** pelo terminal do PyCharm, se ainda não tiver
4. No painel lateral, clique com o botão direito em `memora.py`
5. Clique em **Run 'memora'**
6. O terminal integrado mostra as análises e os gráficos são salvos na pasta `graficos/`

### 🔵 Pelo VS Code

1. Abra o VS Code e clique em **File > Open Folder**
2. Selecione a pasta do projeto (modelagem-matematica)
3. Instale a extensão Python da Microsoft, se ainda não tiver
   (Extensões > pesquisar "Python" > instalar)
4. Instale **numpy** e **matplotlib** pelo terminal, se ainda não tiver
5. Abra o arquivo `memora.py`
6. Clique no botão **▶ Run Python File** no canto superior direito ou use o atalho `Ctrl + F5`

O programa não precisa da Central do Cuidador para rodar. Se não encontrar um `dados.json` na mesma pasta, ele usa o mesmo exemplo da opção [8] da Central (Dona Lurdes e o Remédio A).

---

## 🗂️ Estrutura do repositório

```
modelagem-matematica/
├── memora.py                           ← Código com o modelo, as análises e os gráficos
├── Relatorio.pdf                       ← Relatório técnico com as contas e a interpretação
└── README.md                           ← Este arquivo
```

## 💬 Organização interna do código

```
memora.py
│
├── Parâmetros do modelo (fictícios)
│   ├── T_MEIA_VIDA, K ← meia-vida de 6 h e constante de eliminação
│   ├── C0, TAU, R ← incremento por dose, intervalo de 12 h e fração que sobra
│   ├── LIM_MIN, LIM_MAX, C_MEDIA_MIN ← referências arbitrárias de comparação
│   └── TOLERANCIA_MINUTOS ← 20 min, a mesma da Central do Cuidador
│
├── Funções do modelo
│   ├── concentracao() ← soma o que sobrou de cada dose
│   ├── derivada_exata() / derivada_numerica() ← taxa de variação e conferência
│   ├── auc_exata() / auc_numerica() ← integral calculada à mão e por trapézios
│   ├── extremos_e_tempos() ← picos, vales e horas fora da faixa
│   └── pico_n(), vale_n(), patamar_pico(), patamar_vale() ← série geométrica e limites
│
├── Cenários
│   ├── montar_cenarios() ← no horário, atraso, esquecida, dobrada e compensada
│   └── analisar_cenarios() ← AUC, picos, vales e horas fora da faixa de cada um
│
├── Ligação com a Central do Cuidador
│   ├── exemplo_em_memoria() ← mesmo exemplo da opção [8] do Memora.py
│   ├── carregar_eventos_memora() ← lê o dados.json, se existir
│   └── analisar_dados_memora() ← curva das retiradas registradas x curva ideal
│
├── calcular() ← reúne todos os resultados
├── imprimir() ← mostra as 4 partes no terminal
├── grafico_*() ← gera os 6 gráficos
├── verificar_calculos() ← testes automáticos das propriedades matemáticas
│
└── main() ← executa tudo e salva os resultados
```

## 💾 Arquivos gerados

- **`graficos/`:** pasta com os seis gráficos em PNG, criada automaticamente.
- **`resultados.json`:** todos os valores calculados, para consulta.

Os dois ficam na **mesma pasta do `memora.py`**. Se houver um `dados.json` da Central do Cuidador nessa pasta, o programa usa os registros dele na Parte 4.

---

## ⚠️ Limitações desta versão

Este é um **modelo acadêmico e didático**. Nesta versão:

- Os parâmetros são fictícios e não representam nenhum remédio real.
- O modelo não inclui absorção gradual nem diferenças entre pessoas.
- As referências de 3, 16 e 5 u.a. foram escolhidas para o exercício e não validam os prazos de lembrete e alerta do protótipo.
- A curva dos registros supõe que o remédio foi tomado na hora da retirada, mas o sistema registra só a *retirada*, não a ingestão.
- Os resultados não comprovam efeito clínico e não substituem orientação médica.

## 🛠️ Tecnologias utilizadas

- **Python 3:** cálculos e análises.
- **NumPy:** contas com vetores e integração numérica.
- **Matplotlib:** gráficos.
- **Bibliotecas padrão `json`, `math`, `os`, `sys` e `datetime`:** arquivos, contas e datas.

## 👩‍💻 Equipe

| Nome | RM |
|------|----|
| Eric Kang | 572575 |
| Guilherme Tomé Nogueira | 570144 |
| Lucas de Andrade Astorini | 569119 |
| Sabrina Lopes da Silva | 571870 |
| Sofia Satomi Hagio | 569120 |
