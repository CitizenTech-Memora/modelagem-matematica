import json
import math
import os
import sys
from datetime import date, datetime, timedelta

import numpy as np
import matplotlib

MOSTRAR = '--mostrar' in sys.argv
if not MOSTRAR:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

PASTA = os.path.dirname(os.path.abspath(__file__))
PASTA_GRAFICOS = os.path.join(PASTA, 'graficos')
ARQUIVO_DADOS = os.path.join(PASTA, 'dados.json')

#Parâmetros do modelo (fictícios)
T_MEIA_VIDA = 6.0                  #horas: tempo para a concentração cair pela metade
K = math.log(2) / T_MEIA_VIDA      #constante de eliminação (1/h)
C0 = 10.0                          #aumento de concentração causado por UMA dose (u.a.)
TAU = 12.0                         #intervalo entre doses (h): 2 doses por dia
R = math.exp(-K * TAU)             #fração que sobra de uma dose após um intervalo
LIM_MIN = 3.0                      #referência inferior fictícia (u.a.)
LIM_MAX = 16.0                     #referência superior fictícia (u.a.)
C_MEDIA_MIN = 5.0                  #referência fictícia para a média (u.a.)
TOLERANCIA_MINUTOS = 20            #mesma constante do Memora.py (regra RN03)

#Integral numérica (compatível com numpy 1.x e 2.x)
integral_trapezio = getattr(np, 'trapezoid', None) or np.trapz


#Funções do modelo
def concentracao(t, doses):
    """Soma os resíduos das administrações. No salto, usa o lado direito."""
    t = np.asarray(t, dtype=float)
    c = np.zeros_like(t)
    for instante, incremento in doses:
        #Evita calcular exponenciais de doses que ainda não aconteceram.
        decorrido = np.maximum(t - instante, 0.0)
        c += incremento * np.exp(-K * decorrido) * (t >= instante)
    return c


def derivada_exata(t, doses):
    """Entre doses, C'=-kC. Na descontinuidade, retorna NaN."""
    t = np.asarray(t, dtype=float)
    resultado = -K * concentracao(t, doses)
    for instante, incremento in doses:
        if incremento != 0:
            resultado = np.where(np.isclose(t, instante, rtol=0, atol=1e-10),
                                 np.nan, resultado)
    return resultado


def derivada_numerica(t, doses, h=1e-5):
    """Confere a derivada sem atravessar as descontinuidades."""
    t = np.asarray(t, dtype=float)
    resultado = (concentracao(t+h, doses)-concentracao(t-h, doses))/(2*h)
    for instante, incremento in doses:
        if incremento != 0:
            resultado = np.where(np.abs(t-instante) <= h, np.nan, resultado)
    return resultado


def auc_exata(doses, a, b):
    """Integral de C(t) em [a, b], calculada a mão: d/K * (e^{-K(max(a,ti)-ti)} - e^{-K(b-ti)})."""
    total = 0.0
    for ti, di in doses:
        if ti < b:
            inicio = max(a, ti)
            total += di / K * (math.exp(-K * (inicio - ti)) - math.exp(-K * (b - ti)))
    return total


def intervalos(doses, a, b):
    """Separa a janela nos instantes de administração."""
    if b <= a:
        raise ValueError('O fim da janela deve ser posterior ao início.')
    pontos = sorted({a, b, *(ti for ti, _ in doses if a < ti < b)})
    return list(zip(pontos, pontos[1:]))


def auc_numerica(doses, a, b, passo=0.002):
    """Trapézios por intervalo: não desenha uma rampa sobre os saltos."""
    if passo <= 0:
        raise ValueError('O passo deve ser positivo.')
    total = 0.0
    for inicio, fim in intervalos(doses, a, b):
        ativos = [(ti, di) for ti, di in doses if ti <= inicio]
        n = max(1, math.ceil((fim-inicio)/passo))
        tempos = np.linspace(inicio, fim, n+1)
        total += float(integral_trapezio(concentracao(tempos, ativos), tempos))
    return total


def extremos_e_tempos(doses, a, b):
    """Extremos pelos limites laterais e durações por cruzamentos exatos."""
    extremos = []
    abaixo = acima = 0.0
    for inicio, fim in intervalos(doses, a, b):
        c_inicio = float(concentracao(inicio, doses))
        duracao = fim-inicio
        c_fim_antes = c_inicio * math.exp(-K*duracao)
        extremos.extend([c_inicio, c_fim_antes])
        if c_inicio == 0:
            abaixo += duracao
            continue
        tempo_min = math.log(c_inicio/LIM_MIN)/K
        tempo_max = math.log(c_inicio/LIM_MAX)/K
        abaixo += duracao-max(0.0, min(duracao, tempo_min))
        acima += max(0.0, min(duracao, tempo_max))
    #Inclui eventual administração exatamente no fim da janela fechada.
    extremos.append(float(concentracao(b, doses)))
    return max(extremos), min(extremos), abaixo, acima


def curva_com_saltos(doses, a, b):
    """Inclui os dois lados dos saltos; os segmentos verticais não têm duração."""
    tempos = np.unique(np.concatenate((np.linspace(a,b,2500),
                        [ti for ti,_ in doses if a <= ti <= b])))
    x, y = [], []
    for t in tempos:
        depois = float(concentracao(t,doses))
        incremento = sum(di for ti,di in doses if abs(ti-t)<1e-10)
        if incremento:
            x.append(t); y.append(depois-incremento)
        x.append(t); y.append(depois)
    return np.array(x), np.array(y)


def pico_n(n):
    """Pico logo após a n-ésima dose (soma parcial da série geométrica)."""
    return C0 * (1 - R ** n) / (1 - R)


def vale_n(n):
    """Concentração logo antes da (n+1)-ésima dose."""
    return pico_n(n) * R


def patamar_pico(tau=TAU):
    """Limite do pico quando n -> infinito: C0 / (1 - e^{-K tau})."""
    return C0 / (1 - math.exp(-K * tau))


def patamar_vale(tau=TAU):
    return patamar_pico(tau) * math.exp(-K * tau)


def horas_fora_da_janela(doses, a, b):
    _, _, abaixo, acima = extremos_e_tempos(doses, a, b)
    return abaixo, acima


#Cenários (dose no horário, atrasada, esquecida, dobrada)
J_INICIO, J_FIM = 0.0, 72.0   #janela de análise: 3 dias
T_EVENTO = 36.0               #instante da dose que sofre o problema
HISTORICO = [(TAU * i, C0) for i in range(-5, 6)]  #histórico finito: aproxima o regime periódico; não é exatamente o limite


def montar_cenarios():
    base = list(HISTORICO)

    def sem_dose(lista, t):
        return [(ti, di) for ti, di in lista if abs(ti - t) > 1e-9]

    cenarios = {}
    cenarios['No horário'] = base
    cenarios['Atraso de 3 h'] = sem_dose(base, T_EVENTO) + [(T_EVENTO + 3.0, C0)]
    cenarios['Dose esquecida'] = sem_dose(base, T_EVENTO)
    cenarios['Dose dobrada'] = sem_dose(base, T_EVENTO) + [(T_EVENTO, 2 * C0)]
    cenarios['Esquecida e compensada'] = (sem_dose(sem_dose(base, T_EVENTO), T_EVENTO + TAU)
                                          + [(T_EVENTO + TAU, 2 * C0)])
    return cenarios


def analisar_cenarios():
    cenarios = montar_cenarios()
    ref = auc_exata(cenarios['No horário'], J_INICIO, J_FIM)
    tabela = {}
    for nome, doses in cenarios.items():
        auc = auc_exata(doses, J_INICIO, J_FIM)
        pico, vale, abaixo, acima = extremos_e_tempos(doses, J_INICIO, J_FIM)
        tabela[nome] = {
            'auc': auc,
            'auc_numerica': auc_numerica(doses, J_INICIO, J_FIM),
            'delta_auc_pct': 100 * (auc - ref) / ref,
            'pico_max': pico,
            'vale_min': vale,
            'h_abaixo': abaixo,
            'h_acima': acima,
        }
    return cenarios, tabela


#Ligação com o Memora.py (dados.json)
def exemplo_em_memoria():
    """Reproduz o exemplo de criar_exemplo() do Memora.py (opção 8 do menu)."""
    hoje = date.today()
    eventos = []
    for dias_atras in range(1, 7):
        dia = (hoje - timedelta(days=dias_atras)).isoformat()
        for horario in ['08:00', '14:00']:
            if dias_atras == 3 and horario == '14:00':
                tipo, minutos = 'perdida', None
            elif dias_atras == 2 and horario == '14:00':
                tipo, minutos = 'atrasada', 35
            else:
                tipo, minutos = 'no_horario', 5
            eventos.append({'data': dia, 'horario': horario, 'tipo': tipo,
                            'minutos_apos_horario': minutos})
    return 'Remédio A (exemplo)', ['08:00', '14:00'], eventos, 'exemplo em memória'


def carregar_eventos_memora():
    """Lê um medicamento com eventos; o exemplo é usado só quando não há arquivo."""
    if not os.path.exists(ARQUIVO_DADOS):
        nome, horarios, eventos, origem = exemplo_em_memoria()
        # O início corresponde ao primeiro dia do exemplo original.
        return nome, horarios, eventos, origem, (date.today()-timedelta(days=6)).isoformat()
    with open(ARQUIVO_DADOS, encoding='utf-8') as arquivo:
        dados = json.load(arquivo)
    candidatos = []
    for remedio in dados['remedios']:
        eventos = [e for e in dados['eventos'] if e['remedio_id']==remedio['id']
                   and e['tipo'] in ('no_horario','atrasada','perdida')]
        if eventos:
            candidatos.append((len(eventos),remedio,eventos))
    if not candidatos:
        raise ValueError('Nenhum medicamento com registros de retirada foi encontrado.')
    _, remedio, eventos = max(candidatos,key=lambda item:item[0])
    return remedio['nome'],remedio['horarios'],eventos,'dados.json',remedio.get('inicio')


def analisar_dados_memora():
    nome, horarios, eventos, origem, inicio = carregar_eventos_memora()
    #Só dias completos: não interpreta horários de hoje como dose ausente.
    d1 = datetime.combine(date.today(), datetime.min.time())
    d0 = d1-timedelta(days=7)
    if inicio:
        d0 = max(d0, datetime.strptime(inicio,'%Y-%m-%d'))
    if d0 >= d1:
        raise ValueError('Ainda não há dias completos após o início do acompanhamento.')
    horarios = sorted(set(horarios))
    for horario in horarios:
        datetime.strptime(horario,'%H:%M')
    if not horarios:
        raise ValueError('O medicamento não tem horários cadastrados.')
    horas_total = (d1-d0).total_seconds()/3600
    previstos = {}
    dia = d0
    while dia < d1:
        for horario in horarios:
            h,m = map(int,horario.split(':'))
            dt = dia+timedelta(hours=h,minutes=m)
            previstos[(dia.date().isoformat(),horario)] = (dt-d0).total_seconds()/3600
        dia += timedelta(days=1)
    registros = {}
    for evento in eventos:
        chave = (evento['data'],evento['horario'])
        if chave not in previstos:
            continue
        if chave in registros:
            raise ValueError('Há registros duplicados para a mesma dose prevista.')
        atraso = evento.get('minutos_apos_horario')
        if evento['tipo'] in ('no_horario','atrasada'):
            if not isinstance(atraso,(int,float)) or not math.isfinite(atraso) or atraso < 0:
                raise ValueError('Uma retirada tem atraso ausente ou inválido.')
        registros[chave] = evento
    contagem = {tipo:0 for tipo in ('no_horario','atrasada','perdida')}
    registradas, faltas, desconhecidas = [], [], []
    for chave,t in previstos.items():
        evento = registros.get(chave)
        if evento is None:
            desconhecidas.append(t)
        else:
            contagem[evento['tipo']] += 1
            if evento['tipo']=='perdida':
                faltas.append(t)
            else:
                registradas.append((t+evento['minutos_apos_horario']/60,C0))
    ideais = [(t,C0) for t in previstos.values()]
    sem_registro = len(desconhecidas)
    #Só se compara a curva quando não há horários sem informação.
    auc_hipotetica = None if sem_registro else auc_exata(registradas,0,horas_total)
    auc_ideal = auc_exata(ideais,0,horas_total)
    return {'remedio':nome,'origem':origem,'horas':horas_total,'previstas':len(ideais),
            'contagem':contagem,'sem_registro':sem_registro,'ideais':ideais,
            'registradas':registradas,'faltas':faltas,'desconhecidas':desconhecidas,
            'taxa_memora':100*contagem['no_horario']/len(ideais),
            'auc_ideal':auc_ideal,'auc_hipotetica':auc_hipotetica,
            'indice_exposicao':None if auc_hipotetica is None else 100*auc_hipotetica/auc_ideal,
            'data_inicio':d0,'data_fim':d1,'status':'ok'}


#Cálculo de todos os resultados (usado também no relatório)
def calcular():
    res = {}
    #Parte 1 - derivada
    res['k'] = K
    res['r'] = R
    res['taxa_inicial'] = -K * C0
    res['t_ate_cair_ao_vale'] = TAU
    t_teste = 5.0
    res['deriv_exata'] = float(derivada_exata(np.array([t_teste]), [(0.0, C0)])[0])
    res['deriv_numerica'] = float(derivada_numerica(np.array([t_teste]), [(0.0, C0)])[0])
    vale_ss = patamar_vale()
    res['vale_ss'] = vale_ss
    res['pico_ss'] = patamar_pico()
    linha = []
    for atraso_min in (10, 20, 30, 60, 180):
        d = atraso_min / 60
        exata = vale_ss * (1 - math.exp(-K * d))
        linear = vale_ss * K * d
        linha.append({'min': atraso_min, 'exata': exata, 'linear': linear,
                      'queda_pct': 100 * exata / vale_ss,
                      'erro_rel_pct': 100 * abs(exata - linear) / exata})
    res['diferencial'] = linha
    res['atraso_5pct_min'] = 60 * (-math.log(0.95) / K)
    res['atraso_ate_limite_min'] = 60 * (math.log(vale_ss / LIM_MIN) / K)
    res['intervalo_dose_extra_h'] = math.log(patamar_pico() / (LIM_MAX - C0)) / K
    #Parte 2 - limites
    res['serie'] = [(n, pico_n(n), vale_n(n)) for n in range(1, 9)]
    res['n_99'] = math.ceil(math.log(0.01) / math.log(R))
    res['dobro_para_sempre'] = 2 * patamar_pico()
    res['peak_vs_tau'] = [(tau, patamar_pico(tau)) for tau in (1, 2, 4, 6, 8, 12, 24)]
    #Parte 3 - integrais
    res['auc_dose_unica'] = C0 / K
    res['auc_intervalo_ss'] = C0 / K
    res['auc_intervalo_1dose'] = C0 / K * (1 - R)
    res['c_media_ss'] = C0 / (K * TAU)
    res['adesao_min_pct'] = 100 * C_MEDIA_MIN / res['c_media_ss']
    cenarios, tabela = analisar_cenarios()
    res['cenarios'] = cenarios
    res['tabela'] = tabela
    try:
        res['memora'] = analisar_dados_memora()
    except (OSError, ValueError, KeyError, TypeError) as erro:
        res['memora'] = {'status':'indisponivel', 'motivo':str(erro)}
    return res


#no terminal
def imprimir(res):
    print('=' * 66)
    print(' MEMORA - MODELAGEM MATEMÁTICA: POR QUE O HORÁRIO IMPORTA')
    print('=' * 66)
    print(f'Meia-vida = {T_MEIA_VIDA:.0f} h | k = ln2/t1/2 = {K:.5f} 1/h | C0 = {C0:.0f} u.a.')
    print(f'Intervalo entre doses = {TAU:.0f} h | r = e^(-k*tau) = {R:.4f}')
    print(f'Referências fictícias: {LIM_MIN:.0f} a {LIM_MAX:.0f} u.a. (valores fictícios)')

    print('\n--- PARTE 1: DERIVADA (taxa de variação) ---')
    print(f"Derivada à direita em t=0+ (uma dose): {res['taxa_inicial']:.4f} u.a./h")
    print(f"Conferência em t=5 h: exata = {res['deriv_exata']:.6f} | numérica = {res['deriv_numerica']:.6f}")
    print(f"Vale em regime: {res['vale_ss']:.3f} u.a. | pico em regime: {res['pico_ss']:.3f} u.a.")
    print('Atraso | queda exata do vale | aproximação C\'(t)*dt | erro')
    for l in res['diferencial']:
        print(f"{l['min']:4d} min | {l['exata']:.4f} u.a. ({l['queda_pct']:.1f}%) | "
              f"{l['linear']:.4f} u.a. | {l['erro_rel_pct']:.2f}%")
    print(f"Atraso que derruba o vale em 5%: {res['atraso_5pct_min']:.1f} min")
    print(f"Atraso que leva o vale à referência inferior: {res['atraso_ate_limite_min']:.1f} min "
          f"(tolerância do Memora.py: {TOLERANCIA_MINUTOS} min)")
    print(f"Cruzamento da referência superior para incremento extra: {res['intervalo_dose_extra_h']:.1f} h após o pico limite (hipótese fictícia)")

    print('\n--- PARTE 2: LIMITES ---')
    print('n | pico após n-ésima dose | vale antes da próxima')
    for n, p, v in res['serie']:
        print(f'{n} | {p:8.3f} | {v:8.3f}')
    print(f"lim pico = C0/(1-r) = {res['pico_ss']:.3f} | lim vale = {res['vale_ss']:.3f}")
    print(f"Doses para chegar a 99% do patamar: {res['n_99']}")
    print(f"Se TODA dose fosse dobrada, o patamar iria a {res['dobro_para_sempre']:.2f} (> {LIM_MAX:.0f})")
    print('Pico em regime x intervalo (tau -> 0 faz o pico -> infinito):')
    for tau, p in res['peak_vs_tau']:
        print(f'  tau = {tau:2d} h -> {p:7.2f}')

    print('\n--- PARTE 3: INTEGRAIS (exposição acumulada) ---')
    print(f"Uma dose, de 0 a infinito: C0/k = {res['auc_dose_unica']:.2f} u.a.*h")
    print(f"Concentração média em regime: {res['c_media_ss']:.2f} u.a. | fração administrada teórica p/ média >= "
          f"{C_MEDIA_MIN:.0f}: {res['adesao_min_pct']:.1f}%")
    print(f'Janela de {J_INICIO:.0f} a {J_FIM:.0f} h, problema na dose de t = {T_EVENTO:.0f} h:')
    print(f"{'Cenário':26s} {'AUC':>8s} {'AUC num':>8s} {'dAUC%':>7s} {'Pico':>6s} {'Vale':>6s} "
          f"{'h<min':>6s} {'h>max':>6s}")
    for nome, v in res['tabela'].items():
        print(f"{nome:26s} {v['auc']:8.2f} {v['auc_numerica']:8.2f} {v['delta_auc_pct']:7.2f} "
              f"{v['pico_max']:6.2f} {v['vale_min']:6.2f} {v['h_abaixo']:6.2f} {v['h_acima']:6.2f}")

    m = res['memora']
    print('\n--- PARTE 4: REGISTROS DO MEMORA ---')
    if m['status'] != 'ok':
        print('Análise de registros não realizada:', m['motivo'])
        return
    print(f"Remédio: {m['remedio']} | período: {m['horas']:.0f} h | doses previstas: {m['previstas']}")
    print(f"No horário: {m['contagem']['no_horario']} | atrasadas: {m['contagem']['atrasada']} | "
          f"perdidas: {m['contagem']['perdida']}")
    print(f"Taxa de adesão do Memora.py: {m['taxa_memora']:.1f}%")
    print(f"Doses sem informação: {m['sem_registro']} | origem: {m['origem']}")
    if m['indice_exposicao'] is None:
        print('Curva hipotética não calculada: há doses sem informação.')
    else:
        print(f"AUC hipotética / AUC de referência: {m['indice_exposicao']:.1f}%")
    print('A curva pressupõe ingestão na retirada; o sistema não confirma essa hipótese.')


#gráficos
COR = {'ref': '#1f77b4', 'atraso': '#ff7f0e', 'esq': '#d62728', 'dobro': '#9467bd',
       'comp': '#2ca02c', 'faixa': '#2e7d32', 'risco': '#c62828', 'cinza': '#666666'}
COR_CEN = {'No horário': COR['ref'], 'Atraso de 3 h': COR['atraso'], 'Dose esquecida': COR['esq'],
           'Dose dobrada': COR['dobro'], 'Esquecida e compensada': COR['comp']}


def salvar(fig, nome):
    os.makedirs(PASTA_GRAFICOS, exist_ok=True)
    fig.savefig(os.path.join(PASTA_GRAFICOS, nome), dpi=170, bbox_inches='tight')
    if not MOSTRAR:
        plt.close(fig)


def faixa_alvo(ax, x0, x1):
    ax.axhspan(LIM_MIN, LIM_MAX, color=COR['faixa'], alpha=0.08, label='Referências fictícias')
    ax.axhline(LIM_MIN, color=COR['faixa'], lw=1, ls='--')
    ax.axhline(LIM_MAX, color=COR['risco'], lw=1, ls='--')
    ax.set_xlim(x0, x1)


def grafico_derivada(res):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    t = np.linspace(0, 24, 600)
    c = C0 * np.exp(-K * t)
    ax[0].plot(t, c, color=COR['ref'], lw=2.2, label='C(t) = C0·e^(-kt)')
    ax[0].plot(t, C0 - K * C0 * t, color=COR['atraso'], ls='--', lw=1.4, label="Tangente em t = 0 (C'(0+) = -kC0)")
    for n in (1, 2, 3):
        ax[0].plot(n * T_MEIA_VIDA, C0 / 2 ** n, 'o', color=COR['esq'])
        ax[0].annotate(f'{C0 / 2 ** n:.2f}', (n * T_MEIA_VIDA, C0 / 2 ** n), textcoords='offset points',
                       xytext=(6, 6), fontsize=8)
    ax[0].set_ylim(0, C0 * 1.05)
    ax[0].set_title('Uma dose: a concentração cai pela metade a cada 6 h')
    ax[0].set_xlabel('Tempo após a dose (h)')
    ax[0].set_ylabel('Concentração (u.a.)')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    ax[1].plot(t, -K * c, color=COR['esq'], lw=2.2, label="C'(t) = -k·C(t)")
    ax[1].axhline(0, color='k', lw=0.8)
    ax[1].set_title('Taxa de variação: elimina mais rápido quando há mais remédio')
    ax[1].set_xlabel('Tempo após a dose (h)')
    ax[1].set_ylabel("dC/dt (u.a./h)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    salvar(fig, 'fig1_derivada.png')


def grafico_limite(res):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ns = np.arange(1, 11)
    ax[0].plot(ns, [pico_n(n) for n in ns], 'o-', color=COR['ref'], label='Pico após a n-ésima dose')
    ax[0].plot(ns, [vale_n(n) for n in ns], 's-', color=COR['atraso'], label='Vale antes da próxima')
    ax[0].axhline(res['pico_ss'], color=COR['ref'], ls=':', label=f"Limite do pico = {res['pico_ss']:.2f}")
    ax[0].axhline(res['vale_ss'], color=COR['atraso'], ls=':', label=f"Limite do vale = {res['vale_ss']:.2f}")
    ax[0].set_title('Doses repetidas: picos e vales convergem')
    ax[0].set_xlabel('Número de doses n')
    ax[0].set_ylabel('Concentração (u.a.)')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    taus = np.linspace(1, 24, 200)
    ax[1].plot(taus, [patamar_pico(x) for x in taus], color=COR['dobro'], lw=2)
    ax[1].axvline(TAU, color=COR['cinza'], ls='--')
    ax[1].annotate('Modelo (12 h)', (TAU, res['pico_ss']), textcoords='offset points', xytext=(6, 14), fontsize=8)
    ax[1].axhline(LIM_MAX, color=COR['risco'], ls='--', lw=1, label='Referência superior')
    ax[1].set_ylim(0, 50)
    ax[1].set_title('Pico limite quando τ diminui (τ → 0)')
    ax[1].set_xlabel('Intervalo entre doses τ (h)')
    ax[1].set_ylabel('Pico em regime (u.a.)')
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    salvar(fig, 'fig2_limite.png')


def grafico_cenarios(res):
    cen, tab = res['cenarios'], res['tabela']
    t = np.linspace(J_INICIO, J_FIM, 8000)
    ref = concentracao(t, cen['No horário'])
    fig, axs = plt.subplots(3, 2, figsize=(11, 10), sharex=True, sharey=True)
    axs = axs.ravel()
    for ax, (nome, doses) in zip(axs, cen.items()):
        t, c = curva_com_saltos(doses, J_INICIO, J_FIM)
        ref = concentracao(t, cen['No horário'])
        faixa_alvo(ax, J_INICIO, J_FIM)
        ax.fill_between(t, 0, c, color=COR_CEN[nome], alpha=0.22)
        ax.plot(t, ref, color=COR['cinza'], lw=0.9, ls=':')
        ax.plot(t, c, color=COR_CEN[nome], lw=1.8)
        ax.set_title(f"{nome}  |  área = {tab[nome]['auc']:.0f} ({tab[nome]['delta_auc_pct']:+.1f}%)", fontsize=10)
        ax.grid(alpha=0.3)
        ax.set_ylim(0, 27)
    ax = axs[5]
    faixa_alvo(ax, J_INICIO, J_FIM)
    for nome, doses in cen.items():
        x,y = curva_com_saltos(doses, J_INICIO, J_FIM)
        ax.plot(x,y, color=COR_CEN[nome], lw=1.4, label=nome)
    ax.set_title('Todos os cenários juntos', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(alpha=0.3)
    for ax in axs[4:]:
        ax.set_xlabel('Tempo (h)')
    for ax in axs[::2]:
        ax.set_ylabel('Concentração (u.a.)')
    fig.suptitle('Efeito de uma única dose problemática (t = 36 h) sobre C(t) e sobre a área sob a curva', y=0.995)
    fig.tight_layout()
    salvar(fig, 'fig3_cenarios.png')


def grafico_comparacao(res):
    tab = res['tabela']
    nomes = list(tab.keys())
    cores = [COR_CEN[n] for n in nomes]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].barh(nomes, [tab[n]['auc'] for n in nomes], color=cores)
    for i, n in enumerate(nomes):
        ax[0].text(tab[n]['auc'] + 4, i, f"{tab[n]['auc']:.0f}  ({tab[n]['delta_auc_pct']:+.1f}%)",
                   va='center', fontsize=8)
    ax[0].set_xlim(0, max(tab[n]['auc'] for n in nomes) * 1.3)
    ax[0].invert_yaxis()
    ax[0].set_title('Exposição acumulada (integral) em 72 h')
    ax[0].set_xlabel('Área sob a curva (u.a.·h)')
    y = np.arange(len(nomes))
    ax[1].barh(y - 0.2, [tab[n]['h_abaixo'] for n in nomes], height=0.4, color=COR['atraso'], label='Abaixo do mínimo')
    ax[1].barh(y + 0.2, [tab[n]['h_acima'] for n in nomes], height=0.4, color=COR['risco'], label='Acima do máximo')
    ax[1].set_yticks(y)
    ax[1].set_yticklabels(nomes)
    ax[1].invert_yaxis()
    ax[1].set_title('Horas fora da faixa de referência')
    ax[1].set_xlabel('Horas (em 72 h)')
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    salvar(fig, 'fig4_comparacao.png')


def grafico_tolerancia(res):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    atrasos = np.linspace(0, 180, 400)
    queda = 100 * (1 - np.exp(-K * atrasos / 60))
    linear = 100 * K * atrasos / 60
    ax[0].plot(atrasos, queda, color=COR['ref'], lw=2.2, label='Queda exata do vale')
    ax[0].plot(atrasos, linear, color=COR['atraso'], ls='--', label="Aproximação linear C'(t)·Δt")
    ax[0].axvline(TOLERANCIA_MINUTOS, color=COR['faixa'], lw=1.6, label=f'Tolerância do Memora ({TOLERANCIA_MINUTOS} min)')
    ax[0].axvline(res['atraso_ate_limite_min'], color=COR['risco'], ls=':', lw=1.6,
                  label=f"Vale chega ao mínimo ({res['atraso_ate_limite_min']:.0f} min)")
    ax[0].set_title('Quanto o vale cai a cada minuto de atraso')
    ax[0].set_xlabel('Atraso da dose (min)')
    ax[0].set_ylabel('Queda do vale (%)')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    p = np.linspace(0, 100, 200)
    ax[1].plot(p, res['c_media_ss'] * p / 100, color=COR['dobro'], lw=2.2)
    ax[1].axhline(C_MEDIA_MIN, color=COR['risco'], ls='--', lw=1, label=f'Média mínima ({C_MEDIA_MIN:.0f} u.a.)')
    ax[1].axvline(res['adesao_min_pct'], color=COR['faixa'], ls='--', lw=1,
                  label=f"Fração teórica ≈ {res['adesao_min_pct']:.0f}%")
    ax[1].set_title('Média teórica × fração administrada')
    ax[1].set_xlabel('Fração administrada hipotética (%)')
    ax[1].set_ylabel('Concentração média (u.a.)')
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    salvar(fig, 'fig5_tolerancia_adesao.png')


def grafico_memora(res):
    m = res['memora']
    if m['status'] != 'ok':
        return
    fig,ax = plt.subplots(figsize=(11,4))
    if m['sem_registro']:
        ax.bar(['No horário','Atrasadas','Perdidas','Sem informação'],
               [m['contagem']['no_horario'],m['contagem']['atrasada'],
                m['contagem']['perdida'],m['sem_registro']],
               color=['#1f77b4','#ff7f0e','#d62728','#999999'])
        ax.set_ylabel('Número de registros / horários')
        ax.set_title('Registros incompletos: a concentração não foi simulada')
    else:
        for doses,rotulo,cor,estilo in [
            (m['ideais'],'Referência: todas no horário','#666666',':'),
            (m['registradas'],'Hipótese: ingestão no instante da retirada','#1f77b4','-')]:
            t,c = curva_com_saltos(doses,0,m['horas'])
            ax.plot(t,c,label=rotulo,color=cor,ls=estilo)
        for t in m['faltas']:
            ax.axvline(t,color='#d62728',ls='--',lw=1)
        ax.set(xlabel='Horas desde o início da janela',ylabel='Concentração hipotética (u.a.)',
               title='Exercício com registros de retirada; não é medição no usuário')
        ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    fig.tight_layout()
    salvar(fig,'fig6_dados_memora.png')


def verificar_calculos(res):
    """Confere propriedades independentes das tabelas impressas."""
    assert math.isclose(float(concentracao(6,[(0,C0)])),C0/2)
    assert np.isnan(derivada_exata(36,res['cenarios']['No horário']))
    assert math.isclose(res['deriv_exata'],res['deriv_numerica'],rel_tol=1e-8)
    for n in range(1,8):
        assert math.isclose(pico_n(n+1),C0+R*pico_n(n))
    for v in res['tabela'].values():
        assert math.isclose(v['auc'],v['auc_numerica'],rel_tol=1e-7)
        assert 0 <= v['h_abaixo']+v['h_acima'] <= J_FIM-J_INICIO
    tab = res['tabela']
    assert math.isclose(tab['Dose dobrada']['auc']-tab['No horário']['auc'],
                        tab['No horário']['auc']-tab['Dose esquecida']['auc'])
    assert math.isclose(extremos_e_tempos([(0,10)],0,6)[0],10)
    assert math.isclose(extremos_e_tempos([(0,10)],0,6)[1],5)


def main():
    res = calcular()
    verificar_calculos(res)
    imprimir(res)
    grafico_derivada(res)
    grafico_limite(res)
    grafico_cenarios(res)
    grafico_comparacao(res)
    grafico_tolerancia(res)
    grafico_memora(res)
    with open(os.path.join(PASTA,'resultados.json'),'w',encoding='utf-8') as arquivo:
        json.dump(res,arquivo,ensure_ascii=False,indent=2,default=str,allow_nan=False)
    print('Verificações matemáticas concluídas.')
    print(f'\nGráficos salvos em: {PASTA_GRAFICOS}')
    if MOSTRAR:
        plt.show()


if __name__ == '__main__':
    main()
