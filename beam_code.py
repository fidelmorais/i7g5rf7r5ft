
import math

def gammas(combinacao):
    if combinacao == 'normais':
        gamma_c = 1.4
        gamma_s = 1.15
    elif combinacao == 'especiais':
        gamma_c = 1.2
        gamma_s = 1.15
    elif combinacao == 'excepcionais':
        gamma_c = 1.2
        gamma_s = 1.0
    else:
        raise ValueError("Combinacão de carregamento desconhecida. Use 'normais', 'especiais' ou 'excepcionais'.")

    return gamma_c, gamma_s

def coefs(f_ck):
    """
    Calcula os coeficientes λ, α_c, η_c, μ_lim e ω_lim com base no valor de f_ck.

    :param f_ck: Resistência característica à compressão do concreto (em MPa).
    :return: Tuple contendo os valores de λ, α_c, η_c, μ_lim e ω_lim.
    """

    # Cálculo de η_c
    if f_ck <= 40:
        η_c = 1
    elif 40 < f_ck <= 50:
        η_c = (40 / f_ck) ** (1 / 3)
    else:
        η_c = 1.4 + 23.4 * ((90 - f_ck) / 100) ** 4

    # Cálculo de α_c
    if f_ck <= 50:
        α_c = 0.85
    else:
        α_c = 0.85 * (1 - (f_ck - 50) / 200)

    # Cálculo de λ
    if f_ck <= 50:
        λ = 0.8
    else:
        λ = 0.8 - (f_ck - 50) / 400

    # Tabela de valores de "mi" e "ômega" para f_ck entre 50 e 90 MPa
    tabelas = {
        50: (0.2952, 0.3600),
        55: (0.2376, 0.2756),
        60: (0.2345, 0.2713),
        65: (0.2313, 0.2669),
        70: (0.2280, 0.2625),
        75: (0.2248, 0.2581),
        80: (0.2216, 0.2538),
        85: (0.2183, 0.2494),
        90: (0.2150, 0.2450),
    }

    # Cálculo de μ_lim e ω_lim
    if f_ck <= 50:
        μ_lim, ω_lim = tabelas[50]
    elif f_ck in tabelas:
        μ_lim, ω_lim = tabelas[f_ck]
    else:
        raise ValueError("f_ck fora do intervalo permitido. Deve ser ≤ 90 MPa.")

    return λ, α_c, η_c, μ_lim, ω_lim

def calcular_fcd_fyd(f_ck, f_yk, combinacao):
    λ, α_c, η_c, μ_lim, ω_lim = coefs(f_ck)
    gamma_c, gamma_s = gammas(combinacao)

    f_cd = f_ck / gamma_c
    f_yd = f_yk / gamma_s

    return f_cd, f_yd, η_c, λ, α_c, μ_lim, ω_lim

def calcular_area_aco_ret(M_sd, f_ck, f_yk, b_w, d, d_prime, E_s, combinacao):
    f_cd, f_yd, η_c, λ, α_c, μ_lim, ω_lim = calcular_fcd_fyd(f_ck, f_yk, combinacao)

    M_sd = M_sd * 100   # converte o momento dado em kn.m para kn.cm
    f_cd = f_cd/10      # converte a tensão dada em MPa para Kn/cm²
    f_yd = f_yd/10      # converte a tensão dada em MPa para Kn/cm²

    mu = M_sd / (α_c * η_c * f_cd * b_w * (d**2))

    if mu > μ_lim:
        print(f"Atenção: µ ({mu:.4f}) excede o µ_lim ({μ_lim:.4f}). Usando armadura dupla.")
        Msd_lim = μ_lim * α_c * η_c * f_cd * b_w * (d**2)
        As1 = Msd_lim / (f_yd * d * (1 - ω_lim / 2))
        delta_M_sd = M_sd - Msd_lim
        As2 = delta_M_sd / (f_yd * (d - d_prime))

        # Cálculo da deformação e tensão da armadura comprimida
        ξ_lim = ω_lim / λ
        # ε_s_prime = 0.0035 * (1 - d_prime / (ξ_lim * d))
        ε_s_prime = (0.0035 * ((ξ_lim*d) - d_prime)) / (ξ_lim*d)
        ε_yd = f_yk / E_s

        if ε_s_prime >= ε_yd:
            σ_s_prime = f_yd
        else:
            σ_s_prime = ε_s_prime * E_s

        A_s_prime = delta_M_sd / ((σ_s_prime - α_c * η_c * f_cd) * (d - d_prime))
        As_total = As1 + As2
        return As_total, A_s_prime

    else:
        omega = 1 - math.sqrt(1 - 2 * mu)
        As = (omega * α_c * η_c * f_cd * b_w * d) / f_yd
        return As, 0

def verific_flex(b, h, f_ck, M_sd, A_s, A_s_prime,omega,x, λ):
    """
    Função para verificar a flexão de uma viga de concreto armado, dada a área de aço calculada.

    :param b: Largura da seção transversal da viga (em cm).
    :param h: Altura da seção transversal da viga (em cm).
    :param f_ck: Resistência característica à compressão do concreto (em MPa).
    :param M_sd: Momento solicitante (em kNm).
    :param A_s: Área de aço tracionada calculada (em cm²).
    :param A_s_prime: Área de aço comprimida calculada (em cm²), se houver.
    :return: Resultado da verificação de flexão.
    """

    # Dimensões e resistência do material
    A_c = b * h  # Área da seção de concreto (cm²)
    A_s_min = 0.0015 * A_c  # 0,15% da área da seção
    A_s_max = 0.04 * A_c    # 4% da área da seção

    # Cálculo de f_ctm e f_ctk,sup
    f_ctm = 0.3 * (f_ck ** (2 / 3))
    f_ctk_sup = 1.3 * f_ctm
    f_ctk_sup = f_ctk_sup/10 #convertendo para KN/cm^2

    # Verificação do momento resistente mínimo
    Yt = d - (omega * d) / (λ)
    W_0 = ((b * (h ** 2)) / 12)/Yt  # Momento de inércia dividido pela distância do eixo neutro
    M_sd_min = 0.8 * W_0 * f_ctk_sup  # Mínimo momento resistente

    # Verificação de momento resistente
    if M_sd >= M_sd_min:
        print(f"M_sd ({M_sd:.2f} kNm) ≥ M_sd,min ({M_sd_min:.2f} kNm): Verificação de flexão OK.")
    else:
        print(f"M_sd ({M_sd:.2f} kNm) < M_sd,min ({M_sd_min:.2f} kNm): Flexão insuficiente.")

    # Verificação dos limites de armadura
    if A_s_min <= A_s <= A_s_max:
        print(f"A_s ({A_s:.2f} cm²) está dentro dos limites permitidos ({A_s_min:.2f} cm² - {A_s_max:.2f} cm²).")
    else:
        print(f"A_s ({A_s:.2f} cm²) está fora dos limites permitidos ({A_s_min:.2f} cm² - {A_s_max:.2f} cm²).")

    # Se houver área de aço comprimida (A_s_prime > 0), verifica se é razoável
    if A_s_prime > 0:
        print(f"A_s' (área de aço comprimida) = {A_s_prime:.2f} cm².")
    else:
        print("Sem armadura comprimida (A_s' = 0 cm²).")

    return {
        "M_sd_min": M_sd_min,
        "A_s_min": A_s_min,
        "A_s_max": A_s_max,
        "verificacao_momento": M_sd >= M_sd_min,
        "verificacao_armadura": A_s_min <= A_s <= A_s_max
    }



# Calcular as áreas de aço sem dividir f_ck, f_yk e E_s
As_total, A_s_prime = calcular_area_aco_ret(M_sd, f_ck, f_yk, b_w, d, d_prime, E_s, combinacao)

print(f'Área de aço tracionada (As_total) = {As_total:.2f} cm²')
print(f'Área de aço comprimida (A\'s) = {A_s_prime:.2f} cm²')

import numpy as np
import math

def calcular_area_aco_viga_T_2(M_sd, f_ck, f_yk, b_w, d, bf, hf,E_s,d_prime, combinacao):

    """
    Calcula a área de aço As necessária para uma viga em T.

    :param M_sd: Momento fletor de projeto (em kNm).
    :param f_ck: Resistência característica à compressão do concreto (em MPa).
    :param f_yk: Resistência característica ao escoamento do aço (em MPa).
    :param b_w: Largura da alma da viga (em cm).
    :param d: Altura útil da viga (em cm).
    :param bf: Largura da mesa da viga T (em cm).
    :param hf: Altura da mesa da viga T (em cm).
    :param combinacao: Tipo de combinação de carregamento. Pode ser 'normais', 'especiais', ou 'excepcionais'.
    :return: Área de aço necessária (As_total) e (As2) em cm².
    """

    # Calcular os parâmetros de material e seção
    f_cd, f_yd, η_c, λ, α_c, μ_lim, ω_lim = calcular_fcd_fyd(f_ck, f_yk, combinacao)

    f_cd = f_cd/10      # converte a tensão dada em MPa para Kn/cm²
    f_yd = f_yd/10      # converte a tensão dada em MPa para Kn/cm²

    mu = (M_sd * 100) / (α_c * η_c * f_cd * bf * (d**2))   # converte o momento dado em kn.m para kn.cm
    omega = 1 - ((1 - 2 * mu)**(1/2))
    x = (omega * d) / (λ)                                   # calculando a profundidade da linha neutra.

    if x <= hf:
      return calcular_area_aco_ret(M_sd, f_ck, f_yk, bf, d, d_prime, E_s, combinacao)

    else:

      # Cálculo da força de compressão resistida pelas abas
      F_c1 = α_c * η_c * f_cd * hf * (bf - b_w)

      # Cálculo do braço de alavanca z1
      z1 = d - (hf / 2)

      # Cálculo do momento resistente das abas
      M_sd1 = F_c1 * z1

      # Cálculo da área de aço As1 (contribuição das abas)
      As1 = M_sd1 / (z1 * f_yd)

      # Momento restante a ser resistido pela alma
      M_sd2 = M_sd*100 - M_sd1

      # Cálculo da área de aço As2 (contribuição da alma) e As3 (armadura de compressão necessária para a nervura se for o caso).
      As2,As3 = calcular_area_aco_ret((M_sd2/100), f_ck, f_yk, b_w, d, d_prime, E_s, combinacao)

      # Área de aço total
      As_total = As1 + As2
      As_prime = As3

      return As_total, As_prime


import math
import numpy as np
import pandas as pd

def verif_compr_diag(Vk, bw, d, f_ck, f_yk, combinacao):
    # Calcular os parâmetros de material e seção
    f_cd, f_yd, η_c, λ, α_c, μ_lim, ω_lim = calcular_fcd_fyd(f_ck, f_yk, combinacao)

    f_cd = f_cd / 10  # converte a tensão dada em MPa para kN/cm²
    f_yd = f_yd / 10  # converte a tensão dada em MPa para kN/cm²
    fctm = 0.3 * (f_ck ** (2 / 3))
    fctk_inf = 0.7 * fctm
    gamma_c, gamma_s = gammas(combinacao)
    fctd = fctk_inf / gamma_c  # Considerando gamma_c para a resistência à tração do concreto
    alpha_v2 = 1 - (f_ck / 250)
    vrd2 = 0.27 * alpha_v2 * f_cd * bw * d
    Vd = Vk * 1.4  # Majoração do Vk

    status_compr = 'ok' if Vd <= vrd2 else 'reprovado'

    return {
        'alpha_v2': alpha_v2,
        'Vrd2 (kN)': vrd2,
        'Vd (kN)': Vd,
        'Status': status_compr
    }

def verif_trac_diag(Vk, bw, d, f_ck, f_yk, combinacao):
    # Calcular os parâmetros de material e seção
    f_cd, f_yd, η_c, λ, α_c, μ_lim, ω_lim = calcular_fcd_fyd(f_ck, f_yk, combinacao)

    f_cd = f_cd / 10  # converte a tensão dada em MPa para kN/cm²
    f_yd = f_yd / 10  # converte a tensão dada em MPa para kN/cm²
    fctm = 0.3 * (f_ck ** (2 / 3))
    fctk_inf = 0.7 * fctm
    gamma_c, gamma_s = gammas(combinacao)
    fctd = fctk_inf / gamma_c  # Considerando gamma_c para a resistência à tração do concreto
    Vd = Vk * 1.4  # Majoração do Vk

    Vc = 0.6 * fctd * bw * d
    asw_min_cm = 0.2 * fctm * (bw / f_yk)
    asw_min_m = asw_min_cm * 100
    Vsw_min = asw_min_cm * 0.9 * d * (f_yd)
    vrd3_min = Vsw_min + Vc

    if Vd <= vrd3_min:
        asw_adot = asw_min_m
        status_trac = 'Asw mínima'
    else:
        asw_adot = 100 * (Vd - Vc) / (0.9 * d * f_yd)
        status_trac = 'Asw acima da mínima'

    return {
        'Vc (kN)': Vc,
        'asw_min (cm²/cm)': asw_min_cm,
        'asw_min (cm²/m)': asw_min_m,
        'Vrd3_min (kN)': vrd3_min,
        'Status': status_trac,
        'asw_adot (cm²/m)': asw_adot
    }

def detalhar_asw(asw_adot, stirrup_leg):
    diametros = np.array([5.0, 6.3, 8.0, 10.0, 12.5, 16.0, 20.0, 25.0, 32.0, 40.0])
    areas = stirrup_leg * np.pi * (diametros / 10) ** 2 / 4
    espacamento = (areas / asw_adot) * 100
    espacamento_arred = np.floor(espacamento)

    return pd.DataFrame({
        'Diâmetro (mm)': diametros,
        'Espaçamento (cm)': espacamento_arred
    })

def calcular_viga_cisalhamento(nome, bw, h, Vk, f_ck, f_yk, stirrup_leg, combinacao):
    d = h - 5  # Altura útil da viga
    props_concreto = calcular_fcd_fyd(f_ck, f_yk, combinacao)
    compr_diag = verif_compr_diag(Vk, bw, d, f_ck, f_yk, combinacao)
    trac_diag = verif_trac_diag(Vk, bw, d, f_ck, f_yk, combinacao)
    detalhamento = detalhar_asw(trac_diag['asw_adot (cm²/m)'], stirrup_leg)

    print(f'RESULTADOS - {nome} ({bw} X {h} cm)')
    print('')
    print('PROPRIEDADES DO CONCRETO:')
    print(props_concreto)
    print('')
    print('VERIFICAÇÃO DA DIAGONAL COMPRIMIDA:')
    print(compr_diag)
    print('')
    print('VERIFICAÇÃO DA BIELA TRACIONADA:')
    print(trac_diag)
    print('')
    print('DETALHAMENTO DOS ESTRIBOS:')
    print(detalhamento)


import math
import numpy as np
import pandas as pd

# Função que calcula a quantidade de barras necessárias para uma área de aço
def barrasAco(areaAco):
    # Array de bitolas disponíveis
    bitola = np.array([5.0, 6.3, 8.0, 10.0, 12.5, 16.0, 20.0, 22.0, 25.0, 32.0, 40.0])
    # Calcula a área de cada bitola em cm^2
    areaBitola = np.pi * (bitola / 10) ** 2 / 4
    # Calcula a quantidade de barras necessária para atingir a área de aço
    quantBarras = areaAco / areaBitola
    # Arredonda para o próximo número inteiro de barras
    quantBarrasRound = np.ceil(quantBarras)
    # Exporta o resultado em um Series do pandas com bitolas como índice
    exportResult = pd.Series(quantBarrasRound, index=bitola)

    return exportResult


import numpy as np

# Função que calcula o espaçamento entre as barras de aço de acordo com a bitola e espaço disponível
def alojVigas(bw, cob, diamEstribo, diamAgreg, bitola):
    """
    Calcula o espaçamento entre barras de aço em uma viga.

    :param bw: Largura da viga (cm).
    :param cob: Cobrimento nominal (cm).
    :param diamEstribo: Diâmetro do estribo (mm).
    :param diamAgreg: Diâmetro máximo do agregado (mm).
    :param bitola: Diâmetro das barras de aço (mm).
    :return: Lista com espaçamentos viáveis entre as barras.
    """

    # Calcula a largura disponível para alocar as barras
    dispAloj = bw - (2 * cob) - (2 * diamEstribo / 10)  # Converte diâmetro do estribo de mm para cm

    # Define um intervalo de possíveis quantidades de barras (de 2 a 100 barras)
    quantBarras = np.linspace(2, 100, 99)

    # Calcula o espaçamento entre as barras
    ah = (dispAloj - (quantBarras * bitola / 10)) / (quantBarras - 1)

    # Define as condições mínimas de espaçamento
    condicoes = [2.0, bitola / 10, 1.2 * diamAgreg / 10]  # Converte bitola e agregado para cm
    valorLim = max(condicoes)  # Usa o valor máximo como limitante

    # Seleciona apenas os espaçamentos que atendem à condição mínima
    ahSelect = ah[ah >= valorLim]

    # Exibe os resultados de bitolas e espaçamentos possíveis
    for i in range(len(ahSelect)):
        print(f'{i + 2} \u03C6 {bitola} -> espaçamento = {ahSelect[i]:.2f} cm')

    return ahSelect

# Função que executa o cálculo de espaçamento para várias bitolas
def opcoesAloj(bw, cob, diamEstribo, diamAgreg):
    """
    Exibe opções de espaçamento para diferentes bitolas.

    :param bw: Largura da viga (cm).
    :param cob: Cobrimento nominal (cm).
    :param diamEstribo: Diâmetro do estribo (mm).
    :param diamAgreg: Diâmetro máximo do agregado (mm).
    """

    # Lista de bitolas possíveis
    listaBitolas = [5.0, 6.3, 8.0, 10.0, 12.5, 16.0, 20.0, 22.0, 25.0, 32.0, 40.0]

    # Para cada bitola, calcula o espaçamento
    for bitola in listaBitolas:
        print(f"\nOpções de espaçamento para bitola {bitola:.1f} mm:")
        alojVigas(bw, cob, diamEstribo, diamAgreg, bitola)


# Função principal para dimensionar uma viga retangular, com altura útil fornecida diretamente
def dimensionar_viga_retangular(nome, bw, h, d, M_sd, Vk, f_ck, f_yk, d_prime, E_s, stirrup_leg, combinacao, cob, diamEstribo, diamAgreg):
    """
    Função principal para dimensionar uma viga retangular, verificando flexão, cisalhamento e detalhamento dos estribos.
    A altura útil 'd' é fornecida diretamente como parâmetro.
    """

    print(f'DIMENSIONAMENTO DA VIGA RETANGULAR - {nome} ({bw} cm X {h} cm)')
    print(f'Altura útil da viga: d = {d} cm\n')

    # Cálculo da área de aço necessária para a flexão
    As, A_s_prime = calcular_area_aco_ret(M_sd, f_ck, f_yk, bw, d, d_prime, E_s, combinacao)
    print('1. DIMENSIONAMENTO À FLEXÃO:')
    print(f'Área de aço tracionada necessária (A_s): {As:.2f} cm²')
    if A_s_prime > 0:
        print(f'Área de aço comprimida necessária (A_s_prime): {A_s_prime:.2f} cm²\n')
    else:
        print('Sem necessidade de armadura comprimida (A_s_prime = 0 cm²)\n')

    # Cálculo das propriedades do concreto e aço
    props_concreto = calcular_fcd_fyd(f_ck, f_yk, combinacao)
    print('2. PROPRIEDADES DO CONCRETO E DO AÇO:')
    print(f'Propriedades calculadas: {props_concreto}\n')

    # Verificação à flexão
    λ, α_c, η_c, μ_lim, ω_lim = coefs(f_ck)
    omega = 1 - math.sqrt(1 - 2 * M_sd / (α_c * η_c * props_concreto[0] * bw * (d ** 2)))
    x = omega * d / λ
    resultado_flexao = verific_flex(bw, h, f_ck, M_sd, As, A_s_prime, omega, x, λ)
    print('3. VERIFICAÇÃO À FLEXÃO:')
    print(resultado_flexao, '\n')

    # Verificação do cisalhamento (compressão diagonal)
    compr_diag = verif_compr_diag(Vk, bw, d, f_ck, f_yk, combinacao)
    print('4. VERIFICAÇÃO DA DIAGONAL COMPRIMIDA:')
    print(compr_diag, '\n')

    # Verificação do cisalhamento (biela tracionada)
    trac_diag = verif_trac_diag(Vk, bw, d, f_ck, f_yk, combinacao)
    print('5. VERIFICAÇÃO DA BIELA TRACIONADA:')
    print(trac_diag, '\n')

    # Detalhamento dos estribos
    detalhamento = detalhar_asw(trac_diag['asw_adot (cm²/m)'], stirrup_leg)
    print('6. DETALHAMENTO DOS ESTRIBOS:')
    print(detalhamento, '\n')

    # Cálculo da quantidade de barras e bitolas disponíveis
    barras = barrasAco(As)
    print('7. QUANTIDADE DE BARRAS POR BITOLA:')
    print(barras)

    # Cálculo de espaçamentos possíveis para as bitolas
    print('8. OPÇÕES DE ESPAÇAMENTO ENTRE AS BARRAS:')
    for bitola in barras.index:
        print(f"\nOpções de espaçamento para bitola {bitola:.1f} mm:")
        ah = alojVigas(bw, cob, diamEstribo, diamAgreg, bitola)
        if len(ah) > 0:
            for i, espaçamento in enumerate(ah, start=2):
                print(f"{i} barras -> espaçamento = {espaçamento:.2f} cm")
        else:
            print("Nenhuma opção de espaçamento viável para esta bitola.")

    print('MEMORIAL DE CÁLCULO COMPLETO.')

    # Adicionando as propriedades calculadas no dicionário de resultados
    return {
        'Área de aço tracionada (cm²)': As,
        'Área de aço comprimida (cm²)': A_s_prime,
        'Propriedades calculadas': {
            'fck': f_ck,
            'fyk': f_yk
        },
        'Verificação à flexão': resultado_flexao,
        'Verificação de compressão diagonal': compr_diag,
        'Verificação de biela tracionada': trac_diag,
        'Detalhamento dos estribos': detalhamento,
        'Quantidade de barras': barras
    }

# Dicionário para valores de alfa_fg com base no tipo de seção
secoes = {
    1.2: "Seções T ou duplo T",
    1.3: "Seções I ou T invertido",
    1.5: "Seções retangulares"
}

# Dicionário para valores de eta1 com base no tipo de barra
eta1_table = {
    "lisas": 1.0,
    "dentadas": 1.4,
    "nervuradas": 2.25
}

def obter_alfa_fg(tipo_secao):
    """
    Retorna o valor de alfa_fg com base no tipo de seção.
    """
    for alfa, descricao in secoes.items():
        if descricao == tipo_secao:
            return alfa
    raise ValueError("Tipo de seção não encontrado. Verifique a descrição fornecida.")

def Mr(b, h, fck, Yt, alfa_fg):
    """
    Calcula o momento de fissuração Mr.
    """
    Ic = b * h**3 / 12
    f_ctm  = 0.3 * fck**(2/3)  # MPa
    f_ct = 0.7 * f_ctm / 10  # Conversão para kN/cm²
    Mr = alfa_fg * f_ct * Ic / Yt
    return Mr

def alfa_e(fck, Es):
    """
    Calcula o coeficiente alfa_e com base nos módulos de elasticidade do aço e concreto.
    """
    Eci = 5600 * math.sqrt(fck)
    alfa_i = 0.8 + 0.2 * (fck / 80)
    if alfa_i > 1:
      alfa_i = 1
    Ecs = alfa_i * Eci
    alfa_e = Es / Ecs
    return alfa_e

def calcular_x_II(bw, hf, bf, As, As_prime, d, alfa_e):
    """
    Calcula a posição da linha neutra em vigas "T" no Estado II Puro.
    """
    a1 = bw / 2
    a2 = hf * (bf - bw) + (alfa_e - 1) * As_prime + alfa_e * As
    a3 = -d * (alfa_e - 1) * As_prime - d * alfa_e * As - (hf**2 / 2) * (bf - bw)

    discriminante = a2**2 - 4 * a1 * a3

    if discriminante < 0:
        raise ValueError("Não há solução real para a equação.")
    else:
        xii = (-a2 + math.sqrt(discriminante)) / (2 * a1)
    return xii

def Iii(xii, bw, hf, bf, alfa_e, As, As_prime, d):
    """
    Calcula o momento de inércia no Estado II.
    """
    if xii < hf:
        Iii = (bf * (xii ** 3) / 3
                   + alfa_e * As * (xii - d) ** 2
                   + (alfa_e - 1) * As_prime * (xii - d) ** 2)
    else:
        Iii = ((bf - bw) * (hf ** 3) / 12
                   + bw * (xii ** 3) / 3
                   + (bf - bw) * ((xii - (hf / 2)) ** 2)
                   + alfa_e * As * (xii - d) ** 2
                   + (alfa_e - 1) * As_prime * (xii - d) ** 2)
    return Iii

def sigma_s(Mcf, d, xii, Iii, alfa_e):
    """
    Calcula a tensão no aço no Estádio II.
    """
    sigma_c_star = Mcf * (d - xii) / Iii
    sigma_s = sigma_c_star * alfa_e
    return sigma_s

def Wk(bw, cobrimento, phi_est, n_camadas, phi_long, a_v, As, sigma_si, Es, fctm, tipo_barra):
    """
    Calcula a abertura de fissura (Wk).
    """
    # Seleciona o valor de eta1 com base no tipo de barra
    eta1 = eta1_table.get(tipo_barra, 2.25)  # Padrão: nervurada (2.25)

    A_cri = bw * (cobrimento + phi_est + (n_camadas - 1) * phi_long + (n_camadas - 1) * a_v + phi_long / 2 + 7.5 * phi_long)
    rho_cri = As / A_cri
    wk1 = (phi_long / (12.5 * eta1)) * (sigma_si / Es) * (3 * sigma_si / fctm)
    wk2 = (phi_long / (12.5 * eta1)) * (sigma_si / Es) * (4 / rho_cri + 45)

    # Retorna o menor valor entre wk1 e wk2
    wk = min(wk1, wk2)
    return wk

def verificar_ELS_W(deformacao_calculada, L):
    """
    Verifica o Estado Limite de Serviço para deformação (ELS-W).

    Args:
        deformacao_calculada (float): Deformação obtida em mm pelo Autodesk Robot.
        L (float): Comprimento do vão da viga (em metros).

    Returns:
        dict: Resultados da verificação ELS-W.
    """
    deformacao_maxima_permitida = (L * 1000) / 250  # Deformação permitida em mm (L/250)
    if deformacao_calculada <= deformacao_maxima_permitida:
        conclusao = "Deformação dentro dos limites normativos."
    else:
        conclusao = "Deformação excede os limites normativos."

    return {
        "Deformação calculada (mm)": deformacao_calculada,
        "Deformação máxima permitida (mm)": deformacao_maxima_permitida,
        "Conclusão": conclusao
    }

# Função ELS-F que engloba todas as anteriores
def ELS_F(b, h, fck, Yt, alfa_fg, bw, hf, bf, As, As_prime, d, Es, cobrimento, phi_est, n_camadas, phi_long, a_v, Mcf, tipo_barra):
    """
    Função que engloba todas as funções para verificação do Estado Limite de Fissuração (ELS-F).
    """
    # Calcula o momento de fissuração (Mr)
    momento_fissuracao = Mr(b, h, fck, Yt, alfa_fg)

    # Calcula o coeficiente alfa_e
    alfa_e_value = alfa_e(fck, Es)

    # Calcula a posição da linha neutra no Estado II
    xii = calcular_x_II(bw, hf, bf, As, As_prime, d, alfa_e_value)

    # Calcula o momento de inércia no Estado II
    momento_inercia = Iii(xii, bw, hf, bf, alfa_e_value, As, As_prime, d)

    # Calcula a tensão no aço
    sigma_aco = sigma_s(Mcf, d, xii, momento_inercia, alfa_e_value)

    # Calcula a abertura de fissura (Wk)
    abertura_fissura = Wk(bw, cobrimento, phi_est, n_camadas, phi_long, a_v, As, sigma_aco, Es, fck, tipo_barra)

    return {
        "Momento de Fissuração": momento_fissuracao,
        "Posição da Linha Neutra (xii)": xii,
        "Momento de Inércia (Iii)": momento_inercia,
        "Tensão no Aço (σs)": sigma_aco,
        "Abertura de Fissura (Wk)": abertura_fissura
    }

# Função geral de verificação ELS que engloba ELS-F e ELS-W
def verificar_ELS(b, h, fck, Yt, tipo_secao, bw, hf, bf, As, As_prime, d, Es, cobrimento, phi_est,
                  n_camadas, phi_long, a_v, Mcf, tipo_barra, L, deformacao_calculada, limite_wk):
    """
    Verifica os Estados Limites de Serviço (ELS-F e ELS-W) para a viga.
    """
    memorial = ""
      # Obter o valor de alfa_fg com base no tipo de seção
    alfa_fg = obter_alfa_fg(tipo_secao)
    memorial += f"Tipo de Seção: {tipo_secao} (alfa_fg = {alfa_fg})\n"

    # Verificação do ELS-F (fissuração)
    momento_fissuracao = Mr(b, h, fck, Yt, alfa_fg)
    memorial += f"Momento de Fissuração (Mr): {momento_fissuracao:.2f} kN.m\n"
    memorial += f"Momento Fletor Aplicado (Mcf): {Mcf:.2f} kN.m\n"

    if Mcf < momento_fissuracao:
        memorial += "Não há fissuração (Mcf < Mr). Cálculo está adequado no ELS-F.\n\n"
    else:
        # Caso haja fissuração, prossegue com o cálculo de Wk
        memorial += "Há fissuração (Mcf ≥ Mr). Prosseguindo com o cálculo de abertura de fissura...\n"

        # Calcula o coeficiente alfa_e
        alfa_e_value = alfa_e(fck, Es)

        # Calcula a posição da linha neutra no Estado II
        xii = calcular_x_II(bw, hf, bf, As, As_prime, d, alfa_e_value)

        # Calcula o momento de inércia no Estado II
        momento_inercia = Iii(xii, bw, hf, bf, alfa_e_value, As, As_prime, d)

        # Calcula a tensão no aço
        sigma_aco = sigma_s(Mcf, d, xii, momento_inercia, alfa_e_value)

        # Calcula a abertura de fissura (Wk)
        abertura_fissura = Wk(bw, cobrimento, phi_est, n_camadas, phi_long, a_v, As, sigma_aco, Es, fck, tipo_barra)

        memorial += f"Abertura de Fissura (Wk): {abertura_fissura:.3f} mm\n"

        # Verifica se Wk está dentro do limite permitido
        if abertura_fissura <= limite_wk:
            memorial += f"Abertura de fissura dentro do limite (≤ {limite_wk:.2f} mm). Cálculo adequado no ELS-F.\n\n"
        else:
            memorial += f"Abertura de fissura fora do limite (> {limite_wk:.2f} mm). Cálculo inadequado no ELS-F.\n\n"

    # Verificação do ELS-W (deformação)
    deformacao_maxima_permitida = (L / 250) * 1000  # em mm
    memorial += f"Deformação Calculada: {deformacao_calculada:.2f} mm\n"
    memorial += f"Deformação Máxima Permitida: {deformacao_maxima_permitida:.2f} mm\n"

    if deformacao_calculada <= deformacao_maxima_permitida:
        memorial += "Deformação dentro do limite permitido. Cálculo adequado no ELS-W.\n\n"
    else:
        memorial += "Deformação fora do limite permitido. Cálculo inadequado no ELS-W.\n\n"

    # Verificação final se a viga está adequada em ambos os ELS
    if Mcf < momento_fissuracao or (abertura_fissura <= limite_wk and deformacao_calculada <= deformacao_maxima_permitida):
        memorial += "Conclusão: A viga está adequada nos parâmetros de ELS-F e ELS-W.\n"
    else:
        memorial += "Conclusão: A viga não está adequada nos parâmetros de ELS-F e/ou ELS-W.\n"

    return memorial

import streamlit as st
import math
import pandas as pd

# As funções e cálculos já fornecidos anteriormente devem ser incluídos aqui

# Função para exibir os resultados na interface do Streamlit
def exibir_memorial(memorial):
    st.subheader("Memorial de Cálculo")
    st.text(memorial)

# Função principal para criar o front-end
def main():
    st.title("Dimensionamento de Vigas de Concreto Armado")

    # Coleta de dados de entrada
    nome = st.text_input("Nome da Viga", value="Viga 1")
    bw = st.number_input("Largura da Seção (bw) [cm]", value=30)
    h = st.number_input("Altura da Seção (h) [cm]", value=50)
    d = st.number_input("Altura Útil (d) [cm]", value=45)
    M_sd = st.number_input("Momento Solicitante (M_sd) [kN.m]", value=200)
    Vk = st.number_input("Força Cortante (Vk) [kN]", value=100)
    f_ck = st.number_input("Resistência do Concreto (f_ck) [MPa]", value=30)
    f_yk = st.number_input("Tensão de Escoamento do Aço (f_yk) [MPa]", value=500)
    d_prime = st.number_input("Altura Útil da Armadura Comprimida (d') [cm]", value=5)
    E_s = st.number_input("Módulo de Elasticidade do Aço (E_s) [GPa]", value=210)
    stirrup_leg = st.number_input("Número de Pernas dos Estribos", value=2)
    combinacao = st.selectbox("Combinação de Carga", ["Permanente", "Acidental"])
    cob = st.number_input("Cobrimento Nominal (cob) [cm]", value=3)
    diamEstribo = st.number_input("Diâmetro do Estribo [mm]", value=6.3)
    diamAgreg = st.number_input("Diâmetro Máximo do Agregado [mm]", value=20)
    L = st.number_input("Comprimento do Vão (L) [m]", value=6)
    deformacao_calculada = st.number_input("Deformação Calculada [mm]", value=5.0)
    limite_wk = st.number_input("Limite de Abertura de Fissura (Wk) [mm]", value=0.3)
    tipo_secao = st.selectbox("Tipo de Seção", ["Seções T ou duplo T", "Seções I ou T invertido", "Seções retangulares"])
    tipo_barra = st.selectbox("Tipo de Barra", ["lisas", "dentadas", "nervuradas"])

    # Cálculos e verificação
    if st.button("Dimensionar Viga"):
        # Chama a função de verificação geral de ELS (incluindo ELS-F e ELS-W)
        memorial = verificar_ELS(
            bw, h, f_ck, h/2, tipo_secao, bw, h, bw, M_sd, M_sd, d, E_s, cob, diamEstribo,
            1, 12.5, 2.25, M_sd, tipo_barra, L, deformacao_calculada, limite_wk
        )

        # Exibe o memorial de cálculo
        exibir_memorial(memorial)

# Rodar o aplicativo Streamlit
if __name__ == "__main__":
    main()

