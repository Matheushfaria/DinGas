import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import griddata

def calculate_properties(beta_angles, params):
    """
    Calcula todas as propriedades do difusor para um dado conjunto de ângulos beta.
    Retorna a eficiência global e um dicionário com todas as propriedades.
    """
    b_1_deg, b_2_deg = beta_angles
    b_1_rad, b_2_rad = np.deg2rad(b_1_deg), np.deg2rad(b_2_deg)
    
    props = {}
    
    # Unpack params for easier access
    mach_inicial = params['mach_inicial']
    pressao_inicial = params['pressao_inicial']
    temperatura_inicial = params['temperatura_inicial']
    densidade_inicial = params['densidade_inicial']
    gamma = params['gamma']
    R = params['R']
    cp = params['cp']
    L_1 = params['L_1']
    L_2 = params['L_2']

    try:
        # --- Onda de Choque 1 (Oblíqua) ---
        if mach_inicial * np.sin(b_1_rad) <= 1: return -1e6, None
        
        theta_1_rad = np.arctan(2 * (1 / np.tan(b_1_rad)) * ((mach_inicial**2 * np.sin(b_1_rad)**2 - 1) / 
                              (mach_inicial**2 * (gamma + np.cos(2 * b_1_rad)) + 2)))
        props['theta_1'] = np.rad2deg(theta_1_rad)

        Mn_inf = mach_inicial * np.sin(b_1_rad)
        Mn_1_sq = (Mn_inf**2 + (2 / (gamma - 1))) / (((2 * gamma) / (gamma - 1)) * Mn_inf**2 - 1)
        if Mn_1_sq < 0: return -1e6, None
        Mn_1 = np.sqrt(Mn_1_sq)
        props['M_1'] = Mn_1 / np.sin(b_1_rad - theta_1_rad)
        
        p_r1 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn_inf**2 - 1)
        rho_r1 = ((gamma + 1) * Mn_inf**2) / (2 + (gamma - 1) * Mn_inf**2)
        T_r1 = p_r1 / rho_r1
        
        props['p_1'] = p_r1 * pressao_inicial
        props['rho_1'] = rho_r1 * densidade_inicial
        props['T_1'] = T_r1 * temperatura_inicial
        
        deltaS_1 = cp * np.log(T_r1) - R * np.log(p_r1)
        eff_1 = np.exp(-deltaS_1 / R)

        # --- Onda de Choque 2 (Oblíqua) ---
        if props['M_1'] * np.sin(b_2_rad) <= 1: return -1e6, None

        theta_2_rad = np.arctan(2 * (1 / np.tan(b_2_rad)) * ((props['M_1']**2 * np.sin(b_2_rad)**2 - 1) / 
                              (props['M_1']**2 * (gamma + np.cos(2 * b_2_rad)) + 2)))
        props['theta_2'] = np.rad2deg(theta_2_rad)

        Mn_1_2 = props['M_1'] * np.sin(b_2_rad)
        Mn_2_sq = (Mn_1_2**2 + (2 / (gamma - 1))) / (((2 * gamma) / (gamma - 1)) * Mn_1_2**2 - 1)
        if Mn_2_sq < 0: return -1e6, None
        Mn_2 = np.sqrt(Mn_2_sq)
        props['M_2'] = Mn_2 / np.sin(b_2_rad - theta_2_rad)

        p_r2 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn_1_2**2 - 1)
        rho_r2 = ((gamma + 1) * Mn_1_2**2) / (2 + (gamma - 1) * Mn_1_2**2)
        T_r2 = p_r2 / rho_r2

        props['p_2'] = p_r2 * props['p_1']
        props['rho_2'] = rho_r2 * props['rho_1']
        props['T_2'] = T_r2 * props['T_1']

        deltaS_2 = cp * np.log(T_r2) - R * np.log(p_r2)
        eff_2 = np.exp(-deltaS_2 / R)

        # --- Onda de Choque 3 (Normal) ---
        Mn_2_3 = props['M_2']
        if Mn_2_3 <= 1: return -1e6, None

        Mn_3_sq = (Mn_2_3**2 + (2 / (gamma - 1))) / (((2 * gamma) / (gamma - 1)) * Mn_2_3**2 - 1)
        if Mn_3_sq < 0: return -1e6, None
        props['M_3'] = np.sqrt(Mn_3_sq)

        p_r3 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn_2_3**2 - 1)
        rho_r3 = ((gamma + 1) * Mn_2_3**2) / (2 + (gamma - 1) * Mn_2_3**2)
        T_r3 = p_r3 / rho_r3

        props['p_3'] = p_r3 * props['p_2']
        props['rho_3'] = rho_r3 * props['rho_2']
        props['T_3'] = T_r3 * props['T_2']

        deltaS_3 = cp * np.log(T_r3) - R * np.log(p_r3)
        eff_3 = np.exp(-deltaS_3 / R)

        # --- Eficiência e Propriedades Finais ---
        global_eff = eff_1 * eff_2 * eff_3
        
        # Adiciona velocidades calculadas
        props['veltotal_inicial'] = mach_inicial * (gamma * R * temperatura_inicial)**0.5
        props['v_1'] = props['M_1'] * (gamma * R * props['T_1'])**0.5
        props['V_2'] = props['M_2'] * (gamma * R * props['T_2'])**0.5
        props['v_3'] = props['M_3'] * (gamma * R * props['T_3'])**0.5
        
        # Arrasto
        props['D'] = ((props['p_2'] - props['p_1']) * L_1 * np.sin(np.deg2rad(props['theta_1']))) + \
                     ((props['p_3'] - props['p_1']) * L_2 * np.sin(np.deg2rad(props['theta_1'] + props['theta_2'])))

        if not np.isreal(global_eff) or not np.isreal(props['D']):
            return -1e6, None

        return global_eff, props

    except (ValueError, ZeroDivisionError, RuntimeWarning):
        return -1e6, None

def objective_function(beta_angles, params):
    """
    Função objetivo para o otimizador. Retorna a eficiência negativa.
    """
    global_eff, _ = calculate_properties(beta_angles, params)
    return -global_eff # Queremos minimizar a eficiência negativa (maximizar a eficiência)

def constraints(beta_angles, params):
    """
    Define as restrições físicas para o otimizador.
    A restrição deve ser >= 0 para o método 'SLSQP'.
    """
    _, props = calculate_properties(beta_angles, params)
    
    if props is None:
        return [-1, -1, -1, -1] # Viola todas as restrições se o cálculo falhar

    # c(x) >= 0
    c1 = params['mach_inicial'] - props['M_1'] # M_1 < mach_inicial
    c2 = props['M_1'] - props['M_2']           # M_2 < M_1
    c3 = props['M_2'] - props['M_3']           # M_3 < M_2
    c4 = 1.0 - props['M_3']                    # M_3 <= 1
    
    return [c1, c2, c3, c4]

def generate_contour_data(params):
    """
    Gera dados para o gráfico de contorno através de uma busca exaustiva.
    """
    print("\n--- Gerando dados para o gráfico de contorno (pode levar um momento)... ---")
    thetal_v = []
    theta2_v = []
    eff_v = []
    
    # Usar um passo maior para acelerar a geração do gráfico
    for b_1 in range(1, 90, 2):
        for b_2 in range(1, 90, 2):
            eff, props = calculate_properties([b_1, b_2], params)
            if props is not None and eff > 0:
                # Verifica se as restrições são atendidas para um ponto válido
                if (props['M_1'] < params['mach_inicial'] and 
                    props['M_2'] < props['M_1'] and 
                    props['M_3'] < props['M_2'] and 
                    props['M_3'] <= 1):
                    thetal_v.append(props['theta_1'])
                    theta2_v.append(props['theta_2'])
                    eff_v.append(eff)
    print("--- Geração de dados concluída. ---")
    return thetal_v, theta2_v, eff_v

def main():
    """
    Script principal para otimização de difusor usando busca por gradiente com SciPy.
    """
    plt.close('all')

    # %%%%%%%%%%%%%% Parametros de entrada
    params = {
        'mach_inicial': 3.0,
        'pressao_inicial': 54048.0,
        'temperatura_inicial': 255.69,
        'densidade_inicial': 0.73643,
        'gamma': 1.4,
        'R': 287.0,
        'L_1': 49.53e-3,
        'L_2': 31.73e-3
    }
    params['cp'] = (params['gamma'] * params['R']) / (params['gamma'] - 1)

    # --- Configuração do Otimizador ---
    b0 = [30.0, 35.0]  # Chute inicial para [beta_1, beta_2]
    bounds = [(1, 89), (1, 89)] # Limites para beta_1 e beta_2
    
    # Define o dicionário de restrições para o otimizador
    cons = ({'type': 'ineq', 'fun': constraints, 'args': (params,)})

    print("--- Iniciando Otimização por Gradiente ---")
    result = minimize(
        objective_function, 
        b0, 
        args=(params,),
        method='SLSQP', # Método que suporta bounds e constraints
        bounds=bounds,
        constraints=cons,
        options={'disp': True, 'ftol': 1e-9}
    )
    print("--- Otimização Concluída ---")

    # --- Resultados ---
    if result.success:
        b_opt = result.x
        max_eff = -result.fun
        _, opt_vars = calculate_properties(b_opt, params)

        print("\n~~~~~~ Parte 1 ~~~~~~")
        print("Propriedades Difusor:")
        
        tabela_data = {
            'Estacao': ['Infinito', 'Estacao 1', 'Estacao 2', 'Estacao 3'],
            'Mach': [params['mach_inicial'], opt_vars['M_1'], opt_vars['M_2'], opt_vars['M_3']],
            'Pressao [Pa]': [params['pressao_inicial'], opt_vars['p_1'], opt_vars['p_2'], opt_vars['p_3']],
            'Temperatura [K]': [params['temperatura_inicial'], opt_vars['T_1'], opt_vars['T_2'], opt_vars['T_3']],
            'Velocidade [m/s]': [opt_vars['veltotal_inicial'], opt_vars['v_1'], opt_vars['V_2'], opt_vars['v_3']],
            'Massa especifica [Kg/m3]': [params['densidade_inicial'], opt_vars['rho_1'], opt_vars['rho_2'], opt_vars['rho_3']]
        }
        T = pd.DataFrame(tabela_data)
        print(T.to_string())

        print(f"\nTheta (graus): Theta_1:{opt_vars['theta_1']:.3f} Theta_2:{opt_vars['theta_2']:.3f}")
        print(f"Beta (graus): Beta_1:{b_opt[0]:.3f}; Beta_2:{b_opt[1]:.3f}")
        print(f"Eficiencia [%%]: {max_eff*100:.3f}")
        print(f"Arrasto [N]: {opt_vars['D']:.3f}")

        # --- Geração de Gráficos ---
        # Gráfico 1: Contorno de Eficiência
        thetal_v, theta2_v, eff_v = generate_contour_data(params)
        
        fig1, ax1 = plt.subplots()
        if thetal_v:
            xi = np.linspace(min(thetal_v), max(thetal_v), 100)
            yi = np.linspace(min(theta2_v), max(theta2_v), 100)
            Xq, Yq = np.meshgrid(xi, yi)
            Zq = griddata((thetal_v, theta2_v), eff_v, (Xq, Yq), method='cubic')
            
            contour = ax1.contour(Xq, Yq, Zq, 20, cmap='jet')
            fig1.colorbar(contour, ax=ax1, label=r'$\eta$')
            
            max_value_formatted = f'{max(eff_v):.3f}'
            ax1.plot(opt_vars['theta_1'], opt_vars['theta_2'], 'r*', markersize=10, label=f'Ótimo ({max_eff:.3f})')
            ax1.text(opt_vars['theta_1'], opt_vars['theta_2'], f'  $\\eta_{{max}} = {max_eff:.3f}$', color='red')

        ax1.set_xlabel(r'$\theta_1$ (graus)')
        ax1.set_ylabel(r'$\theta_2$ (graus)')
        ax1.set_title('Superfície de Eficiência do Difusor')
        ax1.grid(True)
        ax1.legend()

        # --- Gráfico 2: Propriedades Termodinâmicas (Temp, Mach, Pressão) ---
        fig2, host2 = plt.subplots(figsize=(9, 7))
        fig2.subplots_adjust(left=0.4) # Aumenta o espaço à esquerda
        x = np.arange(1, 5)
        names = [r'$X_{\infty}$', r'$X_1$', r'$X_2$', r'$X_3$']
        
        axes2 = {"Mach": host2.twinx(), "Pressao": host2.twinx()}
        
        # Move os eixos para a esquerda
        axes2["Mach"].spines['right'].set_visible(False)
        axes2["Mach"].spines['left'].set_visible(True)
        axes2["Mach"].yaxis.tick_left()
        axes2["Mach"].yaxis.set_label_position('left')

        axes2["Pressao"].spines['right'].set_visible(False)
        axes2["Pressao"].spines['left'].set_visible(True)
        axes2["Pressao"].spines['left'].set_position(('outward', 80))
        axes2["Pressao"].yaxis.tick_left()
        axes2["Pressao"].yaxis.set_label_position('left')

        # Plotando os dados
        p1, = host2.plot(x, T['Temperatura [K]'], '-^', color=[0, 0.392, 0], label='Temperatura', drawstyle='steps-post')
        p2, = axes2["Mach"].plot(x, T['Mach'], '-o', color='black', label='Mach', drawstyle='steps-post')
        p3, = axes2["Pressao"].plot(x, T['Pressao [Pa]'], '-s', color=[0.494, 0.184, 0.556], label='Pressao', drawstyle='steps-post')

        # Configurando eixos e labels
        host2.set_xticks(x); host2.set_xticklabels(names); host2.set_xlim(0.5, 4.5)
        host2.set_ylabel('Temperatura [K]'); host2.yaxis.label.set_color(p1.get_color())
        axes2["Mach"].set_ylabel('Mach'); axes2["Mach"].yaxis.label.set_color(p2.get_color())
        axes2["Pressao"].set_ylabel('Pressao [Pa]'); axes2["Pressao"].yaxis.label.set_color(p3.get_color())
        
        host2.legend(handles=[p1, p2, p3], loc='center left')
        host2.grid(True)
        fig2.suptitle('Propriedades Termodinâmicas nas Estações')

        # --- Gráfico 3: Propriedades de Escoamento (Massa Específica, Velocidade) ---
        fig3, host3 = plt.subplots(figsize=(9, 7))
        fig3.subplots_adjust(left=0.3) # Aumenta o espaço à esquerda
        
        axes3 = {"Velocidade": host3.twinx()}
        
        # Move o eixo para a esquerda
        axes3["Velocidade"].spines['right'].set_visible(False)
        axes3["Velocidade"].spines['left'].set_visible(True)
        axes3["Velocidade"].spines['left'].set_position(('outward', 80))
        axes3["Velocidade"].yaxis.tick_left()
        axes3["Velocidade"].yaxis.set_label_position('left')

        # Plotando os dados
        p4, = host3.plot(x, T['Massa especifica [Kg/m3]'], '-d', color=[0, 0.447, 0.741], label='Massa especifica', drawstyle='steps-post')
        p5, = axes3["Velocidade"].plot(x, T['Velocidade [m/s]'], '-*', color=[0.85, 0.325, 0.098], label='Velocidade', drawstyle='steps-post')

        # Configurando eixos e labels
        host3.set_xticks(x); host3.set_xticklabels(names); host3.set_xlim(0.5, 4.5)
        host3.set_ylabel('Massa especifica [kg/m3]'); host3.yaxis.label.set_color(p4.get_color())
        axes3["Velocidade"].set_ylabel('Velocidade [m/s]'); axes3["Velocidade"].yaxis.label.set_color(p5.get_color())
        
        host3.legend(handles=[p4, p5], loc='center left')
        host3.grid(True)
        fig3.suptitle('Propriedades de Escoamento nas Estações')
        
        plt.show()
    else:
        print("\nOtimização falhou.")
        print(result.message)

if __name__ == '__main__':
    main()