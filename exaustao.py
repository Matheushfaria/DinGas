import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def main():
    """
    Tradução do script MATLAB para Python para análise de difusor com ondas de choque.
    """
    # Limpa figuras anteriores, se houver
    plt.close('all')

    # %%%%%%%%%%%%%% Parametros de entrada
    mach_inicial = 3.0
    pressao_inicial = 54048.0  # [Pa]
    temperatura_inicial = 255.69  # [K]
    densidade_inicial = 0.73643  # [kg/m3]
    gamma = 1.4
    R = 287.0  # [J/kgK]
    L_1 = 49.53e-3  # Comprimento estação 1 [m]
    L_2 = 31.73e-3  # Comprimento estação 2 [m]

    # Propriedades calculadas
    velsom_inicial = (gamma * R * temperatura_inicial)**0.5
    veltotal_inicial = mach_inicial * velsom_inicial
    cp = (gamma * R) / (gamma - 1)

    # %%%%%%%%%%%%%%%%%% Início da Otimização por Busca Exaustiva
    max_eff = 0.0
    
    # Variáveis para armazenar os valores ótimos
    opt_vars = {}

    # Listas para o gráfico de contorno
    thetal_v = []
    theta2_v = []
    eff_v = []

    # Loop principal
    for b_1 in range(1, 90):
        # --- Primeira Onda de Choque (Oblíqua) ---
        try:
            # Converte para radianos para os cálculos com numpy
            b_1_rad = np.deg2rad(b_1)
            
            # Condição para existência de choque
            if mach_inicial * np.sin(b_1_rad) <= 1:
                continue

            theta_1_rad = np.arctan(2 * (1 / np.tan(b_1_rad)) * ((mach_inicial*2 * np.sin(b_1_rad)*2 - 1) / 
                                  (mach_inicial**2 * (gamma + np.cos(2 * b_1_rad)) + 2)))
            theta_1 = np.rad2deg(theta_1_rad)

            Mn_inf = mach_inicial * np.sin(b_1_rad)
            Mn_1_sq = (Mn_inf*2 + (2 / (gamma - 1))) / (((2 * gamma) / (gamma - 1)) * Mn_inf*2 - 1)
            if Mn_1_sq < 0: continue
            Mn_1 = np.sqrt(Mn_1_sq)

            M_1 = Mn_1 / np.sin(b_1_rad - theta_1_rad)
            
            p_r1 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn_inf**2 - 1)
            rho_r1 = ((gamma + 1) * Mn_inf*2) / (2 + (gamma - 1) * Mn_inf*2)
            T_r1 = p_r1 / rho_r1

            p_1 = p_r1 * pressao_inicial
            rho_1 = rho_r1 * densidade_inicial
            T_1 = T_r1 * temperatura_inicial
            
            deltaS_1 = cp * np.log(T_r1) - R * np.log(p_r1)
            eff_1 = np.exp(-deltaS_1 / R)
            v_1 = M_1 * (gamma * R * T_1)**0.5

        except (ValueError, ZeroDivisionError):
            continue

        for b_2 in range(1, 90):
            # --- Segunda Onda de Choque (Oblíqua) ---
            try:
                if not np.isreal(M_1) or M_1 > mach_inicial:
                    continue
                
                b_2_rad = np.deg2rad(b_2)
                
                if M_1 * np.sin(b_2_rad) <= 1:
                    continue

                theta_2_rad = np.arctan(2 * (1 / np.tan(b_2_rad)) * ((M_1*2 * np.sin(b_2_rad)*2 - 1) / 
                                      (M_1**2 * (gamma + np.cos(2 * b_2_rad)) + 2)))
                theta_2 = np.rad2deg(theta_2_rad)

                Mn_1_2 = M_1 * np.sin(b_2_rad)
                Mn_2_sq = (Mn_1_2*2 + (2 / (gamma - 1))) / (((2 * gamma) / (gamma - 1)) * Mn_1_2*2 - 1)
                if Mn_2_sq < 0: continue
                Mn_2 = np.sqrt(Mn_2_sq)

                M_2 = Mn_2 / np.sin(b_2_rad - theta_2_rad)

                p_r2 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn_1_2**2 - 1)
                rho_r2 = ((gamma + 1) * Mn_1_2*2) / (2 + (gamma - 1) * Mn_1_2*2)
                T_r2 = p_r2 / rho_r2

                p_2 = p_r2 * p_1
                rho_2 = rho_r2 * rho_1
                T_2 = T_r2 * T_1

                deltaS_2 = cp * np.log(T_r2) - R * np.log(p_r2)
                eff_2 = np.exp(-deltaS_2 / R)
                V_2 = M_2 * (gamma * R * T_2)**0.5

                # --- Terceira Onda de Choque (Normal) ---
                if not np.isreal(M_2) or M_2 > M_1:
                    continue
                
                Mn_2_3 = M_2
                if Mn_2_3 <= 1: continue

                Mn_3_sq = (Mn_2_3*2 + (2 / (gamma - 1))) / (((2 * gamma) / (gamma - 1)) * Mn_2_3*2 - 1)
                if Mn_3_sq < 0: continue
                M_3 = np.sqrt(Mn_3_sq)

                p_r3 = 1 + ((2 * gamma) / (gamma + 1)) * (Mn_2_3**2 - 1)
                rho_r3 = ((gamma + 1) * Mn_2_3*2) / (2 + (gamma - 1) * Mn_2_3*2)
                T_r3 = p_r3 / rho_r3

                p_3 = p_r3 * p_2
                rho_3 = rho_r3 * rho_2
                T_3 = T_r3 * T_2

                deltaS_3 = cp * np.log(T_r3) - R * np.log(p_r3)
                eff_3 = np.exp(-deltaS_3 / R)
                v_3 = M_3 * (gamma * R * T_3)**0.5

                # --- Eficiência e Verificações Finais ---
                global_eff = eff_1 * eff_2 * eff_3

                if np.isreal(global_eff) and M_3 < M_2 < M_1:
                    thetal_v.append(theta_1)
                    theta2_v.append(theta_2)
                    eff_v.append(global_eff)

                    if global_eff > max_eff and M_3 <= 1:
                        max_eff = global_eff
                        opt_vars['theta_1_opt'] = theta_1
                        opt_vars['theta_2_opt'] = theta_2
                        opt_vars['betal_opt'] = b_1
                        opt_vars['beta2_opt'] = b_2
                        opt_vars['numero_mach'] = [mach_inicial, M_1, M_2, M_3]
                        opt_vars['rho'] = [densidade_inicial, rho_1, rho_2, rho_3]
                        opt_vars['pressure'] = [pressao_inicial, p_1, p_2, p_3]
                        opt_vars['temperature'] = [temperatura_inicial, T_1, T_2, T_3]
                        opt_vars['fluid_velocity'] = [veltotal_inicial, v_1, V_2, v_3]
                        opt_vars['D'] = ((p_2 - p_1) * L_1 * np.sin(np.deg2rad(theta_1))) + \
                                        ((p_3 - p_1) * L_2 * np.sin(np.deg2rad(theta_1 + theta_2)))

            except (ValueError, ZeroDivisionError, RuntimeWarning):
                continue

    # %%%%%%%%%%%%%%%% Resultados
    print("~~~~~~ Parte 1 ~~~~~~")
    print("Propriedades Difusor:")
    
    # Criar DataFrame com Pandas
    tabela_data = {
        'Estacao': ['Infinito', 'Estacao 1', 'Estacao 2', 'Estacao 3'],
        'Mach': opt_vars['numero_mach'],
        'Pressao [Pa]': opt_vars['pressure'],
        'Temperatura [K]': opt_vars['temperature'],
        'Velocidade [m/s]': opt_vars['fluid_velocity'],
        'Massa especifica [Kg/m3]': opt_vars['rho']
    }
    T = pd.DataFrame(tabela_data)
    print(T.to_string())

    print(f"\nTheta (graus): Theta_1:{opt_vars['theta_1_opt']:.3f} Theta_2:{opt_vars['theta_2_opt']:.3f}")
    print(f"Beta (graus): Beta_1:{opt_vars['betal_opt']:.3f}; Beta_2:{opt_vars['beta2_opt']:.3f}")
    print(f"Eficiencia [%%]: {max_eff*100:.3f}")
    print(f"Arrasto [N]: {opt_vars['D']:.3f}")

    # %%%%%%%%%%%%%%%% Gráficos
    # --- Gráfico de Contorno ---
    fig1, ax1 = plt.subplots()
    if thetal_v: # Apenas plota se houver dados
        xi = np.linspace(min(thetal_v), max(thetal_v), 100)
        yi = np.linspace(min(theta2_v), max(theta2_v), 100)
        Xq, Yq = np.meshgrid(xi, yi)
        Zq = griddata((thetal_v, theta2_v), eff_v, (Xq, Yq), method='cubic')
        
        contour = ax1.contour(Xq, Yq, Zq, 20, cmap='jet')
        cbar = fig1.colorbar(contour, ax=ax1)
        cbar.set_label(r'$\eta$', fontsize=14)
        
        max_value_formatted = f'{max(eff_v):.3f}'
        ax1.text(opt_vars['theta_1_opt'], opt_vars['theta_2_opt'], rf'$\eta_{{max}} = {max_value_formatted}$', 
                 color='black', fontsize=12, ha='center')
        
    ax1.set_xlabel(r'$\theta_1$ (graus)')
    ax1.set_ylabel(r'$\theta_2$ (graus)')
    ax1.grid(True)
    ax1.set_title('Gráfico de Eficiência')

    # --- Gráfico de Múltiplos Eixos ---
    fig2, host = plt.subplots(figsize=(9, 7))
    fig2.subplots_adjust(right=0.6) # Ajusta o espaço para os eixos Y

    x = np.arange(1, 5)
    names = ['Estacao Inf', 'Estacao 1', 'Estacao 2', 'Estacao 3']
    
    # Dicionário para os eixos
    axes = {
        "host": host,
        "Velocidade": host.twinx(),
        "Pressao": host.twinx(),
        "Massa": host.twinx(),
        "Temperatura": host.twinx()
    }
    
    # Posiciona os eixos Y para não sobrepor
    axes["Pressao"].spines['right'].set_position(('outward', 60))
    axes["Massa"].spines['right'].set_position(('outward', 120))
    axes["Temperatura"].spines['right'].set_position(('outward', 180))

    # Plotando os dados
    p1, = axes["host"].plot(x, opt_vars['numero_mach'], '-o', color='black', label='Mach')
    p2, = axes["Velocidade"].plot(x, opt_vars['fluid_velocity'], '-*', color=[0.8500, 0.3250, 0.0980], label='Velocidade')
    p3, = axes["Pressao"].plot(x, opt_vars['pressure'], '-s', color=[0.4940, 0.1840, 0.5560], label='Pressao')
    p4, = axes["Massa"].plot(x, opt_vars['rho'], '-d', color=[0, 0.4470, 0.7410], label='Massa especifica')
    p5, = axes["Temperatura"].plot(x, opt_vars['temperature'], '-^', color=[0, 0.3922, 0], label='Temperatura')

    # Configurando eixos e labels
    axes["host"].set_xticks(x)
    axes["host"].set_xticklabels(names)
    axes["host"].set_xlim(0.5, 4.5)
    axes["host"].set_ylabel('Mach')
    axes["host"].yaxis.label.set_color(p1.get_color())

    axes["Velocidade"].set_ylabel('Velocidade [m/s]')
    axes["Velocidade"].yaxis.label.set_color(p2.get_color())
    
    axes["Pressao"].set_ylabel('Pressao [Pa]')
    axes["Pressao"].yaxis.label.set_color(p3.get_color())

    axes["Massa"].set_ylabel('Massa especifica [kg/m3]')
    axes["Massa"].yaxis.label.set_color(p4.get_color())

    axes["Temperatura"].set_ylabel('Temperatura [K]')
    axes["Temperatura"].yaxis.label.set_color(p5.get_color())

    # Legenda
    host.legend(handles=[p1, p2, p3, p4, p5], loc='best')
    
    host.grid(True)
    host.set_aspect('auto') # Equivalente a 'axis square' é mais complexo, 'auto' é o padrão
    
    plt.show()

if __name__ == "__main__":
    main()
