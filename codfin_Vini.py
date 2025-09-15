import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

#===============================================
# Cálculos do difusor
#===============================================

# Constantes do ar
gamma = 1.4
c_p = 1005.0
R = 287.0

# Dados iniciais
props_inf = {
    'mach': 3.0,
    'pressure': 54048.0 / 1e6, # Convertido para MPa
    'temperature': 255.69,
    'rho': 0.73643,
}
# -------------------------------------------------
# Função que calcula a onda de choque (beta em graus)
# ----------------------------------------------

def shock_wave(beta, props):
    beta_rad = np.deg2rad(beta)
    m1 = props['mach']
    p1 = props['pressure']
    T1 = props['temperature']
    rho1 = props['rho']

    if np.tan(beta_rad) == 0 or m1 * np.sin(beta_rad) <= 1:
        return 0.0, props, 0.0

    theta = np.arctan(
        (2 / np.tan(beta_rad)) *
        ((m1**2) * (np.sin(beta_rad)**2) - 1) /
        (m1**2 * (gamma + np.cos(2 * beta_rad)) + 2)
    )
    m_n1 = m1 * np.sin(beta_rad)
    m_n2 = np.sqrt((1 + ((gamma - 1) / 2) * m_n1**2) / (gamma * m_n1**2 - (gamma - 1) / 2))
    delta_s = c_p * np.log(
        (1 + (2 * gamma * (m_n1**2 - 1)) / (gamma + 1)) *
        (2 + (gamma - 1) * m_n1**2) /
        ((gamma + 1) * m_n1**2)
    ) - R * np.log(1 + (2 * gamma * (m_n1**2 - 1)) / (gamma + 1))
    stc_pressure_ratio = np.exp(-delta_s / R)
    m2 = m_n2 / np.sin(beta_rad - theta)

    rho2 = (((gamma + 1) * m_n1**2) / (2 + (gamma - 1) * m_n1**2)) * rho1
    p2 = (1 + ((2 * gamma) / (gamma + 1)) * (m_n1**2 - 1)) * p1
    T2 = (p2 / p1) * (rho1 / rho2) * T1

    props_out = {
        'mach': m2,
        'pressure': p2,
        'temperature': T2,
        'rho': rho2
    }

    return stc_pressure_ratio, props_out, np.rad2deg(theta)

# -------------------------
# Parte da otimização (busca por gradiente)
# -------------------------

# Históricos
eff_history = []
theta_history = []
path_history = []

# Função objetivo
def objective(betas):
    beta1, beta2 = betas
    eff1, props1, theta1 = shock_wave(beta1, props_inf)
    eff2, props2, theta2 = shock_wave(beta2, props1)
    eff3, _, _ = shock_wave(90, props2)
    return -(eff1 * eff2 * eff3)

# Callback
def callback(xk):
    beta1, beta2 = xk
    eff1, props1, theta1 = shock_wave(beta1, props_inf)
    eff2, props2, theta2 = shock_wave(beta2, props1)
    eff3, _, _ = shock_wave(90, props2)
    eff_total = eff1 * eff2 * eff3

    eff_history.append(eff_total)
    theta_history.append((theta1, theta2))
    path_history.append((theta1, theta2))

# Otimização
beta_guess = [30.0, 30.0]
bounds = [(1, 89), (1, 89)]
res = minimize(objective, beta_guess, method='L-BFGS-B', bounds=bounds, callback=callback)

opt_beta1, opt_beta2 = res.x
max_eff_total = -res.fun

# Resultados finais
eff1, props1, theta1 = shock_wave(opt_beta1, props_inf)
eff2, props2, theta2 = shock_wave(opt_beta2, props1)
eff3, props3, _ = shock_wave(90, props2)

#=====================================
# Visualização dos resultados
#====================================

# Função para velocidade total
def velocidade_total(M, T):
    return M * np.sqrt(gamma * R * T)

# Lista com os nomes das posições
posicoes = ['X_inf', 'X_1', 'X_2', 'X_3']

# Propriedades organizadas
mach_vals = [props_inf['mach'], props1['mach'], props2['mach'], props3['mach']]
pressao_vals = [props_inf['pressure'], props1['pressure'], props2['pressure'], props3['pressure']]
temperatura_vals = [props_inf['temperature'], props1['temperature'], props2['temperature'], props3['temperature']]
densidade_vals = [props_inf['rho'], props1['rho'], props2['rho'], props3['rho']]
velocidade_vals = [velocidade_total(M, T) for M, T in zip(mach_vals, temperatura_vals)]

# Tabela com os dados
tabela_data = {
    'Posicao': posicoes,
    'Mach': mach_vals,
    'Pressao [MPa]': pressao_vals,
    'Temperatura [K]': temperatura_vals,
    'Velocidade [m/s]': velocidade_vals,
    'Massa especifica [kg/m³]': densidade_vals
}

T = pd.DataFrame(tabela_data)
print(T.to_string(index=False, float_format='{:.3f}'.format))

print(f"Angulo otimo beta1: {opt_beta1:.2f} graus -> theta1 = {theta1:.2f} graus")
print(f"Angulo otimo beta2: {opt_beta2:.2f} graus -> theta2 = {theta2:.2f} graus")
print(f"Eficiencia total maxima: {100*max_eff_total:.3f}%")

# ---------------------------
# 1. Convergencia da eficiencia
# ---------------------------
plt.figure(figsize=(8, 5))
plt.plot(eff_history, marker='o', color='green')
plt.title("Convergencia da Eficiencia Total")
plt.xlabel("Iteracao")
plt.ylabel("Eficiencia Total")
plt.grid(True)
plt.tight_layout()
plt.show()

# ---------------------------
# 2. Curvas de nivel com trajetoria (theta1 x theta2) - estilo preto e branco
# ---------------------------

# Gera malha de beta1 x beta2
beta1_range = np.linspace(1, 89, 100)
beta2_range = np.linspace(1, 89, 100)
B1, B2 = np.meshgrid(beta1_range, beta2_range)
Z_eff = np.zeros_like(B1)
Theta1 = np.zeros_like(B1)
Theta2 = np.zeros_like(B2)

for i in range(B1.shape[0]):
    for j in range(B1.shape[1]):
        b1 = B1[i, j]
        b2 = B2[i, j]
        eff1, props1, th1 = shock_wave(b1, props_inf)
        eff2, props2, th2 = shock_wave(b2, props1)
        eff3, _, _ = shock_wave(90, props2)

        eff_total = eff1 * eff2 * eff3
        Z_eff[i, j] = eff_total
        Theta1[i, j] = th1
        Theta2[i, j] = th2

# Prepara interpolação para o espaço theta1 x theta2
theta1_vals = Theta1.flatten()
theta2_vals = Theta2.flatten()
eff_vals = Z_eff.flatten()

mask = ~np.isnan(theta1_vals) & ~np.isnan(theta2_vals) & (eff_vals > 0)
theta1_vals = theta1_vals[mask]
theta2_vals = theta2_vals[mask]
eff_vals = eff_vals[mask]

th1_grid = np.linspace(min(theta1_vals), max(theta1_vals), 200)
th2_grid = np.linspace(min(theta2_vals), max(theta2_vals), 200)
TH1, TH2 = np.meshgrid(th1_grid, th2_grid)

Z_interp = griddata(
    (theta1_vals, theta2_vals),
    eff_vals,
    (TH1, TH2),
    method='cubic'
)

# Geração do gráfico com linhas pretas e níveis rotulados
plt.figure(figsize=(8, 6))
contour = plt.contour(TH1, TH2, Z_interp, levels=15, colors='black', linewidths=1)
plt.clabel(contour, inline=True, fontsize=9, fmt="%.3f")

# Trajetória do otimizador
path = np.array(path_history)
plt.plot(path[:, 0], path[:, 1], color='gray', linestyle='--', marker='o', label='Trajetoria')

# Ponto ótimo
theta1_opt, theta2_opt = path[-1]
plt.scatter(theta1_opt, theta2_opt, s=100, color='red', edgecolor='black', zorder=5, label='Otimo')

plt.title("Curvas de Nivel - Eficiencia Total (theta1 x theta2)")
plt.xlabel("theta1 (graus)")
plt.ylabel("theta2 (graus)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#=============================================
# 3. Gráficos de propriedades ao longo do difusor
#============================================

# ================== Preparação dos dados ==================
x = [1, 2, 3, 4]
names = [r'$X_{\infty}$', r'$X_1$', r'$X_2$', r'$X_3$']

# converte pressão para MPa (já está em MPa, então não é necessário)
# pressao_vals_mpa = [p for p in pressao_vals]

# função para estender após X3 (passo final "flat")
def extend(x_list, y_list, dx=1):
    return x_list + [x_list[-1] + dx], y_list + [y_list[-1]]

# séries estendidas (para "steps-post")
xT, yT = extend(x, temperatura_vals)
xM, yM = extend(x, mach_vals)
xP, yP = extend(x, pressao_vals)
xr, yr = extend(x, densidade_vals)
xV, yV = extend(x, velocidade_vals)

# ================== Figura 1: T, M, P ==================
fig1, axs1 = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

# Temperatura
color_T = 'red'
idx_T = 0
axs1[idx_T].plot(xT, yT, '-', color=color_T, drawstyle='steps-post', lw=2.5)
axs1[idx_T].set_ylabel('Temperatura [K]')
axs1[idx_T].grid(True, linestyle='--', alpha=0.3)
axs1[idx_T].spines['bottom'].set_visible(False)
axs1[idx_T].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
axs1[idx_T].tick_params(axis='y', colors=color_T)
axs1[idx_T].yaxis.label.set_color(color_T)
axs1[idx_T].spines['left'].set_color(color_T)
axs1[idx_T].yaxis.set_major_locator(MaxNLocator(nbins=10))  

# Mach
color_M = 'black'
idx_M = 1
axs1[idx_M].plot(xM, yM, '-', color=color_M, drawstyle='steps-post', lw=2.5)
axs1[idx_M].set_ylabel('Mach')
axs1[idx_M].grid(True, linestyle='--', alpha=0.3)
axs1[idx_M].spines['bottom'].set_visible(False)
axs1[idx_M].spines['top'].set_visible(False)
axs1[idx_M].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
axs1[idx_M].tick_params(axis='y', colors=color_M)
axs1[idx_M].yaxis.label.set_color(color_M)
axs1[idx_M].spines['left'].set_color(color_M)
axs1[idx_M].yaxis.set_major_locator(MaxNLocator(nbins=10))

# Pressão
color_P = 'blue'
idx_P = 2
axs1[idx_P].plot(xP, yP, '-', color=color_P, drawstyle='steps-post', lw=2.5)
axs1[idx_P].set_ylabel('Pressão [MPa]')
axs1[idx_P].grid(True, linestyle='--', alpha=0.3)
axs1[idx_P].spines['top'].set_visible(False)
axs1[idx_P].tick_params(axis='y', colors=color_P)
axs1[idx_P].yaxis.label.set_color(color_P)
axs1[idx_P].spines['left'].set_color(color_P)
axs1[idx_P].yaxis.set_major_locator(MaxNLocator(nbins=10))

# Linhas verticais e ajustes eixo x
for ax in axs1:
    for xc in [2, 3, 4]:
        ax.axvline(xc, color='lightgray', linestyle=':', lw=1)

axs1[-1].set_xticks(x)
axs1[-1].set_xticklabels(names)
axs1[-1].set_xlim(0.5, 5.2)
axs1[-1].set_xlabel('Estações')

fig1.suptitle('Temperatura, Mach e Pressão ao longo do difusor')
fig1.tight_layout(rect=[0, 0, 1, 0.95])
fig1.subplots_adjust(hspace=0)

# ================== Figura 2: ρ, V ==================
fig2, axs2 = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# Massa específica
color_rho = 'green'
idx_rho = 0
axs2[idx_rho].plot(xr, yr, '-', color=color_rho, drawstyle='steps-post', lw=2.5)
axs2[idx_rho].set_ylabel('Massa específica [kg/m³]')
axs2[idx_rho].grid(True, linestyle='--', alpha=0.3)
axs2[idx_rho].spines['bottom'].set_visible(False)
axs2[idx_rho].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
axs2[idx_rho].tick_params(axis='y', colors=color_rho)
axs2[idx_rho].yaxis.label.set_color(color_rho)
axs2[idx_rho].spines['left'].set_color(color_rho)
axs2[idx_rho].yaxis.set_major_locator(MaxNLocator(nbins=8))

# Velocidade
color_V = 'black'
idx_V = 1
axs2[idx_V].plot(xV, yV, '-', color=color_V, drawstyle='steps-post', lw=2.5)
axs2[idx_V].set_ylabel('Velocidade [m/s]')
axs2[idx_V].grid(True, linestyle='--', alpha=0.3)
axs2[idx_V].spines['top'].set_visible(False)
axs2[idx_V].tick_params(axis='y', colors=color_V)
axs2[idx_V].yaxis.label.set_color(color_V)
axs2[idx_V].spines['left'].set_color(color_V)
axs2[idx_V].yaxis.set_major_locator(MaxNLocator(nbins=10))

# Linhas verticais e ajustes eixo x
for ax in axs2:
    for xc in [2, 3, 4]:
        ax.axvline(xc, color='lightgray', linestyle=':', lw=1)

axs2[-1].set_xticks(x)
axs2[-1].set_xticklabels(names)
axs2[-1].set_xlim(0.5, 5.2)
axs2[-1].set_xlabel('Estações')

fig2.suptitle('Massa específica e Velocidade ao longo do difusor')
fig2.tight_layout(rect=[0, 0, 1, 0.95])
fig2.subplots_adjust(hspace=0)

plt.show()