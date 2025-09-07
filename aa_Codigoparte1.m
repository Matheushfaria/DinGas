%% Dinamica dos Gases 1/2025 - Work 1 Parte 1

clear;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parte 1 - Difusor
%%%%%%%%%%%%%% Parametros de entrada

mach_inicial = 3; % Mach
pressao_inicial = 54048; % [Pa]
temperatura_inicial = 255.69; % [K]
densidade_inicial = 0.73643; % [kg/m3]
gamma = 1.4;
R = 287; % [J/kgK]
L_1 = 49.53e-3; %Comprimento estação 1 [m]
L_2 = 31.73e-3; %Comprimento estação 2 [m]

% Velocidade do som [m/s]
velsom_inicial = (gamma * R * temperatura_inicial)^(1/2); 

% Pressao total
pressaotot_inicial = pressao_inicial * (1+((gamma - 1) /2) * mach_inicial^2)^(gamma/(gamma - 1));

% Velocidade [m/s]
veltotal_inicial = mach_inicial * velsom_inicial;
cp = (gamma * R)/(gamma - 1);

%%%%%%%%%%%%%%%%%% Valores apos a primeira onda de choque (obliqua)
max_eff = 0;
i = 0;
for b_1 = 1:1:89
    theta_1 = atand (2 * cotd (b_1)*(((mach_inicial^2) * (sind(b_1)^2) - 1) / ...
    ((mach_inicial^2) * (gamma + cosd (2* b_1)) + 2)));

    % Mn = Mach normal
    Mn_inf = mach_inicial * sind (b_1);
    Mn_1 = sqrt ((Mn_inf^2 + (2/(gamma-1)))/((((2* gamma)/(gamma - 1)) *Mn_inf ^2)-1));

    % Razao de pressoes
    p_r1 = 1 + ((2 * gamma) / (gamma + 1)) * ((Mn_inf^2) - 1);

    % Pressao
    p_1 = p_r1 * pressao_inicial;

    % Razão de massa especifica
    rho_r1 = ((gamma + 1) * (Mn_inf^2)) / (2 + (gamma - 1) * (Mn_inf^2));

    % Massa especifica
    rho_1 = rho_r1 * densidade_inicial;

    % Razao de temperaturas
    T_r1 = p_r1 / rho_r1;

    % Temperatura
    T_1 = T_r1 * temperatura_inicial;

    % Numero de Mach
    M_1 = Mn_1 / sind (b_1 - theta_1);
    deltaS_1 = cp*log(T_r1) - (R*log(p_r1));

    % Eficiencia
    eff_1 = exp((-(deltaS_1)/R));

    % Velocidade do som
    a_1 = (gamma*R*T_1)^(1/2);

    % Velocidade do escoamento
    v_1 = M_1*a_1;

    %%%%%%%%%%%%%%%%%% Valores apos segunda onda de choque (obliqua)
    for b_2 = 1:1:89
        % Verificar se M1 real e menor que mach_inicial
        if isreal (M_1) == false && M_1 > mach_inicial
            continue
        end
        theta_2 = atand ((2 * cotd (b_2) * (((M_1^2) *sind (b_2)^2) - 1))/((M_1 ^2) * (gamma + cosd (2*b_2)) + 2));
        Mn_1 = M_1 * sind (b_2);
        
        % Verificar se Mn_I real
        if isreal (Mn_1) == false
            continue
        end
        Mn_2 = sqrt((Mn_1 ^2+(2 /(gamma-1)))/((((2* gamma) / (gamma - 1)) *Mn_1 ^2) -1));

        % Verificar se Mn_2 real
        if isreal (Mn_2) == false
            continue
        end

        % Numero de Mach
        M_2 = Mn_2 / sind (b_2 - theta_2);

        % Razao de massa especifica
        rho_r2 = ((gamma + 1) * Mn_1^2)/(2 + ((gamma - 1) * Mn_1^2));

        % Massa especifica
        rho_2 = rho_1*rho_r2;

        % Razao de pressoes
        p_r2 = 1+(((2*gamma)/(gamma+1))*((Mn_1^2)-1));

        % Pressao
        p_2 = p_1*p_r2;

        % Razao de temperaturas
        T_r2 = p_r2/rho_r2;

        % Temperatura
        T_2 = T_1*T_r2;
        deltas_2 = cp*log(T_r2) - (R*log(p_r2));

        % Eficiencia
        eff_2 = exp((-(deltas_2)) /R);

        % Velocidade do som
        a_2 = (gamma*R*T_2)^(1/2);

        % Velocidade do escoamento
        V_2 = M_2*a_2;

        % Verificar se M_2 real e menor que M_1
        if isreal (M_2)==false && M_2 > M_1
            continue
        end

        %%%%%%%%%%%%%%%%%% Valores apos a terceira onda de choque (normal)
        % Angulo da onda de choque (normal)
        b_3 = 90;
        theta_3 = atand ((2 * cotd (b_3) * (((M_2^2) *sind (b_3)^2) -1)) ...
        /((M_2^2)*(gamma+cosd(2*b_3))+2));
        Mn_2 = M_2*sind(b_3);

        % Verificar se Mn_2 real
        if isreal (Mn_2) == false
            continue
        end

        % Mach normal
        Mn_3 = sqrt((Mn_2^2 + (2 / (gamma - 1)))/((((2 * gamma) / (gamma - 1)) * Mn_2^2)-1));
        
        % Verificar se Mn_3 real
        if isreal (Mn_3) == false
            continue
        end

        % Numero de Mach
        M_3 = Mn_3 / sind (b_3 - theta_3);

        % Verificar se M_3 real e menor que M_2
        if isreal (M_3)==false && M_3 > M_2
            continue
        end

        % Razao de massa especifica
        rho_r3 = ((gamma + 1) * Mn_2^2)/ (2 + ((gamma - 1) * Mn_2^2));

        % Massa especifica
        rho_3 = rho_2*rho_r3;

        % Razao de pressoes
        p_r3 = 1+(((2*gamma)/(gamma+1))*((Mn_2^2)-1));

        % Pressao
        p_3 = p_2*p_r3;

        % Razao de temperatura
        T_r3 = p_r3 / rho_r3;

        % Temperatura
        T_3 = T_2*T_r3;
        deltas_3 = cp*log(T_r3)-(R*log(p_r3));

        % Eficiencia depois da onda de choque normal
        eff_3 = exp((-(deltas_3))/R);

        % Velocidade do som
        a_3 = (gamma*R*T_3)^(1/2);

        % Velocidade do escoamento
        v_3 = M_3*a_3;

        % Eficiência total
        global_eff = eff_1 * eff_2 * eff_3;

        % Verificação para o gráfico de superfície
        if isreal (theta_1) == true && isreal (theta_2) == true && isreal (global_eff) == true && M_3<M_2 && M_2 <M_1
            i = i+1;
            thetal_v(i) = theta_1;
            theta2_v(i) = theta_2;
            eff_v(i) = global_eff;
        end

        %%%%%%%%%%%%%%%%%% Otimiza
        if global_eff > max_eff
            if M_1 < mach_inicial && M_2 < M_1 && M_3 < M_2 && M_3 <= 1
                max_eff = global_eff; % Efici neia total
                deltas_2_max = deltas_2;
                angulo_theta = [theta_1 theta_2 theta_3]; % Theta optimum
                theta_1_opt = theta_1; theta_2_opt = theta_2;
                theta_3_opt = theta_3;
                angulo_beta = [b_1 b_2 b_3]; % Beta optimum
                betal_opt = b_1; beta2_opt = b_2; beta3_opt = b_3;
                numero_mach = [mach_inicial M_1 M_2 M_3]; % Mach optimum
                mach_inicial_opt = mach_inicial; M1_opt = M_1; M2_opt = M_2; M3_opt = M_3;
                rho = [densidade_inicial rho_1 rho_2 rho_3]; % Massa espec fica (optimum)
                densidade_inicial_opt = densidade_inicial; rhol_opt = rho_1; rho2_opt = rho_2; rho3_opt = rho_3;
                pressure = [pressao_inicial p_1 p_2 p_3]; % Presso (optimum)
                pinf_opt = pressao_inicial; p1_opt = p_1; p2_opt = p_2; p3_opt = p_3;
                temperature = [temperatura_inicial T_1 T_2 T_3]; % Temperatura (optimum)
                tempinf_opt = temperatura_inicial; T1_opt = T_1; T2_opt = T_2; T3_opt = T_3;
                fluid_velocity = [veltotal_inicial v_1 V_2 v_3]; % Velocidade do escoamento (optimum)
                vinf_opt = veltotal_inicial; v1_opt = v_1; v2_opt = V_2; v3_opt = v_3;
                velocity_sound = [velsom_inicial a_1 a_2 a_3]; % Velocidade do som (optimum)
                ainf_opt = velsom_inicial; a1_opt = a_1; a2_opt = a_2; a3_opt = a_3;
                PO = [0 eff_1 eff_2 eff_3];
                D = ((p_2-p_1)*L_1*sind (theta_1))+((p_3-p_1) *L_2 *sind(theta_1_opt+theta_2)); % Arrasto [N]
            end
        end
    end
end

%%%%%%%%%%%%%%%% Resultados
%% Tabela com resultados otimizados
T = table({'Infinito'; 'Estacao 1'; 'Estacao 2'; 'Estacao 3'}, ...
    [mach_inicial_opt; M1_opt; M2_opt; M3_opt], ...
    [pinf_opt; p1_opt; p2_opt; p3_opt], ...
    [tempinf_opt; T1_opt; T2_opt; T3_opt], ...
    [vinf_opt; v1_opt; v2_opt; v3_opt], ...
    [densidade_inicial_opt; rhol_opt; rho2_opt; rho3_opt], ...
    'VariableNames', {'Estacao', 'Mach', 'Pressao [Pa]', 'Temperatura [K]', 'Velocidade [m/s]', 'Massa especifica [Kg/m3]'});
disp('~~~~~~ Parte 1 ~~~~~~')
disp('Propriedades Difusor:')
disp(T)
fprintf('Theta (graus): Theta_1:%.3f Theta_2:%.3f\n', theta_1_opt, theta_2_opt);
fprintf('Beta (graus): Beta_1:%.3f; Beta_2:%.3f\n', betal_opt, beta2_opt);
fprintf('Eficiencia [%%]: %.3f \n', max_eff*100);
fprintf('Arrasto [N]: %.3f \n', D);

% Cria o do gr fico de efic ncia
[Xq, Yq] = meshgrid (linspace (min (thetal_v), max (thetal_v), 100), linspace (min (theta2_v), max(theta2_v), 100));
Zq = griddata (thetal_v, theta2_v, eff_v, Xq, Yq, 'cubic'); % interpola
max_value = max(eff_v);
figure (1)
hold on
contour (Xq, Yq, Zq, 20);
xlabel('\theta_1 (graus)');
ylabel('\theta_2 (graus)');
colormap (jet);
cbar = colorbar;
cbar.Label.String = '\eta';
cbar.Label.FontSize = 14;
max_value_formatted = sprintf('%.3f', max_value);
text (theta_1_opt, theta_2_opt, ['\eta_{max} = ' max_value_formatted], 'Color', 'black', 'FontSize', 12);
grid on;
hold off

%% Gráfico de todas as propriedades
x = [1 2 3 4];
figure (2)
plot(x, numero_mach, '-o', 'Color', [0, 0, 0]);
names = {'Estacao Inf'; 'Estacao 1'; 'Estacao 2'; 'Estacao 3'};
set (gca, 'xtick', 1:4, 'xticklabel', names)
grid on

% Adicionando os eixos de P, T, rho, M e V
addaxis (x, fluid_velocity, '-*', 'Color', [0.8500, 0.3250, 0.0980]);
addaxis(x, pressure, '-s', 'Color', [0.4940, 0.1840, 0.5560]);
addaxis (x, rho, '-d', 'Color', [0 0.4470 0.7410]);
addaxis(x, temperature, '-^', 'Color', [0, 0.3922, 0]);
addaxislabel (1, 'Mach');
addaxislabel (2, 'Velocidade [m/s]');
addaxislabel (3, 'Pressao [Pa]');
addaxislabel (4, 'Massa especifica [kg/m3]');
addaxislabel (5, 'Temperatura [K]');

% Deixar o fundo da legenda trasnparente
legend ('Mach', 'Velocidade', 'Pressao', 'Massa especifica', 'Temperatura');
lgd = legend;
set (lgd, 'Color', 'none');

% Tornar o gráfico quadrado
axis square;
set (gcf, 'Position', [100, 100, 800, 500]);