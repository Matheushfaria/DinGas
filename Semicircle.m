function [x,y] = Semicircle(x0, y0, raio, angulo_inicial, angulo_final)
    % Número de pontos para a discretização
    if y0==157.20e-3
    num_points = 400;
    else
    num_points = 50;
    end
    % Gere os ângulos para o semicírculo no sentido horário
    theta = linspace(angulo_final, angulo_inicial, num_points);  % Ângulos em ordem decrescente

    % Calcule as coordenadas dos pontos do semicírculo
    x = raio * cosd(theta);
    y = raio * sind(theta);

    % Aplicar o deslocamento
    x = x + x0;
    y = y + y0;
end
