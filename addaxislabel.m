function varargout = addaxislabel(axnum, label, fontSize)
% ADDAXISLABEL adds axis labels to axes made with ADDAXIS.m
% handle_to_text = addaxislabel(axis_number, label, fontSize);

% Verifique se o tamanho da fonte foi especificado
if nargin < 3
    fontSize = 12;  % Tamanho da fonte padrão
end

% Obtenha o identificador da figura atual
cah = gca;

% Obtenha os identificadores dos eixos adicionais criados com ADDAXIS
axh = getaddaxisdata(cah, 'axisdata');

% Obtenha os identificadores dos eixos e suas posições
axhand = cah;
postot(1, :) = get(cah, 'position');
for I = 1:length(axh)
    axhand(I+1) = axh{I}(1);
    postot(I+1, :) = get(axhand(I+1), 'position');
end

% Defina os eixos atuais para os eixos a serem rotulados
axes(axhand(axnum));
htxt = ylabel(label);
set(htxt, 'color', get(axhand(axnum), 'ycolor'));
set(htxt, 'FontSize', fontSize);  % Defina o tamanho da fonte

% Restaure os eixos atuais para os eixos principais
axes(cah);

if nargout == 1
    varargout{1} = htxt;
end
