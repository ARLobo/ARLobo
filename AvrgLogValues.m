% Este array tem que ser populado com os valores de BS por sector angular determinado pela resolução angular definida pelo utilizador (angRes)
% Array de valores em dB a calcular a média
db_values = [-19.9, -19.9, -19.9];

% Conversão dos dBs numa escala de valores lineares
linear_values = 10.^(db_values/20);
%linear_values = linear_values*10;

% Calculo da média dos valores lineares
average_linear = mean(linear_values);
disp(['Average of linear values: ', num2str(average_linear)]);

% Conversão da média dos valores lineares novamente em dBs
average_db = 20 * log10(average_linear(1,:));
disp(['Average of dB values: ', num2str(average_db)]);

