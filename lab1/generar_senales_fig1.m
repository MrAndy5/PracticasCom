function [s1, s2, t] = generar_senales_fig1(T, Ts)
    L = T / Ts;                             % Número de muestras (20)
    t = (0:L-1) * Ts;                       % Vector de tiempo de un símbolo
    
    s1 = ones(1, L);                        % Pulso constante a 1
    s2 = [ones(1, L/2), -ones(1, L/2)];     % Primera mitad a 1, segunda a -1
end