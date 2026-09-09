function [y1, y2] = correlatorType(T, Ts, r)
% CORRELATORTYPE Demodulador correlacionando con las bases ortonormales phi1 y phi2.
%
% ENTRADAS:
%   T  - Tiempo de símbolo (10 ms = 0.01 s).
%   Ts - Periodo de muestreo (0.5 ms).
%   r  - Señal recibida (L muestras).
%
% SALIDAS:
%   y1, y2 - Salidas temporales de los correladores de ambas ramas.

    % 1. Formas de onda base s1 y s2
    [s1, s2, ~] = generar_senales_fig1(T, Ts);
    
    % 2. Bases ortonormales phi = s / sqrt(T)
    phi1 = s1 / sqrt(T);
    phi2 = s2 / sqrt(T);
    
    % 3. Correlación acumulada en el tiempo (integral continua discretizada)
    y1 = cumsum(r .* phi1) * Ts;
    y2 = cumsum(r .* phi2) * Ts;
end