function [y1, y2] = demoduladorp12(r_t, phi1, phi2, Ts)
% DEMODULADORP12 Demodulador por correlación sobre bases phi1 y phi2 para N símbolos (Parte 2).
%
% ENTRADAS:
%   r_t        - Señal recibida en el tiempo (longitud N * L).
%   phi1, phi2 - Vectores de las funciones base ortonormales (longitud L).
%   Ts         - Periodo de muestreo en segundos.
%
% SALIDAS:
%   y1, y2     - Muestras de las ramas 1 y 2 en t = k*T (vectores 1 x N).
%   (Si se pide una sola salida: y = [y1; y2] matriz de 2 x N).

    L = length(phi1);
    N = length(r_t) / L;
    
    y1 = zeros(1, N);
    y2 = zeros(1, N);
    
    % Correlación e integración para cada uno de los N periodos de símbolo
    for k = 1:N
        idx_ini = (k - 1) * L + 1;
        idx_fin = k * L;
        
        r_k = r_t(idx_ini:idx_fin);
        
        % Integral discreta con las bases ortonormales
        y1(k) = sum(r_k .* phi1) * Ts;
        y2(k) = sum(r_k .* phi2) * Ts;
    end
    
    % Permite llamar tanto como [y1, y2] = demoduladorp12(...)
    % como y = demoduladorp12(...) (matriz 2 x N para el Detector)
    if nargout <= 1
        y1 = [y1; y2];
    end
end