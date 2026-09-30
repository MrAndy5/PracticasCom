function s_det = detector(y1, y2, C)
% DETECTOR Detector de mínima distancia euclídea genérico para cualquier constelación.
%
% ENTRADAS:
%   y1, y2 - Muestras demoduladas a la salida de las ramas 1 y 2.
%   C      - Matriz de constelación (M x 2): cada fila i contiene [c_i1, c_i2].
%
% SALIDA:
%   s_det  - Vector con los índices de los símbolos detectados (1 a M).

    if nargin < 3 || isempty(C)
        % Matriz por defecto de la práctica si no se pasa
        C = [ 0.1,    0;   ... % s1
                0,  0.1;   ... % s2
             -0.1,    0;   ... % s3
                0, -0.1];      % s4
    end

    M = size(C, 1);
    N = length(y1);
    
    % Matriz de distancias al cuadrado (M x N)
    D = zeros(M, N);
    for i = 1:M
        D(i, :) = (y1 - C(i, 1)).^2 + (y2 - C(i, 2)).^2;
    end

    % Decisión de mínima distancia euclídea
    [~, s_det] = min(D, [], 1);
end