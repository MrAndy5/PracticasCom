function [y1, y2] = correlatorNoise(T, Ts, r, sigma)
% CORRELATORNOISE Aplica ruido WGN (media 0) a la señal y la demodula por correlación.
%
% ENTRADAS:
%   T     - Tiempo de símbolo (10 ms).
%   Ts    - Periodo de muestreo (0.5 ms).
%   r     - Señal pura enviada (L muestras).
%   sigma - Desviación típica del ruido WGN (media = 0, varianza = sigma^2).
%
% SALIDAS:
%   y1, y2 - Salidas temporales de los correladores 1 y 2.

    if nargin < 4 || isempty(sigma)
        sigma = 0.5; % Valor por defecto de sigma
    end

    % 1. Generación del ruido blanco gaussiano (WGN) con media mu = 0
    ruido = sigma * randn(size(r));
    r_contaminada = r + ruido;

    % 2. Reutilizamos correlatorType pasándole la señal con ruido
    [y1, y2] = correlatorType(T, Ts, r_contaminada);
end