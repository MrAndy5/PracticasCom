function [s_t, s] = modulador(s_in, phi1, phi2)
% MODULADOR Genera la señal modulada continua en banda base para M=4 símbolos.
% Reutiliza las bases ortonormales phi1 y phi2 obtenidas en la P0/P1.
%
% ENTRADAS:
%   s_in - Vector con la secuencia de símbolos (1 a 4) O escalar N con el número de símbolos.
%   phi1 - Función base 1 (obtenida de s1/sqrt(T)).
%   phi2 - Función base 2 (obtenida de s2/sqrt(T)).
%
% SALIDAS:
%   s_t  - Señal modulada concatenada (1 x (N*L) muestras).
%   s    - Vector de símbolos transmitidos (1 a 4).

    % Si solo se pasa un número N, generamos N símbolos aleatorios equiprobables
    if isscalar(s_in)
        N = s_in;
        s = randi([1, 4], 1, N);
    else
        s = s_in;
        N = length(s);
    end

    % Coeficientes de la base c_ij calculados teóricamente
    % s_i(t) = c(i,1)*phi1 + c(i,2)*phi2
    c = [ 0.1,    0;   ... % s1
            0,  0.1;   ... % s2
         -0.1,    0;   ... % s3
            0, -0.1];      % s4

    % Construimos las 4 formas de onda a partir de la base
    ondas = cell(1, 4);
    for i = 1:4
        ondas{i} = c(i, 1) * phi1 + c(i, 2) * phi2;
    end

    % Concatenamos los N símbolos
    s_t = cell2mat(ondas(s));
end