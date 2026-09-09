function [r1, r2] = demodulador_muestreado(r_total, T, Ts)
% DEMODULADOR_MUESTREADO Muestrea al final de cada periodo sobre bases ortonormales.

    L = round(T / Ts);
    Nsymb = length(r_total) / L;
    
    [s1, s2, ~] = generar_senales_fig1(T, Ts);
    
    % Bases ortonormales phi = s / sqrt(T)
    phi1 = s1 / sqrt(T);
    phi2 = s2 / sqrt(T);
    
    r1 = zeros(1, Nsymb);
    r2 = zeros(1, Nsymb);
    
    for k = 1:Nsymb
        idx_ini = (k - 1) * L + 1;
        idx_fin = k * L;
        
        r_k = r_total(idx_ini:idx_fin);
        
        % Integral en [0, T] con la base
        r1(k) = sum(r_k .* phi1) * Ts;
        r2(k) = sum(r_k .* phi2) * Ts;
    end
end