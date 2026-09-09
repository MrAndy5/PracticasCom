function [phi1, phi2, a1, a2] = calcular_bases(s1, s2, Ts)
    % La energía es la integral discreta: sum(s.^2) * Ts
    %a1 = 1 / sqrt(sum(s1.^2) * Ts);
    %a2 = 1 / sqrt(sum(s2.^2) * Ts);
    
    phi1 = s1 / sqrt(Ts);
    phi2 = s2 / sqrt(Ts);
end