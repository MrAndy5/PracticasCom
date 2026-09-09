function [s_tx, t_tx] = modular_secuencia(X, s1, s2, Ts)
    s_tx = [];
    for x = X
        if x == 1
            s_tx = [s_tx, s1];
        else
            s_tx = [s_tx, s2];
        end
    end
    t_tx = (0:length(s_tx)-1) * Ts;
end