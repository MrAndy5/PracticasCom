function SNR_obj = SNR_para_BER(SNR_dB, BER, BER_obj)
% SNR (dB) a la que la curva BER cruza BER_obj, interpolando en log10(BER).
% Devuelve NaN si la curva no llega a cruzar BER_obj en el rango simulado.
    SNR_obj = NaN;
    v = BER > 0 & ~isnan(BER);
    x = SNR_dB(v);
    y = log10(BER(v));
    i = find(y < log10(BER_obj), 1, 'first');
    if isempty(i) || i == 1
        return
    end
    SNR_obj = x(i-1) + (log10(BER_obj) - y(i-1)) * (x(i) - x(i-1)) / (y(i) - y(i-1));
end