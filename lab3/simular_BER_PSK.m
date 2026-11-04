function [BER_sim, BER_teor, SNR_dB] = simular_BER_PSK(tipoMod, EbN0_dB_vect, Nbits, offsetFase_deg)
%% SIMULAR_BER_PSK Simulador Monte Carlo modular para modulaciones M-PSK en canal AWGN
%  Permite evaluar prestaciones de BER frente a Eb/N0 tanto en condiciones coherentes
%  nominales como ante desajustes de fase del oscilador local del receptor.
%
%  ENTRADAS:
%    tipoMod        : Cadena identificadora ('QPSK' o 'DQPSK')
%    EbN0_dB_vect   : Vector con valores de Eb/N0 en dB (ej. -5:2:10)
%    Nbits          : Número de bits a transmitir (debe ser par para 4-PSK)
%    offsetFase_deg : Rotación constante de fase introducida en canal (grados).
%                     Por defecto 0° (sin error de fase).
%
%  SALIDAS:
%    BER_sim        : Vector con la tasa de error de bit (BER) simulada
%    BER_teor       : Vector con la BER teórica en canal AWGN sin rotación
%    SNR_dB         : Vector con los valores de SNR en dB correspondientes
%
%  EJEMPLO DE USO:
%    [ber_qpsk, ber_th, snr] = simular_BER_PSK('QPSK', -5:2:10, 2e5, 0);
%    [ber_dqpsk, ~, ~]       = simular_BER_PSK('DQPSK', -5:2:10, 2e5, 15);
%
%  Máster en Ingeniería de Telecomunicación - Sistemas de Comunicación I

    if nargin < 4 || isempty(offsetFase_deg)
        offsetFase_deg = 0;
    end
    
    % Validación de longitud par
    if rem(Nbits, 2) ~= 0
        error('simular_BER_PSK: Nbits debe ser un número par.');
    end

    k = 2; % Bits por símbolo para QPSK / DQPSK (M = 4)
    numPuntos = length(EbN0_dB_vect);
    
    % Relación entre Eb/N0 y SNR para señal compleja en banda base
    SNR_dB   = EbN0_dB_vect + 10 * log10(k);
    BER_sim  = zeros(1, numPuntos);
    BER_teor = zeros(1, numPuntos);
    
    % Factor complejo de rotación de fase del canal
    rotacionCanal = exp(1j * (offsetFase_deg * pi / 180));
    
    % Generación de la secuencia de información (vector fila)
    txBits = randi([0 1], 1, Nbits);
    
    % Modulación en banda base
    switch upper(tipoMod)
        case 'QPSK'
            xMod = moduladorQPSK(txBits);
        case 'DQPSK'
            xMod = moduladorDQPSK(txBits);
        otherwise
            error('simular_BER_PSK: tipoMod no soportado. Use ''QPSK'' o ''DQPSK''.');
    end
    
    % Bucle de Monte Carlo sobre el vector de Eb/N0
    for idx = 1:numPuntos
        % Adición de AWGN midiendo la potencia de la señal transmitida
        rx_awgn = awgn(xMod, SNR_dB(idx), 'measured');
        
        % Introducción del offset de fase del oscilador local (canal)
        rx_rotada = rx_awgn * rotacionCanal;
        
        % Demodulación
        switch upper(tipoMod)
            case 'QPSK'
                rxBits = demoduladorQPSK(rx_rotada);
            case 'DQPSK'
                rxBits = demoduladorDQPSK(rx_rotada);
        end
        
        % Conteo de errores de bit y cálculo de BER empírica
        numErrores = sum(abs(txBits - rxBits));
        BER_sim(idx) = numErrores / Nbits;
    end
    
    % Cálculo analítico de BER teórica en canal AWGN no perturbado
    EbN0_lin = 10 .^ (EbN0_dB_vect / 10);
    switch upper(tipoMod)
        case 'QPSK'
            % QPSK con codificación Gray: Pb = Q(sqrt(2*Eb/N0))
            BER_teor = qfunc(sqrt(2 * EbN0_lin));
        case 'DQPSK'
            % Expresión teórica de referencia para DQPSK en AWGN:
            % Pb = 1.13 * Q(sqrt(1.2 * EbN0))
            BER_teor = 1.13 * qfunc(sqrt(1.2 * EbN0_lin));
    end
    
    % Criterio del enunciado para valores asintóticamente nulos
    BER_teor(BER_teor < 1e-5) = NaN;
end
