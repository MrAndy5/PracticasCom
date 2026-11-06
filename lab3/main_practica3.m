%% =========================================================================
%  PRÁCTICA 3: MODULACIONES DIGITALES (QPSK, DQPSK, M-QAM, APSK)
%  Laboratorio de Sistemas de Comunicación I - Máster ICAI
%  =========================================================================
%  Script modular completo con métricas avanzadas:
%  - Representación de constelaciones y regiones de Voronoi
%  - Análisis de ruido AWGN, EVM y probabilidades teóricas
%  - Curvas de BER frente a SNR (dB) y distancias al Límite de Shannon
%  - Inmunidad y penalización angular de fase en QPSK vs DQPSK
%  - Comparativa de modulación multinivel M-QAM frente a APSK y PAPR
% =========================================================================

clear; clc; close all;

% Configuración de rutas y guardado
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
cd(scriptDir);
figDir = fullfile(scriptDir, 'figuras');
if ~exist(figDir, 'dir'), mkdir(figDir); end

fprintf('===================================================================\n');
fprintf('  SISTEMAS DE COMUNICACIÓN I - PRÁCTICA 3: MODULACIONES DIGITALES  \n');
fprintf('===================================================================\n\n');

%% =========================================================================
%% EJERCICIO 2.1: Modulación QPSK y Caracterización de Constelación
%% =========================================================================
fprintf('--- [EJERCICIO 2.1] Modulación QPSK y Mapeo Gray ---\n');

parejas_bits = [1 1; 1 0; 0 1; 0 0];
txBits_21 = reshape(parejas_bits', 1, []); 
xMod_21 = moduladorQPSK(txBits_21);

% 1. Tabla de mapeo y símbolos
Bits_Transmitidos = ["1 1"; "1 0"; "0 1"; "0 0"];
I_Real            = real(xMod_21(:));
Q_Imag            = imag(xMod_21(:));
Modulo            = abs(xMod_21(:));
Fase_deg          = mod(angle(xMod_21(:)) * 180 / pi, 360);

tabla_ej21 = table(Bits_Transmitidos, I_Real, Q_Imag, Modulo, Fase_deg, ...
    'VariableNames', {'Bits', 'I', 'Q', 'Modulo', 'Fase_deg'});
disp(tabla_ej21);

% 2. Métricas de Capa Física y Espacio de Señales
d_min       = 2.0;                                   % Distancia euclídea mínima entre símbolos contiguos
Es          = mean(abs(xMod_21).^2);                 % Energía media de símbolo (Es = 2)
Eb          = Es / 2;                               % Energía media por bit (Eb = 1 -> 0 dB)
d_min_norm  = d_min / sqrt(Eb);                     % Distancia normalizada d_min / sqrt(Eb) = 2
PAPR_dB     = 10*log10(max(abs(xMod_21).^2) / Es);   % PAPR de la constelación (0 dB -> envolvente constante)
Eficiencia  = 2;                                    % bits/s/Hz teórica en banda base (Nyquist alpha=0)

tabla_metricas_21 = table(d_min, d_min_norm, Es, Eb, PAPR_dB, Eficiencia, ...
    'VariableNames', {'d_min', 'd_min_norm_Eb', 'Es', 'Eb', 'PAPR_dB', 'Eficiencia_bps_Hz'});
disp(tabla_metricas_21);

% 3. Gráfica de constelación en plano I-Q
figure('Color', 'w', 'Position', [150 150 550 500]);
ax = gca;
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.85 0.85 0.85]);
plot(real(xMod_21), imag(xMod_21), 'b*', 'MarkerSize', 10, 'LineWidth', 1.8);
hold on; grid on; axis square;
axis([-2.2 2.2 -2.2 2.2]);

title('Constelación QPSK y Mapeo Gray (2 bits/símbolo)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('En Fase, I', 'FontSize', 11);
ylabel('Cuadratura, Q', 'FontSize', 11);

xline(0, 'k--', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 1);
th = linspace(0, 2*pi, 200);
plot(sqrt(2)*cos(th), sqrt(2)*sin(th), 'r:', 'LineWidth', 1.2);

etiquetas = {'(11) 45°', '(10) 315°', '(01) 135°', '(00) 225°'};
offset_x = [0.12,  0.12, -0.55, -0.55];
offset_y = [0.12, -0.18,  0.12, -0.18];
for idx = 1:length(xMod_21)
    text(real(xMod_21(idx)) + offset_x(idx), imag(xMod_21(idx)) + offset_y(idx), ...
        etiquetas{idx}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0 0.2 0.7]);
end

legend({'Símbolos QPSK', 'Frontera I = 0', 'Frontera Q = 0', 'Envolvente r = \surd 2'}, ...
    'Location', 'northeastoutside', 'Color', 'w');
saveas(gcf, fullfile(figDir, 'ejercicio2_1_constelacion_qpsk.png'));

%% =========================================================================
%% EJERCICIO 2.2: QPSK en Canal AWGN con Métrica EVM
%% =========================================================================
fprintf('\n--- [EJERCICIO 2.2] Transmisión QPSK sobre Canal AWGN ---\n');

Nbits = 4000;
rng(42);
txBits_rand = randi([0 1], 1, Nbits);
xMod_rand = moduladorQPSK(txBits_rand);

k = log2(4);
SNR_vect_dB = [5, 10, 15];
simb_ideales = [1+1j, -1+1j, -1-1j, 1-1j];

EbN0_dB_calc    = (SNR_vect_dB - 10*log10(k))';
EbN0_lin        = 10.^(EbN0_dB_calc / 10);
sigma_ruido     = sqrt(10.^(-SNR_vect_dB(:) / 10));
margen_frontera = 1.0 ./ sigma_ruido;

% Métricas de EVM y Probabilidad de Error
EVM_porc        = (1 ./ sqrt(10.^(SNR_vect_dB(:)/10))) * 100;
BER_teorica     = qfunc(sqrt(2 * EbN0_lin));
SER_teorica     = 2*BER_teorica - BER_teorica.^2;

tabla_ej22 = table(SNR_vect_dB', EbN0_dB_calc, sigma_ruido, margen_frontera, ...
    EVM_porc, BER_teorica, SER_teorica, ...
    'VariableNames', {'SNR_dB', 'EbN0_dB', 'Sigma_Ruido', 'Margen_d_sigma', ...
                      'EVM_porc', 'BER_Teorica', 'SER_Teorica'});
disp(tabla_ej22);

figure('Color', 'w', 'Position', [80 200 1350 420]);
for i = 1:length(SNR_vect_dB)
    snr_i = SNR_vect_dB(i);
    rx_noisy = awgn(xMod_rand, snr_i, 'measured');
    
    subplot(1, 3, i);
    ax = gca;
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.85 0.85 0.85]);
    
    plot(real(rx_noisy), imag(rx_noisy), '.', 'Color', [0.2 0.45 0.85], 'MarkerSize', 5);
    hold on; grid on; axis square;
    axis([-3.5 3.5 -3.5 3.5]);
    
    xline(0, 'k--', 'LineWidth', 1.1);
    yline(0, 'k--', 'LineWidth', 1.1);
    plot(real(simb_ideales), imag(simb_ideales), 'ro', 'MarkerSize', 7, 'LineWidth', 2);
    
    title(sprintf('SNR = %d dB (E_b/N_0 = %.2f dB) | EVM = %.1f%%', ...
        snr_i, EbN0_dB_calc(i), EVM_porc(i)), 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('En Fase, I'); ylabel('Cuadratura, Q');
    
    if i == 1
        legend({'Muestras recibidas', 'Fronteras de decisión', '', 'Constelación nominal'}, ...
            'Location', 'northwest', 'FontSize', 8, 'Color', 'w');
    end
end
sgtitle('Dispersión de la Constelación QPSK y Evolución de EVM frente a SNR', ...
    'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, fullfile(figDir, 'ejercicio2_2_dispersion_qpsk_evm.png'));

%% =========================================================================
%% EJERCICIO 2.3: Cadena Completa QPSK vs DQPSK (Canal Ideal)
%% =========================================================================
fprintf('\n--- [EJERCICIO 2.3] Cadena Completa QPSK vs DQPSK en Canal Ideal ---\n');

% Rama QPSK Coherente
xMod_QPSK   = moduladorQPSK(txBits_rand);
rxBits_QPSK = demoduladorQPSK(xMod_QPSK);
numErr_QPSK = sum(abs(txBits_rand - rxBits_QPSK));
BER_QPSK    = numErr_QPSK / length(txBits_rand);

% Rama DQPSK Diferencial
xMod_DQPSK   = moduladorDQPSK(txBits_rand);
rxBits_DQPSK = demoduladorDQPSK(xMod_DQPSK);
numErr_DQPSK = sum(abs(txBits_rand - rxBits_DQPSK));
BER_DQPSK    = numErr_DQPSK / length(txBits_rand);

Esquema           = ["QPSK Coherente"; "DQPSK Diferencial"];
Bits_Tx           = [length(txBits_rand); length(txBits_rand)];
Errores_Bit       = [numErr_QPSK; numErr_DQPSK];
BER_Obtenida      = [BER_QPSK; BER_DQPSK];
Condicion_Canal   = ["Canal ideal (sin ruido)"; "Canal ideal (sin ruido)"];
Sincronizacion    = ["Bucle de Costas / PLL obligatorio"; "No coherente (detección de fase relativa)"];

tabla_ej23 = table(Esquema, Bits_Tx, Errores_Bit, BER_Obtenida, Condicion_Canal, Sincronizacion, ...
    'VariableNames', {'Esquema', 'Bits_Tx', 'Errores', 'BER', 'Canal', 'Requisito_Sincronismo'});
disp(tabla_ej23);

%% =========================================================================
%% EJERCICIO 2.4: BER frente a SNR para QPSK y DQPSK
%% =========================================================================
fprintf('\n--- [EJERCICIO 2.4] Curvas de BER frente a SNR en Canal AWGN ---\n');

EbN0_dB_vect = -5:2:11;                % Incluye 11 dB para capturar el cruce por BER = 1e-3 en DQPSK
Nbits_BER    = 4e5;

[BER_sim_QPSK,  BER_teor_QPSK, SNR_dB_vect] = simular_BER_PSK('QPSK',  EbN0_dB_vect, Nbits_BER, 0);
[BER_sim_DQPSK, BER_teor_DQPSK]             = simular_BER_PSK('DQPSK', EbN0_dB_vect, Nbits_BER, 0);

% Curvas teóricas con resolución fina
SNR_fino  = SNR_dB_vect(1):0.1:SNR_dB_vect(end);
EbN0_fino = 10.^((SNR_fino - 10*log10(2)) / 10);
teor_QPSK_fino  = qfunc(sqrt(2*EbN0_fino));
teor_DQPSK_fino = 1.13*qfunc(sqrt(1.2*EbN0_fino));
teor_QPSK_fino(teor_QPSK_fino < 1e-5)   = NaN;
teor_DQPSK_fino(teor_DQPSK_fino < 1e-5) = NaN;

% Ratio blindado contra división por cero
Ratio_BER = BER_sim_DQPSK ./ max(BER_sim_QPSK, 1/Nbits_BER);

tabla_ej24 = table(SNR_dB_vect', EbN0_dB_vect', BER_sim_QPSK', BER_teor_QPSK', ...
    BER_sim_DQPSK', BER_teor_DQPSK', Ratio_BER', ...
    'VariableNames', {'SNR_dB', 'EbN0_dB', 'BER_QPSK_sim', 'BER_QPSK_teor', ...
                      'BER_DQPSK_sim', 'BER_DQPSK_teor', 'BER_DQPSK_div_QPSK'});
disp(tabla_ej24);

% Interpolación de SNR requerida y penalización de DQPSK
BER_obj   = [1e-1; 1e-2; 1e-3];
SNR_QPSK  = arrayfun(@(b) SNR_para_BER(SNR_dB_vect, BER_sim_QPSK,  b), BER_obj);
SNR_DQPSK = arrayfun(@(b) SNR_para_BER(SNR_dB_vect, BER_sim_DQPSK, b), BER_obj);
tabla_SNR_24 = table(BER_obj, SNR_QPSK, SNR_DQPSK, SNR_DQPSK - SNR_QPSK, ...
    'VariableNames', {'BER_objetivo', 'SNR_QPSK_dB', 'SNR_DQPSK_dB', 'Penalizacion_DQPSK_dB'});
disp(tabla_SNR_24);

% Gráfica BER vs SNR
figure('Color', 'w', 'Position', [150 150 750 520]);
semilogy(SNR_dB_vect, BER_sim_QPSK,  'bo', 'MarkerSize', 7, 'LineWidth', 1.6); hold on;
semilogy(SNR_fino,    teor_QPSK_fino,  'b-', 'LineWidth', 1.3);
semilogy(SNR_dB_vect, BER_sim_DQPSK, 'rs', 'MarkerSize', 7, 'LineWidth', 1.6);
semilogy(SNR_fino,    teor_DQPSK_fino, 'r-', 'LineWidth', 1.3);
for b = BER_obj'
    yline(b, 'k:', 'HandleVisibility', 'off');
end
grid on; ylim([1e-5 1]);
title('BER frente a SNR en canal AWGN: QPSK vs DQPSK', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('SNR (dB) [Nota: E_b/N_0 = SNR - 3.01 dB]', 'FontSize', 11);
ylabel('BER (Escala Logarítmica)', 'FontSize', 11);
legend({'QPSK simulada', 'QPSK teórica', 'DQPSK simulada', 'DQPSK teórica'}, ...
    'Location', 'southwest', 'Color', 'w');
saveas(gcf, fullfile(figDir, 'ejercicio2_4_ber_vs_snr_qpsk_dqpsk.png'));

% Gráfica del cociente de BER
figure('Color', 'w', 'Position', [150 150 650 420]);
plot(SNR_dB_vect, Ratio_BER, 'ko-', 'LineWidth', 1.6, 'MarkerSize', 7);
grid on;
title('Cociente BER_{DQPSK} / BER_{QPSK} frente a SNR', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'FontSize', 11);
ylabel('BER_{DQPSK} / BER_{QPSK}', 'FontSize', 11);
yline(2.0, 'r--', 'Factor asintótico de propagación \approx 2', 'LineWidth', 1.2);
saveas(gcf, fullfile(figDir, 'ejercicio2_4_ratio_ber.png'));

%% =========================================================================
%% EJERCICIO 2.5: Influencia del Error de Fase en Recepción
%% =========================================================================
fprintf('\n--- [EJERCICIO 2.5] Error de Fase del Oscilador Local (0°, 10°, 20°, 30°) ---\n');

offsets_deg = [0 10 20 30];
nOff        = length(offsets_deg);
EbN0_dB_25  = -5:1:14;
Nbits_25    = 5e5;                     % Optimizado para precisión y velocidad
colores     = lines(nOff);

BER_rot_QPSK  = zeros(nOff, length(EbN0_dB_25));
BER_rot_DQPSK = zeros(nOff, length(EbN0_dB_25));
BER_rot_teor  = zeros(nOff, length(EbN0_dB_25));
EbN0_lin_25   = 10.^(EbN0_dB_25 / 10);

for i = 1:nOff
    phi = offsets_deg(i);
    [BER_rot_QPSK(i,:), ~, SNR_dB_25] = simular_BER_PSK('QPSK',  EbN0_dB_25, Nbits_25, phi);
    BER_rot_DQPSK(i,:)                = simular_BER_PSK('DQPSK', EbN0_dB_25, Nbits_25, phi);
    
    % Expresión analítica exacta de QPSK con rotación angular constante
    BER_rot_teor(i,:) = 0.5*( qfunc(sqrt(2*EbN0_lin_25)*(cosd(phi) - sind(phi))) + ...
                              qfunc(sqrt(2*EbN0_lin_25)*(cosd(phi) + sind(phi))) );
end
BER_rot_teor(BER_rot_teor < 1e-5) = NaN;

% 1. Tablas de SNR necesaria y penalización
BER_obj_25 = [1e-2 1e-3 1e-4];
SNRreq_Q = zeros(nOff, 3);
SNRreq_D = zeros(nOff, 3);
for i = 1:nOff
    for j = 1:3
        SNRreq_Q(i,j) = SNR_para_BER(SNR_dB_25, BER_rot_QPSK(i,:),  BER_obj_25(j));
        SNRreq_D(i,j) = SNR_para_BER(SNR_dB_25, BER_rot_DQPSK(i,:), BER_obj_25(j));
    end
end
Perdida_asint_dB = -20*log10(cosd(offsets_deg') - sind(offsets_deg'));

tabla_ej25_QPSK = array2table([offsets_deg', SNRreq_Q, SNRreq_Q - SNRreq_Q(1,:), Perdida_asint_dB], ...
    'VariableNames', {'phi_deg', 'SNR_1e2', 'SNR_1e3', 'SNR_1e4', ...
    'Penal_1e2', 'Penal_1e3', 'Penal_1e4', 'Penal_asint_dB'});
tabla_ej25_DQPSK = array2table([offsets_deg', SNRreq_D, SNRreq_D - SNRreq_D(1,:)], ...
    'VariableNames', {'phi_deg', 'SNR_1e2', 'SNR_1e3', 'SNR_1e4', ...
    'Penal_1e2', 'Penal_1e3', 'Penal_1e4'});

disp('QPSK (Penalización severa por desajuste de fase):'); disp(tabla_ej25_QPSK);
disp('DQPSK (Inmunidad a error de fase estático):'); disp(tabla_ej25_DQPSK);

% 2. Constelaciones recibidas con error de fase (SNR = 15 dB)
SNR_const  = 15;
bits_const = randi([0 1], 1, 4000);
xQ_c = moduladorQPSK(bits_const);
xD_c = moduladorDQPSK(bits_const);
figure('Color', 'w', 'Position', [50 80 1400 650]);
for i = 1:nOff
    rot = exp(1j*offsets_deg(i)*pi/180);
    rQ  = awgn(xQ_c, SNR_const, 'measured') * rot;
    rD  = awgn(xD_c, SNR_const, 'measured') * rot;
    zD  = rD(2:end) ./ exp(1j*angle(rD(1:end-1)));

    subplot(2, nOff, i);
    plot(real(rQ), imag(rQ), '.', 'Color', [0.2 0.45 0.85], 'MarkerSize', 4); hold on;
    plot(real(simb_ideales), imag(simb_ideales), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    xline(0, 'k--'); yline(0, 'k--');
    grid on; axis square; axis([-2.5 2.5 -2.5 2.5]);
    title(sprintf('QPSK, \\phi = %d°', offsets_deg(i)), 'FontWeight', 'bold');
    xlabel('I'); ylabel('Q');

    subplot(2, nOff, nOff + i);
    plot(real(zD), imag(zD), '.', 'Color', [0.85 0.33 0.1], 'MarkerSize', 4); hold on;
    plot(real(simb_ideales), imag(simb_ideales), 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
    xline(0, 'k--'); yline(0, 'k--');
    grid on; axis square; axis([-2.5 2.5 -2.5 2.5]);
    title(sprintf('DQPSK (Demodulado), \\phi = %d°', offsets_deg(i)), 'FontWeight', 'bold');
    xlabel('I'); ylabel('Q');
end
sgtitle(sprintf('Constelaciones con Error de Fase en Canal (SNR = %d dB)', SNR_const), ...
    'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, fullfile(figDir, 'ejercicio2_5_constelaciones_rotadas.png'));

% 3. Comparación QPSK vs DQPSK para cada ángulo
figure('Color', 'w', 'Position', [40 200 1450 420]);
for i = 1:nOff
    subplot(1, nOff, i);
    semilogy(SNR_dB_25, BER_rot_QPSK(i,:),  'bo-', 'LineWidth', 1.4, 'MarkerSize', 5); hold on;
    semilogy(SNR_dB_25, BER_rot_DQPSK(i,:), 'rs-', 'LineWidth', 1.4, 'MarkerSize', 5);
    grid on; ylim([1e-5 1]);
    title(sprintf('\\phi = %d°', offsets_deg(i)), 'FontWeight', 'bold');
    xlabel('SNR (dB)'); ylabel('BER');
    legend({'QPSK', 'DQPSK'}, 'Location', 'southwest', 'Color', 'w');
end
sgtitle('Rendimiento QPSK vs DQPSK ante Distintos Errores de Fase', 'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, fullfile(figDir, 'ejercicio2_5_comparacion_angulos.png'));

% 4. Umbral de cruce angular (SNR fija = 10 dB)
SNR_fija      = 10;
theta_barrido = 0:2.5:45;
BER_th_Q      = zeros(size(theta_barrido));
BER_th_D      = zeros(size(theta_barrido));
for i = 1:length(theta_barrido)
    BER_th_Q(i) = simular_BER_PSK('QPSK',  SNR_fija - 10*log10(2), 1e5, theta_barrido(i));
    BER_th_D(i) = simular_BER_PSK('DQPSK', SNR_fija - 10*log10(2), 1e5, theta_barrido(i));
end
EbN0_fija = 10^((SNR_fija - 10*log10(2))/10);
BER_th_teor = 0.5*( qfunc(sqrt(2*EbN0_fija)*(cosd(theta_barrido) - sind(theta_barrido))) + ...
                    qfunc(sqrt(2*EbN0_fija)*(cosd(theta_barrido) + sind(theta_barrido))) );
theta_cruce_sim = interp1(log10(BER_th_Q) - log10(BER_th_D), theta_barrido, 0);

figure('Color', 'w', 'Position', [150 150 750 500]);
semilogy(theta_barrido, BER_th_Q,    'bo-', 'LineWidth', 1.5); hold on;
semilogy(theta_barrido, BER_th_teor, 'b--', 'LineWidth', 1.1);
semilogy(theta_barrido, BER_th_D,    'rs-', 'LineWidth', 1.5);
xline(theta_cruce_sim, 'k:', sprintf('Cruce \\phi = %.1f°', theta_cruce_sim), ...
    'LineWidth', 1.2, 'FontSize', 10);
grid on;
title(sprintf('BER frente al Error de Fase para SNR Fija = %d dB', SNR_fija), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Error de fase del Oscilador Local \phi (°)', 'FontSize', 11);
ylabel('BER', 'FontSize', 11);
legend({'QPSK simulada', 'QPSK teórica', 'DQPSK (Inmune)'}, 'Location', 'northwest', 'Color', 'w');
saveas(gcf, fullfile(figDir, 'ejercicio2_5_cruce_angular.png'));

%% =========================================================================
%% EJERCICIO 3.1: BER frente a SNR para M-QAM y APSK (Con Límite de Shannon)
%% =========================================================================
fprintf('\n--- [EJERCICIO 3.1] Comparativa M-QAM y APSK (Límite de Shannon y PAPR) ---\n');

EbN0_dB_3  = 0:2:24;
Modulacion = ["4-QAM"; "16-QAM"; "64-QAM"; "256-QAM"; "16-APSK"; "32-APSK"];
Bits_simb  = [2; 4; 6; 8; 4; 5];
M_QAM      = [4 16 64 256];
M_16APSK   = [4 12];      R_16APSK = [1 2.5];
M_32APSK   = [4 12 16];   R_32APSK = [1 2.5 4.3];

BER_3       = zeros(6, length(EbN0_dB_3));
BER_QAMteor = zeros(4, length(EbN0_dB_3));
for i = 1:4
    [BER_3(i,:), BER_QAMteor(i,:)] = BER_m_ary_QAM(M_QAM(i), EbN0_dB_3);
end
BER_3(5,:) = BER_APSK(M_16APSK, R_16APSK, EbN0_dB_3);
BER_3(6,:) = BER_APSK(M_32APSK, R_32APSK, EbN0_dB_3);
SNR_3      = EbN0_dB_3 + 10*log10(Bits_simb);

% 1. Métricas: Distancia normalizada, PAPR y Límite de Shannon
consts = {qammod((0:3)', 4), qammod((0:15)', 16), qammod((0:63)', 64), qammod((0:255)', 256), ...
          apskmod((0:15)', M_16APSK, R_16APSK), apskmod((0:31)', M_32APSK, R_32APSK)};
SNR_1e2       = zeros(6,1);
SNR_1e3       = zeros(6,1);
d_min_norm_Eb = zeros(6,1);
PAPR_dB       = zeros(6,1);

for i = 1:6
    c  = consts{i};
    Es = mean(abs(c).^2);
    D  = abs(c - c.');  D(D == 0) = Inf;
    d_min_norm_Eb(i) = min(D(:)) / sqrt(Es / Bits_simb(i));
    PAPR_dB(i)       = 10*log10(max(abs(c).^2) / Es);
    SNR_1e2(i)       = SNR_para_BER(SNR_3(i,:), BER_3(i,:), 1e-2);
    SNR_1e3(i)       = SNR_para_BER(SNR_3(i,:), BER_3(i,:), 1e-3);
end

% Capacidad de Shannon en canal AWGN: C/B = eta -> SNR_min = 2^eta - 1
SNR_Shannon_dB = 10 * log10(2.^Bits_simb - 1);
Shannon_Gap_dB = SNR_1e3 - SNR_Shannon_dB;

tabla_ej31 = table(Modulacion, Bits_simb, SNR_1e2, SNR_1e3, SNR_Shannon_dB, Shannon_Gap_dB, ...
    d_min_norm_Eb, PAPR_dB, ...
    'VariableNames', {'Modulacion', 'Bits_simbolo', 'SNR_BER1e2', 'SNR_BER1e3', ...
                      'Shannon_Limit_dB', 'Shannon_Gap_dB', 'd_min_norm_Eb', 'PAPR_dB'});
disp(tabla_ej31);

% 2. Constelaciones
figure('Color', 'w', 'Position', [40 150 1450 400]);
idx_c = [2 3 5 6];
th_c  = linspace(0, 2*pi, 200);
for j = 1:4
    i = idx_c(j);
    subplot(1, 4, j);
    plot(real(consts{i}), imag(consts{i}), 'b.', 'MarkerSize', 12); hold on;
    if i == 5, R = R_16APSK; elseif i == 6, R = R_32APSK; else, R = []; end
    for r = R
        plot(r*cos(th_c), r*sin(th_c), 'r:', 'LineWidth', 1);
    end
    grid on; axis square; axis equal;
    title(Modulacion(i), 'FontWeight', 'bold');
    xlabel('I'); ylabel('Q');
end
sgtitle('Constelaciones M-QAM y APSK (Concéntricas)', 'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, fullfile(figDir, 'ejercicio3_1_constelaciones.png'));

% 3. Curvas globales de BER frente a SNR
colores = lines(6);
marcas  = {'o-', 'o-', 'o-', 'o-', '^-', '^-'};
figure('Color', 'w', 'Position', [150 150 800 540]);
for i = 1:6
    semilogy(SNR_3(i,:), BER_3(i,:), marcas{i}, 'Color', colores(i,:), 'LineWidth', 1.5); hold on;
end
yline(1e-2, 'k:', 'HandleVisibility', 'off');
yline(1e-3, 'k:', 'HandleVisibility', 'off');
grid on; ylim([1e-5 1]);
title('BER frente a SNR: M-QAM y APSK', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('SNR (dB)', 'FontSize', 11); ylabel('BER', 'FontSize', 11);
legend(cellstr(Modulacion), 'Location', 'northeast', 'Color', 'w');
saveas(gcf, fullfile(figDir, 'ejercicio3_1_ber_todas.png'));

% 4. Compromiso Eficiencia vs SNR y PAPR
figure('Color', 'w', 'Position', [100 150 1200 450]);
subplot(1, 2, 1);
plot(Bits_simb(1:4), SNR_1e3(1:4), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(Bits_simb(5:6), SNR_1e3(5:6), 'r^-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Bits_simb, SNR_Shannon_dB, 'k--', 'LineWidth', 1.3);
text(Bits_simb + 0.15, SNR_1e3, cellstr(Modulacion));
grid on; xlim([1 9]);
title('SNR requerida para BER = 10^{-3} vs Límite de Shannon', 'FontWeight', 'bold');
xlabel('Eficiencia Espectral (bits/símbolo)'); ylabel('SNR (dB)');
legend({'M-QAM', 'APSK', 'Límite de Shannon'}, 'Location', 'northwest', 'Color', 'w');

subplot(1, 2, 2);
bar(categorical(Modulacion, Modulacion), PAPR_dB, 'FaceColor', [0.2 0.5 0.8]);
grid on;
title('Relación de Potencia Pico a Media (PAPR)', 'FontWeight', 'bold');
ylabel('PAPR (dB)');
saveas(gcf, fullfile(figDir, 'ejercicio3_1_eficiencia_papr.png'));

fprintf('\n¡Simulación completada con éxito! Todas las figuras han sido guardadas en %s\n', figDir);
