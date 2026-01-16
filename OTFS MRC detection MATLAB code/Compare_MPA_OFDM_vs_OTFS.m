% =========================================================================
% Script: Compare_MPA_OFDM_vs_OTFS.m
% Mục đích: So sánh hiệu năng phát hiện MPA trên hai hệ thống OFDM và OTFS
%           trong môi trường kênh truyền có Doppler cao.
% =========================================================================

close all
clear all
rng(1) % Seed để tái lập kết quả

%% 1. Thiết lập tham số hệ thống
N = 8;             % Số ký hiệu thời gian (Giảm nhỏ để chạy MPA nhanh)
M = 16;            % Số sóng mang con (Giảm nhỏ để chạy MPA nhanh)
M_mod = 4;         % 4-QAM (QPSK)
M_bits = log2(M_mod);

% Thông số vật lý (Tạo kênh Doppler cao)
car_fre = 4*10^9;  % 4 GHz
delta_f = 15*10^3; % 15 kHz
T = 1/delta_f;
max_speed = 300;   % 300 km/h (Doppler mạnh gây nhiễu ICI cho OFDM)

% Cấu hình mô phỏng
SNR_dB = 0:5:25;   
N_fram = 100;       % Số frame mô phỏng
n_ite_MPA = 10;    % Số vòng lặp MPA

% Khởi tạo biến lưu kết quả
ber_otfs = zeros(1, length(SNR_dB));
ber_ofdm = zeros(1, length(SNR_dB));

% Ma trận DFT cho OFDM
Fn = dftmtx(M);       % DFT size M (cho sóng mang con)
Fn = Fn./norm(Fn);    % Chuẩn hóa
In = eye(N);          % Ma trận đơn vị size N

%% 2. Vòng lặp chính
fprintf('Bắt đầu mô phỏng so sánh MPA-OFDM và MPA-OTFS...\n');

for iesn0 = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(iesn0)/10);
    % Năng lượng trung bình của 4-QAM
    eng_sqrt = 1/sqrt(2)*sqrt(2); % Normalized
    sigma_2 = 1/SNR;
    
    total_err_otfs = 0;
    total_err_ofdm = 0;
    total_bits = 0;
    
    for ifram = 1:N_fram
        % --- A. Tạo dữ liệu ---
        N_syms = N*M;
        trans_bit = randi([0,1], N_syms*M_bits, 1);
        data = qammod(reshape(trans_bit, M_bits, N_syms), M_mod, 'gray', 'InputType', 'bit');
        
        % Mapping lên lưới (Full grid)
        x_in = reshape(data, M, N); % Dữ liệu đầu vào
        
        % --- B. Kênh truyền (Chung cho cả 2) ---
        [chan_coef, delay_taps, Doppler_taps, taps] = Generate_delay_Doppler_channel_parameters(N, M, car_fre, delta_f, T, max_speed);
        [G, ~] = Gen_time_domain_channel(N, M, taps, delay_taps, Doppler_taps, chan_coef);
        
        % Nhiễu Gaussian
        noise_vec = sqrt(sigma_2/2) * (randn(N*M, 1) + 1i*randn(N*M, 1));
        
        % =================================================================
        % HỆ THỐNG 1: OTFS + MPA
        % =================================================================
        % 1. Điều chế OTFS (DD -> Time)
        % Trong code gốc: X_tilda = X_DD * Fn_time' -> s = vec(X_tilda)
        % (Lưu ý: Code gốc dùng ISFFT theo chiều Doppler)
        Fn_time = dftmtx(N); Fn_time = Fn_time./norm(Fn_time);
        data_grid_otfs = x_in; % Coi input là miền Delay-Doppler
        
        % Modulation
        X_tilda = data_grid_otfs * Fn_time'; 
        s_otfs = reshape(X_tilda, N*M, 1);
        
        % Channel
        r_otfs = G * s_otfs + noise_vec;
        
        % Demodulation (Time -> DD)
        Y_tilda_otfs = reshape(r_otfs, M, N);
        Y_dd = Y_tilda_otfs * Fn_time;
        y_vec_otfs = reshape(Y_dd.', N*M, 1);
        
        % Tạo ma trận hiệu dụng H_OTFS (Delay-Doppler)
        [H_OTFS, ~, ~] = Gen_DD_and_DT_channel_matrices(N, M, G, Fn_time);
        
        % MPA Detection OTFS
        % Cần tạo data_grid logic để hàm MPA hiểu vị trí dữ liệu
        data_grid_logic = ones(M, N); 
        [est_bits_otfs, ~, ~] = MPA_detector(N, M, M_mod, sigma_2, data_grid_logic, y_vec_otfs, H_OTFS, n_ite_MPA);
        
        total_err_otfs = total_err_otfs + sum(xor(est_bits_otfs, trans_bit));
        
        % =================================================================
        % HỆ THỐNG 2: OFDM + MPA
        % =================================================================
        % 1. Điều chế OFDM (Freq-Time -> Time)
        % Input x_in coi là miền Tần số - Thời gian (Freq-Time)
        % Mỗi cột là 1 symbol OFDM -> IFFT trên từng cột (trục M)
        % s_ofdm = vec( IFFT(x_in) )
        
        X_time_ofdm = Fn' * x_in; % IFFT cột (Fn là DFT matrix size M)
        s_ofdm = reshape(X_time_ofdm, N*M, 1);
        
        % Channel
        r_ofdm = G * s_ofdm + noise_vec;
        
        % Demodulation (Time -> Freq-Time)
        % Y_freq = FFT( reshape(r) )
        R_mat = reshape(r_ofdm, M, N);
        Y_freq = Fn * R_mat; % FFT cột
        y_vec_ofdm = reshape(Y_freq, N*M, 1); % Vector quan sát
        
        % 2. Xây dựng Ma trận Hệ thống H_OFDM cho MPA
        % Quan hệ: y_vec = H_OFDM * x_vec + noise
        % H_OFDM = (I_N kron F) * G * (I_N kron F')
        
        Full_F = kron(In, Fn);       % Biến đổi FFT toàn khung
        Full_IF = kron(In, Fn');     % Biến đổi IFFT toàn khung
        H_OFDM = Full_F * G * Full_IF;
        
        % Làm thưa ma trận H_OFDM để MPA chạy nổi (Thresholding)
        % Nhiễu ICI trong OFDM làm ma trận dày, cần lọc bỏ giá trị nhỏ
        H_OFDM(abs(H_OFDM) < 1e-4) = 0;
        
        % MPA Detection OFDM
        % Lưu ý: Hàm MPA_detector được thiết kế tổng quát (y = Hx + n), nên dùng được cho cả OFDM
        [est_bits_ofdm, ~, ~] = MPA_detector(N, M, M_mod, sigma_2, data_grid_logic, y_vec_ofdm, H_OFDM, n_ite_MPA);
        
        total_err_ofdm = total_err_ofdm + sum(xor(est_bits_ofdm, trans_bit));
        
        total_bits = total_bits + length(trans_bit);
    end
    
    % Tính BER
    ber_otfs(iesn0) = total_err_otfs / total_bits;
    ber_ofdm(iesn0) = total_err_ofdm / total_bits;
    
    fprintf('SNR: %2d dB | BER OTFS: %.2e | BER OFDM: %.2e\n', SNR_dB(iesn0), ber_otfs(iesn0), ber_ofdm(iesn0));
end

%% 3. Vẽ đồ thị
figure;
semilogy(SNR_dB, ber_ofdm, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OFDM + MPA');
hold on;
semilogy(SNR_dB, ber_otfs, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OTFS + MPA');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title(['Comparison of MPA Detection: OFDM vs OTFS (Speed: ' num2str(max_speed) ' km/h)']);
legend('show');