% =========================================================================
% Script: Compare_MPA_QAM_orders.m
% Mục đích: So sánh hiệu năng thuật toán MPA với 4-QAM, 16-QAM, 64-QAM
% =========================================================================

close all
clear all
rng(1) % Đặt seed để tái lập kết quả

%% 1. Thiết lập tham số chung
N = 16;            % Số ký hiệu thời gian
M = 64;            % Số sóng mang con tần số
length_ZP = M/16;  % Độ dài Zero Padding
M_data = M - length_ZP;
data_grid = zeros(M,N);
data_grid(1:M_data,1:N) = 1; % Lưới dữ liệu hoạt động

% Thông số vật lý
car_fre = 4*10^9;  % Tần số sóng mang 4GHz
delta_f = 15*10^3; % Khoảng cách sóng mang 15kHz
T = 1/delta_f;
max_speed = 500;   % Tốc độ xe (km/h) -> Doppler cao

% Cấu hình chạy mô phỏng
M_mod_list = [4, 16, 64]; % Các chế độ so sánh
SNR_dB = 0:5:35;          % Dải SNR (dB)
N_fram = 100;             % Số frame (Giảm xuống để test nhanh 64-QAM)
n_ite_MPA = 10;           % Số vòng lặp MPA

% Ma trận lưu kết quả BER
ber_results = zeros(length(M_mod_list), length(SNR_dB));

%% 2. Vòng lặp chính qua từng chế độ QAM
for idx_mod = 1:length(M_mod_list)
    M_mod = M_mod_list(idx_mod);
    M_bits = log2(M_mod);
    
    % Chuẩn hóa năng lượng (Quan trọng khi đổi QAM)
    if M_mod == 2
        eng_sqrt = 1;
    else
        eng_sqrt = sqrt((M_mod-1)/6*(2^2));
    end
    
    % Tính toán phương sai nhiễu
    SNR = 10.^(SNR_dB/10);
    sigma_2 = (abs(eng_sqrt)^2)./SNR;
    
    fprintf('Đang chạy mô phỏng cho %d-QAM...\n', M_mod);
    
    % Chuẩn bị ma trận DFT
    Fn = dftmtx(N);
    Fn = Fn./norm(Fn);
    
    % Vòng lặp SNR
    for iesn0 = 1:length(SNR_dB)
        total_errors = 0;
        total_bits = 0;
        
        % Vòng lặp Frame
        for ifram = 1:N_fram
            % -------------------------------------------------------------
            % A. Tạo dữ liệu và điều chế
            % -------------------------------------------------------------
            N_syms_perfram = sum(sum(data_grid));
            trans_info_bit = randi([0,1], N_syms_perfram*M_bits, 1);
            
            data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), M_mod, 'gray', 'InputType', 'bit');
            X = Generate_2D_data_grid(N, M, data, data_grid);
            
            % OTFS Modulation
            X_tilda = X * Fn';
            s = reshape(X_tilda, N*M, 1);
            
            % -------------------------------------------------------------
            % B. Kênh truyền (Time-Varying Channel)
            % -------------------------------------------------------------
            [chan_coef, delay_taps, Doppler_taps, taps] = Generate_delay_Doppler_channel_parameters(N, M, car_fre, delta_f, T, max_speed);
            [G, gs] = Gen_time_domain_channel(N, M, taps, delay_taps, Doppler_taps, chan_coef);
            [H, H_tilda, P] = Gen_DD_and_DT_channel_matrices(N, M, G, Fn);
            
            % Tín hiệu thu + Nhiễu
            noise = sqrt(sigma_2(iesn0)/2) * (randn(size(s)) + 1i*randn(size(s)));
            r_vec = G * s + noise; % Tính r dùng ma trận G cho chính xác
            
            % Demodulation chuẩn bị cho MPA
            Y_tilda = reshape(r_vec, M, N);
            Y = Y_tilda * Fn;
            y = reshape(Y.', N*M, 1); % Vector tín hiệu miền Delay-Doppler
            
            % -------------------------------------------------------------
            % C. Thuật toán MPA
            % -------------------------------------------------------------
            % Gọi hàm MPA_detector có sẵn
            [est_info_bits_MPA, ~, ~] = MPA_detector(N, M, M_mod, sigma_2(iesn0), data_grid, y, H, n_ite_MPA);
            
            % -------------------------------------------------------------
            % D. Đếm lỗi
            % -------------------------------------------------------------
            errors = sum(xor(est_info_bits_MPA, trans_info_bit));
            total_errors = total_errors + errors;
            total_bits = total_bits + length(trans_info_bit);
        end
        
        % Tính BER trung bình cho mức SNR này
        avg_ber = total_errors / total_bits;
        ber_results(idx_mod, iesn0) = avg_ber;
        
        fprintf('   SNR: %2d dB | BER: %.2e\n', SNR_dB(iesn0), avg_ber);
    end
    fprintf('Hoàn thành %d-QAM.\n\n', M_mod);
end

%% 3. Vẽ đồ thị so sánh
figure('Name', 'MPA Performance Comparison');
semilogy(SNR_dB, ber_results(1,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'MPA - 4QAM');
hold on;
semilogy(SNR_dB, ber_results(2,:), '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'MPA - 16QAM');
semilogy(SNR_dB, ber_results(3,:), '-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'MPA - 64QAM');

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Performance of MPA Detector with Different Modulation Orders');
legend('show', 'Location', 'southwest');
axis([min(SNR_dB) max(SNR_dB) 1e-6 1]);