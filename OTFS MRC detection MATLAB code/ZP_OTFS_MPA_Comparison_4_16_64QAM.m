% =========================================================================
% File: ZP_OTFS_MPA_Comparison_4_16_64QAM.m
% Dựa trên: ZP_OTFS_MRC_system.m
% Mục đích: So sánh hiệu năng MPA Detector với 4-QAM, 16-QAM, 64-QAM
% =========================================================================

close all
clear all
rng(1)

%% 1. Thiết lập tham số cơ bản
N = 16;             % Số ký hiệu thời gian
M = 64;             % Số sóng mang con
length_ZP = M/16;   % Zero Padding
M_data = M - length_ZP;

% Tạo lưới dữ liệu
data_grid = zeros(M,N);
data_grid(1:M_data, 1:N) = 1;
N_syms_perfram = sum(sum(data_grid));

% Thông số vật lý
car_fre = 4*10^9;
delta_f = 15*10^3;
T = 1/delta_f;
max_speed = 300;    % Tốc độ cao (km/h) để thử thách MPA

% Cấu hình mô phỏng
SNR_dB = 0:5:30;    % Dải SNR rộng để bao quát cả 64-QAM
SNR = 10.^(SNR_dB/10);

% Danh sách các bậc điều chế cần so sánh
M_mod_list = [4, 16, 64];
colors = {'b-o', 'r-s', 'k-^'}; % Màu vẽ cho 3 đường

% Biến lưu kết quả BER
ber_results = zeros(length(M_mod_list), length(SNR_dB));

% Số frame (Giảm xuống để chạy nhanh, tăng lên 1000 khi cần kết quả đẹp)
N_fram = 50; 

%% 2. Vòng lặp chính qua các bậc điều chế (4 -> 16 -> 64)
fprintf('Bắt đầu so sánh MPA với 4-QAM, 16-QAM, và 64-QAM...\n');

for idx_mod = 1:length(M_mod_list)
    M_mod = M_mod_list(idx_mod);
    M_bits = log2(M_mod);
    N_bits_perfram = N_syms_perfram * M_bits;
    
    % Chuẩn hóa năng lượng (Quan trọng khi đổi QAM)
    if M_mod == 2
        eng_sqrt = 1;
    else
        eng_sqrt = sqrt((M_mod-1)/6*(2^2));
    end
    sigma_2 = (abs(eng_sqrt)^2)./SNR;
    
    fprintf('--- Đang chạy %d-QAM ---\n', M_mod);
    
    Fn = dftmtx(N);
    Fn = Fn./norm(Fn);
    
    % Vòng lặp SNR
    for iesn0 = 1:length(SNR_dB)
        total_errors = 0;
        
        for ifram = 1:N_fram
            % A. Tạo dữ liệu
            trans_info_bit = randi([0,1], N_syms_perfram*M_bits, 1);
            data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), M_mod, 'gray', 'InputType', 'bit');
            X = Generate_2D_data_grid(N, M, data, data_grid);
            
            % OTFS Modulation
            X_tilda = X * Fn';
            s = reshape(X_tilda, N*M, 1);
            
            % B. Kênh truyền (Delay-Doppler)
            [chan_coef, delay_taps, Doppler_taps, taps] = Generate_delay_Doppler_channel_parameters(N, M, car_fre, delta_f, T, max_speed);
            [G, gs] = Gen_time_domain_channel(N, M, taps, delay_taps, Doppler_taps, chan_coef);
            
            % Tạo ma trận kênh H miền DD cho MPA (QUAN TRỌNG)
            [H, ~, ~] = Gen_DD_and_DT_channel_matrices(N, M, G, Fn);
            
            % Tín hiệu thu
            noise = sqrt(sigma_2(iesn0)/2) * (randn(size(s)) + 1i*randn(size(s)));
            r = G*s + noise;
            
            % Demodulation (về miền Delay-Doppler cho MPA)
            Y_tilda = reshape(r, M, N);
            Y = Y_tilda * Fn;
            y = reshape(Y.', N*M, 1); % Vector quan sát y
            
            % C. MPA Detection
            n_ite_MPA = 10; % Số vòng lặp (Với 64-QAM có thể cần tăng lên 15-20)
            
            % Gọi hàm MPA
            [est_info_bits_MPA, ~, ~] = MPA_detector(N, M, M_mod, sigma_2(iesn0), data_grid, y, H, n_ite_MPA);
            
            % D. Đếm lỗi
            errors = sum(xor(est_info_bits_MPA, trans_info_bit));
            total_errors = total_errors + errors;
        end
        
        % Tính BER
        ber_results(idx_mod, iesn0) = total_errors / (N_fram * N_bits_perfram);
        fprintf('   SNR: %2d dB | BER: %.5f\n', SNR_dB(iesn0), ber_results(idx_mod, iesn0));
    end
end

%% 3. Vẽ đồ thị so sánh
figure;
for idx_mod = 1:length(M_mod_list)
    semilogy(SNR_dB, ber_results(idx_mod, :), colors{idx_mod}, 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', sprintf('MPA - %d-QAM', M_mod_list(idx_mod)));
    hold on;
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Performance Comparison of MPA Detector (4/16/64-QAM)');
legend('show');
axis([min(SNR_dB) max(SNR_dB) 1e-5 1]);