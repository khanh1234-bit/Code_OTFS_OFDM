close all
clear all
rng(1)

%% --- 1. OTFS Parameters ---
N = 16;             % Số ký hiệu thời gian
M = 64;             % Số sóng mang con
M_mod = 4;         % 64-QAM (Lưu ý: MPA chạy rất chậm với 64-QAM)
M_bits = log2(M_mod);

% Năng lượng trung bình
if M_mod == 2
    eng_sqrt = 1;
else
    eng_sqrt = sqrt((M_mod-1)/6*(2^2));
end

% Delay-Doppler Grid
length_ZP = M/16;
M_data = M - length_ZP;
data_grid = zeros(M,N);
data_grid(1:M_data, 1:N) = 1;
N_syms_perfram = sum(sum(data_grid));
N_bits_perfram = N_syms_perfram * M_bits;

% Thông số vật lý
car_fre = 4*10^9;
delta_f = 15*10^3;
T = 1/delta_f;

% SNR Setup
SNR_dB = 5:5:30; % Dải SNR
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;

%% --- 2. Simulation Setup ---
% CẢNH BÁO: Với 64-QAM, MPA chạy rất lâu. 
% Tôi giảm N_fram xuống 10 để bạn test code nhanh. 
% Hãy tăng lên 100-1000 khi chạy thật để có đồ thị mượt.
N_fram = 10; 
fprintf('Đang chạy mô phỏng MPA với %d-QAM và %d Frames...\n', M_mod, N_fram);

% Biến lưu kết quả
avg_ber_MPA = zeros(1, length(SNR_dB));
Fn = dftmtx(N);      
Fn = Fn./norm(Fn);   

%% --- 3. Main Loop ---
for iesn0 = 1:length(SNR_dB)
    total_errors_MPA = 0;
    
    for ifram = 1:N_fram
        % 3.1 Tạo dữ liệu
        trans_info_bit = randi([0,1], N_syms_perfram*M_bits, 1);
        data = qammod(reshape(trans_info_bit, M_bits, N_syms_perfram), M_mod, 'gray', 'InputType', 'bit');
        X = Generate_2D_data_grid(N, M, data, data_grid);
        
        % 3.2 OTFS Modulation
        X_tilda = X * Fn';
        s = reshape(X_tilda, N*M, 1);
        
        % 3.3 Channel Generation (Delay-Doppler)
        max_speed = 500; % km/h
        [chan_coef, delay_taps, Doppler_taps, taps] = Generate_delay_Doppler_channel_parameters(N, M, car_fre, delta_f, T, max_speed);
        [G, gs] = Gen_time_domain_channel(N, M, taps, delay_taps, Doppler_taps, chan_coef);
        
        % Quan trọng: MPA cần ma trận kênh H miền Delay-Doppler
        [H, ~, ~] = Gen_DD_and_DT_channel_matrices(N, M, G, Fn);
        
        % 3.4 Received Signal
        noise = sqrt(sigma_2(iesn0)/2) * (randn(size(s)) + 1i*randn(size(s)));
        
        % Tạo tín hiệu thu r (Time domain) chính xác bằng ma trận G
        r = G * s + noise;
        
        % Demodulation về miền Delay-Doppler (biến y) cho MPA
        Y_tilda = reshape(r, M, N);
        Y = Y_tilda * Fn;
        y = reshape(Y.', N*M, 1);
        
        % 3.5 MPA Detection
        n_ite_MPA = 15; % Số vòng lặp MPA
        
        % Gọi hàm MPA Detector
        % Hàm này phải có sẵn trong thư mục của bạn
        [est_info_bits_MPA, ~, ~] = MPA_detector(N, M, M_mod, sigma_2(iesn0), data_grid, y, H, n_ite_MPA);
        
        % 3.6 Đếm lỗi
        errors = sum(xor(est_info_bits_MPA, trans_info_bit));
        total_errors_MPA = total_errors_MPA + errors;
    end
    
    % Tính BER trung bình tại mức SNR này
    avg_ber_MPA(iesn0) = total_errors_MPA / (N_fram * N_bits_perfram);
    fprintf('SNR: %2d dB | BER: %.5f\n', SNR_dB(iesn0), avg_ber_MPA(iesn0));
end

%% --- 4. Plotting Results ---
figure(1)
semilogy(SNR_dB, avg_ber_MPA, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b');
grid on
xlabel('SNR (dB)')
ylabel('Bit Error Rate (BER)')
title(['MPA Performance (' num2str(M_mod) '-QAM)'])
legend('MPA Detector')