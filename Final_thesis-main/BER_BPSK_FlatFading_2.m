clc;
clear all;
close all;

% --- BẮT ĐẦU ĐẾM THỜI GIAN ---
tic;
fprintf('--------------------------------------------------\n');
fprintf('Bắt đầu mô phỏng BER BPSK trong Flat Fading...\n');
fprintf('--------------------------------------------------\n');

% --- Tham số mô phỏng ---
N = 10^6;               % Số lượng bit
snr_dB = -10:2:30;      % Dải SNR (dB). Giảm bước nhảy để chạy nhanh hơn (từ 0.1 lên 2)
% Nếu bạn muốn mịn hơn có thể để 0.5 hoặc 1
snr_lin = 10.^(snr_dB/10);
len = length(snr_dB);

% --- Tạo dữ liệu phát ---
x = randi([0,1], 1, N); 
x_in = 2*x - 1;         % Điều chế BPSK: 0 -> -1, 1 -> 1

BER_flat_fading_BPSK = zeros(1, len);

% --- Vòng lặp SNR ---
for i = 1:len
    fprintf('Đang chạy SNR = %2d dB... ', snr_dB(i));
    
    % 1. Tạo kênh Rayleigh Flat Fading
    % Chuẩn hóa công suất kênh về 1 (chia cho sqrt(2))
    h = (randn(1,N) + 1i*randn(1,N)) / sqrt(2); 
    
    % 2. Truyền qua kênh
    c = h .* x_in;
    
    % 3. Thêm nhiễu AWGN
    % Lưu ý: awgn 'measured' sẽ đo công suất tín hiệu c để thêm nhiễu
    y = awgn(c, snr_dB(i), 'measured');
    
    % 4. Bộ thu (Equalizer - Zero Forcing)
    % Chia tín hiệu nhận được cho hệ số kênh để loại bỏ pha và biên độ kênh
    y_rec = y ./ h;
    
    % 5. Quyết định cứng (Hard Decision) - Vector hóa (Không dùng vòng lặp)
    % Nếu phần thực > 0 -> bit 1, ngược lại -> bit -1
    % Dùng hàm sign: sign(x) trả về 1 nếu x>0, -1 nếu x<0
    y_dec = sign(real(y_rec)); 
    % Xử lý trường hợp = 0 (hiếm gặp nhưng cần thiết) gán về 1 hoặc -1 tùy ý
    y_dec(y_dec == 0) = 1; 
    
    % 6. Đếm lỗi
    err_count = sum(y_dec ~= x_in);
    
    % Tính BER
    BER_flat_fading_BPSK(i) = err_count / N;
    
    fprintf('Xong. BER = %.4e\n', BER_flat_fading_BPSK(i));
end

% --- Tính đường lý thuyết để so sánh ---
% Lý thuyết AWGN: 0.5 * erfc(sqrt(SNR))
BER_theoretical_AWGN = 0.5 * erfc(sqrt(snr_lin));

% Lý thuyết Rayleigh: 0.5 * (1 - sqrt(SNR / (1 + SNR)))
BER_theoretical_Rayleigh = 0.5 * (1 - sqrt(snr_lin ./ (1 + snr_lin)));

% --- Kết thúc đếm thời gian ---
total_time = toc;
fprintf('--------------------------------------------------\n');
fprintf('Mô phỏng hoàn tất trong %.4f giây.\n', total_time);
fprintf('--------------------------------------------------\n');

% --- Vẽ đồ thị ---
figure;

% 1. Đường AWGN Lý thuyết (Màu xanh)
semilogy(snr_dB, BER_theoretical_AWGN, 'b-', 'LineWidth', 2);hold on;

% 2. Đường Fading Lý thuyết (Màu đen nét đứt)
semilogy(snr_dB, BER_theoretical_Rayleigh, 'k--', 'LineWidth', 2);hold on;

% 3. Đường Mô phỏng Fading (Màu đỏ tròn)
semilogy(snr_dB, BER_flat_fading_BPSK, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);hold on;

grid on;
legend('Theoretical AWGN', 'Theoretical Rayleigh', 'Simulated Flat Fading');
title('BER of BPSK in Flat Fading vs AWGN');
ylabel("Bit Error Rate (BER)");
xlabel("SNR (dB)");
ylim([1e-5 1])