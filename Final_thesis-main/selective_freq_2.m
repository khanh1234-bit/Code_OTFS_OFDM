%BER of OFDM in Frequency selective channel
clc;
clear all;
close all;

% --- Cấu hình tham số ---
N = 128;                             % Số lượng sóng mang con (Subcarriers)
Ncp = 16;                            % Độ dài Cyclic prefix
Ts = 1e-3;                           % Chu kỳ lấy mẫu
Fd = 0;                              % Độ dịch Doppler tối đa (0 = kênh tĩnh trong 1 frame)
Np = 4;                              % Số lượng pilot
M = 2;                               % Số mức điều chế (BPSK)
Nframes = 10^3;                      % Số lượng khung OFDM mô phỏng

% --- Tạo dữ liệu phát ---
D = round((M-1)*rand((N-2*Np),Nframes));
const = pskmod([0:M-1],M);
Dmod = pskmod(D,M);
Data = [zeros(Np,Nframes); Dmod ; zeros(Np,Nframes)];   % Chèn Pilot

% --- OFDM Transmitter ---
IFFT_Data = (N/sqrt(120))*ifft(Data,N);
TxCy = [IFFT_Data((N-Ncp+1):N,:); IFFT_Data];  % Thêm Cyclic prefix
[r, c] = size(TxCy); % r = độ dài 1 symbol OFDM (có CP), c = số lượng frame
Tx_Data = TxCy;

% --- Cấu hình kênh Rayleigh (Sử dụng comm.RayleighChannel thay cho rayleighchan) ---
tau = [0 1e-5 3.5e-5 12e-5];         % Độ trễ các đường (Path delays)
pdb = [0 -1 -1 -3];                  % Công suất trung bình các đường (dB)

% Khởi tạo đối tượng kênh
chan = comm.RayleighChannel(...
    'SampleRate', 1/Ts, ...
    'PathDelays', tau, ...
    'AveragePathGains', pdb, ...
    'MaximumDopplerShift', Fd, ...
    'PathGainsOutputPort', true);    % Cho phép xuất hệ số kênh để dùng cho Equalizer

% --- Chuẩn bị ma trận tính toán đáp ứng tần số (Cho Equalizer) ---
% Công thức: H(k) = sum( path_gain * exp(-j*2*pi*k*tau/T_symbol) )
k_indices = 0:N-1; 
% Lưu ý: Chu kỳ của symbol OFDM (không kể CP) là T_fft = N * Ts
% Tuy nhiên trong công thức chuẩn hóa FFT, ta thường dùng tỷ lệ tau/Ts
fft_exp_matrix = exp(-1i * 2 * pi * k_indices' * (tau / (N*Ts))); 

% --- SNR Configuration ---
EbNo = 0:5:30;
EsNo= max(0, EbNo + 10*log10(120/128)+ 10*log10(128/144));
snr= EsNo - 10*log10(128/144);

% --- Vòng lặp mô phỏng ---
berofdm = zeros(1,length(snr));
Rx_Data = zeros((N-2*Np),Nframes);

fprintf('Đang chạy mô phỏng...\n');

for i = 1:length(snr)
    for j = 1:c % Duyệt qua từng frame
        
        % 1. Lọc qua kênh Fading
        % Reset kênh mỗi frame để giả lập Block Fading độc lập (tương tự ResetBeforeFiltering cũ)
        reset(chan); 
        
        tx_signal = Tx_Data(:,j); % Lấy 1 frame (vector cột)
        [hx, path_gains] = chan(tx_signal); 
        
        % 2. Tính toán đáp ứng tần số kênh (Channel Frequency Response) cho frame này
        % Lấy mẫu path_gains tại thời điểm đầu tiên (do Fd=0 nên kênh không đổi trong frame)
        g_snapshot = path_gains(1, :); 
        
        % Tính G(j,:) = DFT của đáp ứng xung kim kênh
        H_freq = (g_snapshot * fft_exp_matrix.').'; 
        G_current = H_freq; 
        
        % 3. Thêm nhiễu AWGN
        y = awgn(hx, snr(i), 'measured');
        
        % --- Receiver ---
        % Loại bỏ Cyclic Prefix
        Rx = y(Ncp+1:r);                                
        
        % Biến đổi FFT và Cân bằng kênh (FDE - Frequency Domain Equalization)
        FFT_Data = (sqrt(120)/128) * fft(Rx,N) ./ G_current;
        
        % Loại bỏ Pilot và Giải điều chế
        Rx_Data(:,j) = pskdemod(FFT_Data(5:124), M);     
    end
    
    % Tính BER cho mức SNR hiện tại
    berofdm(i) = sum(sum(Rx_Data~=D))/((N-2*Np)*Nframes);
    fprintf('SNR: %d dB completed. BER: %e\n', EbNo(i), berofdm(i));
end

% --- Vẽ đồ thị ---
figure;
semilogy(EbNo, berofdm, '--or', 'linewidth', 2);
hold on;
grid on;
title('OFDM BER vs SNR in Frequency selective Rayleigh fading channel');
xlabel('SNR (dB)'); % Thực tế đây là Eb/No theo trục vẽ
ylabel('BER');
legend('Simulation');

% Tùy chọn: Vẽ đường lý thuyết (nếu cần)
% ber_theo = berawgn(EbNo, 'psk', M, 'nondiff');
% semilogy(EbNo, ber_theo, 'b-');