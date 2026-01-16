clc; clear; close all;

% 1. Tạo trục tần số (độ phân giải cao để đường cong mượt)
f = -10:0.001:10;

% 2. Tạo 3 hàm Sinc (Lưu ý: dùng abs() để lấy giá trị tuyệt đối như hình)
% Sinc 1 (Đỏ): Dịch về -1
y1 = abs(sinc(f + 1));
% Sinc 2 (Xanh lá): Ở giữa (0)
y2 = abs(sinc(f));
% Sinc 3 (Xanh dương): Dịch về 1
y3 = abs(sinc(f - 1));

% 3. Vẽ đồ thị
figure('Color', 'w'); % Nền trắng
hold on;
plot(f, y1, 'r', 'LineWidth', 2); % Vẽ đường Đỏ
plot(f, y2, 'g', 'LineWidth', 2); % Vẽ đường Xanh lá
plot(f, y3, 'b', 'LineWidth', 2); % Vẽ đường Xanh dương

% 4. Trang trí giống hệt hình mẫu
grid on;
axis([-10 10 0 1.2]); % Giới hạn trục X và Y
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
title('Fourier Transforms of Rectangular Functions (Sinc Functions)', 'FontWeight', 'bold');
legend('Sinc 1', 'Sinc 2', 'Sinc 3', 'Location', 'northeast');
hold off;