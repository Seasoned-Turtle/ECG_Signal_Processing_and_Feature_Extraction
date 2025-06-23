%%
clc; clear; close all;

%%
% 读取心电信号文件（严格按照peakdetection for experiment example.txt中的代码）
fid = fopen('s0379lre.dat', 'r');
a1 = fread(fid, [12, inf], 'int16');
fclose(fid);
b1 = a1';
fid = fopen('s0379lre.xyz', 'r');
a2 = fread(fid, [3, inf], 'int16');
fclose(fid);
b2 = a2';

% 选取前15000个数据点进行分析（根据任务书要求）
ST = 15000;
X11 = 1/2000 * b1(1:ST, 7); % 选取V1导联
X12 = 1/2000 * b1(1:ST, 8); % 选取V2导联
X13 = 1/2000 * b1(1:ST, 10); % 选取V4导联
X14 = 1/2000 * b1(1:ST, 9); % 选取V3导联
X2 = 1/2000 * b2(1:ST, :); % 选取其他导联

% 2.2 心电信号的频谱和噪音分析
fs = 1000; % 采样频率1000Hz
N = length(X11); % 信号长度
t = (0:N-1)/fs; % 时间向量

% 绘制时域信号
figure;
%subplot(2, 2, 1);
plot(t, X11);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal in Time Domain (V1)');

figure;
%subplot(2, 2, 2);
plot(t, X12);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal in Time Domain (V2)');

figure;
%subplot(2, 2, 3);
plot(t, X13);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal in Time Domain (V4)');

figure;
%subplot(2, 2, 4);
plot(t, X14);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal in Time Domain (V3)');
%% 

% 计算频谱
Y11 = fft(X11);
P2_11 = abs(Y11/N);
P1_11 = P2_11(1:N/2+1);
P1_11(2:end-1) = 2*P1_11(2:end-1);
f = fs*(0:(N/2))/N;

% 绘制频谱
figure;
plot(f, P1_11);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of ECG Signal (V1)');
%%
% 2.3 IIR和FIR滤波器技术指标的确定
% 根据频谱分析结果确定滤波器指标
% 需要滤除0 Hz（基线漂移）和50 Hz（工频干扰）
f_cutoff_IIR_low = 0.1; % 高通滤波器截止频率，去除基线漂移
f_cutoff_IIR_high = 40; % 低通滤波器截止频率，去除50 Hz工频干扰
order_IIR = 4; % 滤波器阶数

f_cutoff_FIR_low = 6; % 高通滤波器截止频率
f_cutoff_FIR_high = 40; % 低通滤波器截止频率
order_FIR = 50; % 滤波器阶数

% 2.3 IIR和FIR滤波器设计
% IIR滤波器设计
[b_IIR_low, a_IIR_low] = butter(order_IIR, f_cutoff_IIR_low/(fs/2), 'high'); % 高通滤波器
[b_IIR_high, a_IIR_high] = butter(order_IIR, f_cutoff_IIR_high/(fs/2), 'low'); % 低通滤波器

% FIR滤波器设计
b_FIR_low = fir1(order_FIR, f_cutoff_FIR_low/(fs/2), 'high'); % 高通滤波器
b_FIR_high = fir1(order_FIR, f_cutoff_FIR_high/(fs/2), 'low'); % 低通滤波器
%%
% 2.4 基于IIR和FIR心电信号滤波
% IIR滤波
filtered_IIR = filtfilt(b_IIR_high, a_IIR_high, filtfilt(b_IIR_low, a_IIR_low, X11));

% FIR滤波
filtered_FIR = filtfilt(b_FIR_high, 1, filtfilt(b_FIR_low, 1, X11));

% 绘制滤波后的信号
figure;

%subplot(3, 1, 1);
plot(t, X11);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal in Time Domain (V1)');

figure;
%subplot(3, 1, 2);
plot(t, filtered_IIR);
xlabel('Time (s)');
ylabel('Amplitude');
title('IIR Filtered ECG Signal (V1)');

figure;
%subplot(3, 1, 3);
plot(t, filtered_FIR);
xlabel('Time (s)');
ylabel('Amplitude');
title('FIR Filtered ECG Signal (V1)');
%%
% 2.5 心电信号峰值检测和特征提取
% 峰值检测
% 使用Pan-Tompkins算法检测R波峰值
% 计算微分信号
diff_signal = diff(filtered_IIR);
% 平方信号
squared_signal = diff_signal.^2;
% 移动平均滤波
window_size = 30; % 移动平均窗口大小
integrated_signal = movmean(squared_signal, window_size);
% 阈值检测
threshold = max(integrated_signal) * 0.4; % 设置阈值为信号最大值的40%
peaks = find(integrated_signal > threshold);

% 绘制峰值检测结果
figure;
plot(t(1:end-1), integrated_signal);
hold on;
plot(peaks/fs, integrated_signal(peaks), 'r*');
xlabel('Time (s)');
ylabel('Amplitude');
title('Detected Peaks of IIR Filtered ECG Signal (V1)');

% 特征提取
% 提取信号的能量、功率谱等特征
energy_IIR = sum(filtered_IIR.^2);
energy_FIR = sum(filtered_FIR.^2);

% 计算功率谱
[P_IIR, f_IIR] = periodogram(filtered_IIR, [], [], fs);
[P_FIR, f_FIR] = periodogram(filtered_FIR, [], [], fs);

% 绘制功率谱
figure;
%subplot(2, 1, 1);
plot(f_IIR, 10*log10(P_IIR));
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectrum of IIR Filtered ECG Signal (V1)');

figure
%subplot(2, 1, 2);
plot(f_FIR, 10*log10(P_FIR));
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectrum of FIR Filtered ECG Signal (V1)');

% 输出特征值
fprintf('Energy of IIR Filtered Signal: %f\n', energy_IIR);
fprintf('Energy of FIR Filtered Signal: %f\n', energy_FIR);
%%
% 2.6 滤波结束后的频谱分析
% 计算滤波后的频谱
Y_IIR = fft(filtered_IIR);
P2_IIR = abs(Y_IIR/N);
P1_IIR = P2_IIR(1:N/2+1);
P1_IIR(2:end-1) = 2*P1_IIR(2:end-1);

Y_FIR = fft(filtered_FIR);
P2_FIR = abs(Y_FIR/N);
P1_FIR = P2_FIR(1:N/2+1);
P1_FIR(2:end-1) = 2*P1_FIR(2:end-1);
%%
% 绘制滤波后的频谱
figure;
%subplot(3, 1, 1);
plot(f, P1_11);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of ECG Signal (V1)');

figure;
%subplot(3, 1, 2);
plot(f, P1_IIR);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of IIR Filtered ECG Signal (V1)');

figure;
%subplot(3, 1, 3);
plot(f, P1_FIR);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of FIR Filtered ECG Signal (V1)');