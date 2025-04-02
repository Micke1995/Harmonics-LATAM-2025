clear all
clc
% load caso_2
% load('senales_kundur_18_may.mat')
load caso_Zamora_01 %load caso_19may
t = simout.Time;
signals = simout.Data;
F1 = 60;
Ts = t(2) - t(1);
Fs = 1/Ts; % Sample frequency
No = round(Fs/F1); %samples per period
% Ts = 1/Fs;
figure;
Vabc = signals(end-12*No:end,19:21);
plot(t(end-12*No:end), Vabc)
xlabel('Time (s)')
ylabel('Voltage (pu)')
legend('v_a', 'v_b', 'v_c', 'Orientation', 'Horizontal')
grid on
set(gca, 'LooseInset', [0,0,0,0]);
%
Iabc = signals(end-12*No:end,19);
N = length(Iabc)-1;
fun = Iabc(1:N)';
fun2 = Iabc(1:N+1)';
k = 12;
r = round(N/2) - 1;
H0 = hankel(fun(1:r),fun(r:N-2));
H1 = hankel(fun(2:r+1),fun(r+1:N-1));
[U,S,V] = svds(H0,k);
A = (S^-(1/2))*U'*H1*V*(S^-(1/2));
z = eig(A);
pot = 0:N;%-N/2:N/2;%
for m = 1:length(z)
    ZZ(:,m) = ( z(m) ).^pot; % normal
end
B = pinv(ZZ)*fun2';
landa = log(z)/Ts;
sigma = real(landa);
omega = imag(landa);
Frec = (omega/(2*pi));
damp=sigma;
damp_ratio = 100*damp ./ omega;
Amp = 2*abs(B);
theta = angle(B);
rows = find(Frec > 0);
ampli = Amp(rows);
freq = Frec(rows);
phase = theta(rows);
dampr = damp_ratio(rows);
rows2 = find(dampr>0);
ampli2 = ampli;%(rows2);
freq2 = freq;%(rows2);
phase2 = phase;%(rows2);
dampr2 = dampr;%(rows2);
todo1 = [freq2 ampli2 phase2 dampr2];
todo_orden = sortrows(todo1);


 disp('------------------ERA------------------')
disp('    Freq       Amp       Phase    DampR')
todo2 = todo_orden;
disp(sortrows(todo2))


% IHD(:,k) = 100*(amp_h(2:end,k))./(amp_h(1,k));
amp_era = [todo2(6,2) todo2(6,2) todo2(5,2) todo2(5,2) todo2(2,2) todo2(2,2) todo2(1,2) todo2(1,2) todo2(4,2) todo2(4,2) todo2(3,2) todo2(3,2)];
suma = sum(todo2(3:end,2).^2) - todo2(4,2).^2;
fundamental = max(todo2(:,2));
THD_ERA = 100*sqrt(suma) ./ fundamental
harmonics_ERA = zeros(10,1);
harmonics_ERA([1 5 7 11]) = 100*[todo2(2:3,2) ; todo2(5:6,2)]./max(todo2(:,2));

%% FFT
load FFT_Kundur_harmonics
specT = spec((N+1)/2+12:12:end);
harmonics_FFT = zeros(10,1);
harmonics_FFT([1 5 7 11]) = 100*[specT(1) specT(5) specT(7) specT(11)]./max(spec2);
suma = sum(spec2(9:end).^2) - spec2(10).^2;
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental

figure;
harmonics_Va = [harmonics_FFT harmonics_ERA];
bar(harmonics_Va)
ylabel('% of fundamental')
xlabel('Harmonic number, #')
% ylim([0 5])
xlim([0 9.5])
legend('FFT', 'ERA')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

I_era_Va = 0;
I_fft_Va = 0;
for i = 1:6
    I_era_Va = todo2(i,2).*cos(todo2(i,1)*2*pi*t(1:N) + todo2(i,3) ) + I_era_Va;
    I_fft_Va = spec2(i+6).*cos( fre_fft(i+6)*2*pi*t(1:N) + aspec2(i+6) ) + I_fft_Va;
end

figure;
Va = fun';
tv = t(1:N);
plot(t(1:N), Va, t(1:N), I_fft_Va, '-.', t(1:N), I_era_Va, '--')
ylabel('Voltage (pu)')
xlabel('Time (s)')
legend('Actual', 'FFT', 'ERA', 'Orientation', 'Horizontal')
xlim([0 0.2])
grid on
set(gca, 'LooseInset', [0,0,0,0]);

RMSE_FFT = sqrt(mean((fun' - I_fft_Va).^2)) 
RMSE_ERA = sqrt(mean((fun' - I_era_Va).^2))
a = norm(fun',2)^2;
b_fft = norm(fun,2)^2 - norm(I_fft_Va,2)^2;
b_era = norm(fun,2)^2 - norm(I_era_Va,2)^2;
SNR_FFT = 10*log10(a/b_fft);
SNR_ERA = 10*log10(a/b_era);
SSE_FFT = norm((fun' - I_fft_Va),2)^2;
% SSE_FFT = sum((fun' - I_fft).^2);
SSE_ERA = sum((fun' - I_era_Va).^2);
MSE_FFT = (1/N)*SSE_FFT;
MSE_ERA = (1/N)*SSE_ERA;
RMSE_fft = sqrt(MSE_FFT);
RMSE_era = sqrt(MSE_ERA);
Errors_FFT = [SNR_FFT ; SSE_FFT ; MSE_FFT ; RMSE_FFT]
Errors_ERA = [SNR_ERA ; SSE_ERA ; MSE_ERA ; RMSE_ERA]

% % % %% Frequency Response
tam = N;%2^15;
fn = Fs*(-tam/2:tam/2)/(tam*F1);%/(tam*F1);%linspace(-Fs/2,Fs/2,tam);
figure;
subaxis(2, 1, 1, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.11, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.2)
Zp = pinv(ZZ);
for kk = 1:k
    ir = Zp(kk,:);
    IR(kk,:) = fftshift(fft(ir));
    plot( fn, amp_era(kk)*abs( IR(kk,:) ));%/max(abs(max(IR(kk,:)))) ) );
    hold on
end
xlim([-12 12])
ylim([0 0.5])
set(gca,'XTickLabel',[])
ylabel('Magnitude')
% xlabel('Normalised Frequency u=f_1T_1')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

subaxis(2, 1, 2, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.11, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.2)
plot(fn, spec);
xlim([-12 12])
ylim([0 0.5])
% ylabel('Magnitude')
xlabel('Normalised Frequency u=f_1T_1')
grid on




