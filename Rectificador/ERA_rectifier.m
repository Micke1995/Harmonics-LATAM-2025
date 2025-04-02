% function [Frec, Amp,theta,thd] = ERA(N,k,Ts,fun)
clear all
clc
load 'Irectifier'
F1 = 60;
Ts = 1e-4;
Fs = 1/Ts; % Sample frequency
No = round(Fs/F1); %samples per period
% Ts = 1/Fs;
figure;
plot(t(1:12*No), Irs(1:12*No,:))
xlabel('Time (s)')
ylabel('Current (A)')
legend('i_a', 'i_b', 'i_c')
grid on
set(gca, 'LooseInset', [0,0,0,0]);
%
Iabc = Irs(1:12*No+1,1);
N = length(Iabc)-1;
fun = Iabc(1:N)';
fun2 = Iabc(1:N+1)';
k = 20;
r = round(N/2) - 1;
H0 = hankel(fun(1:r),fun(r:N-2));
H1 = hankel(fun(2:r+1),fun(r+1:N-1));
[U,S,V] = svds(H0,k);
A = (S^-(1/2))*U'*H1*V*(S^-(1/2));
z = eig(A);
pot = 0:N;%-N/2:N/2;
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
disp(todo2)

algo = [1 5 7 11 13 17 19 23 25 29];

% IHD(:,k) = 100*(amp_h(2:end,k))./(amp_h(1,k));
harmonics_ERA(algo) = 100*todo2(:,2)./max(todo2(:,2));
suma = sum(todo2(2:end,2).^2);
fundamental = max(todo2(:,2));
THD_ERA = 100*sqrt(suma) ./ fundamental

%% FFT
load harmonics_FFT_rectifier
aspec = angle(Esp_sf);
ph_fft = aspec(N/2+13:12:end);
fre_fft = f_sf(N/2+13:12:end)';
amp_fft = spec2;
harmonics_FFT(algo) = 100*spec2(algo)./max(spec2(algo));
% harmonics2_FFT = spec2(:,2);
suma = sum( spec2(algo(2:end)).^2 );
fundamental = max(spec2(algo));
THD_FFT = 100*sqrt(suma) ./ fundamental
% harmonics_FFT(algo2) = 0;
figure;
% bar(harmonics(:,1))
harmonics = [harmonics_FFT' harmonics_ERA'];
bar(harmonics)
% hold on
% bar(harmonics_ERA, 0.5)
ylabel('% of fundamental')
xlabel('Harmonic number, #')
% ylim([0 5])
xlim([-0.5 30])
legend('FFT', 'ERA')
grid on
set(gca, 'LooseInset', [0,0,0,0]);


I_era = 0;
I_fft = 0;
for i = 1:10
    I_era = todo2(i,2).*cos(todo2(i,1)*2*pi*t(1:N) + todo2(i,3) ) + I_era;
    I_fft = amp_fft(algo(i)).*cos(fre_fft(algo(i))*2*pi*t(1:N) + ph_fft(algo(i)) ) + I_fft;
end

figure;
plot(t(1:N), fun', t(1:N), I_fft, '-.', t(1:N), I_era, '--')
ylabel('Current (A)')
xlabel('Time (s)')
legend('Actual', 'FFT', 'ERA')
xlim([0 0.2])
grid on
set(gca, 'LooseInset', [0,0,0,0]);

RMSE_FFT = sqrt(mean((fun' - I_fft).^2)) 
RMSE_ERA = sqrt(mean((fun' - I_era).^2))

a = norm(fun',2)^2;
b_fft = norm(fun,2)^2 - norm(I_fft,2)^2;
b_era = norm(fun,2)^2 - norm(I_era,2)^2;
SNR_FFT = 10*log10(a/b_fft);
SNR_ERA = 10*log10(a/b_era);
SSE_FFT = norm((fun' - I_fft),2)^2;
% SSE_FFT = sum((fun' - I_fft).^2);
SSE_ERA = sum((fun' - I_era).^2);
MSE_FFT = (1/N)*SSE_FFT;
MSE_ERA = (1/N)*SSE_ERA;
RMSE_fft = sqrt(MSE_FFT);
RMSE_era = sqrt(MSE_ERA);
Errors_FFT = [SNR_FFT ; SSE_FFT ; MSE_FFT ; RMSE_FFT]
Errors_ERA = [SNR_ERA ; SSE_ERA ; MSE_ERA ; RMSE_ERA]


% % % %% Frequency Response
tam = N;%2^15;
fn = Fs*(-tam/2:tam/2)/(tam*F1);%linspace(-Fs/2,Fs/2,tam);
amp_era = [todo2(10,2) todo2(10,2) todo2(9,2) todo2(9,2) todo2(8,2) todo2(8,2) todo2(7,2) todo2(7,2) todo2(6,2) todo2(6,2) todo2(1,2) todo2(1,2) todo2(5,2) todo2(5,2) todo2(4,2) todo2(4,2) todo2(2,2) todo2(2,2) todo2(3,2) todo2(3,2)];
figure;
Zp = pinv(ZZ);
subaxis(2, 1, 1, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.13, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.2)
for kk = 1:20
    ir = Zp(kk,:);
    IR(kk,:) = fftshift(fft(ir));
    plot( fn, amp_era(kk)*abs( IR(kk,:)/max(abs(max(IR(kk,:)))) ) );
    hold on
end
xlim([-30 30])
set(gca,'XTickLabel',[])
% title('Frequency responses of filters')
ylabel('Magnitude')
% xlabel('Normalised Frequency u=f_1T_1')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

subaxis(2, 1, 2, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.13, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.2)
plot(f_sf/F1, spec);
xlim([-30 30])
% title('Frequency responses of filters')
% ylabel('Magnitude')
xlabel('Normalised Frequency u=f_1T_1')
grid on


% % % fnf = Fs*(-N/2:N/2-1)/(N*F1);%linspace(-Fs/2,Fs/2,tam);
% % % % fnf(2^15) = 0;
% % % figure;
% % % Zf = inv(dftmtx(N));
% % % for kk = 1:1
% % %     irf = Zf(13,:);
% % % %     irf(2^15) = 0;
% % %     IRf(kk,:) = fftshift(fft(irf));
% % %     plot( fnf, ( abs(IRf(kk,:)) ) )%/max(abs(max(IRf(kk,:)))) ) );
% % %     hold on
% % % end
% % % 
% % % 

