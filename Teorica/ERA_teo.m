% function [Frec, Amp,theta,thd] = ERA(N,k,Ts,fun)
clear all
clc
close all
Fs = 4*1920;
F1 = 60;
Ts = 1/Fs;
t = 0:Ts:1;
phi = 1 - exp(-t/1).*cos(2*pi*t/5); % pi/4; %
xt = cos(2*pi*F1*t + phi);
st = tanh(3*xt);
No = round(Fs/F1); %samples per period
figure;
plot(t(1:12*No), st(1:12*No))
xlabel('Time (s)')
ylabel('XXXX (A)')
% legend('i_a', 'i_b', 'i_c')
grid on
set(gca, 'LooseInset', [0,0,0,0]);
%
Iabc = st(1:12*No+1);
N = length(Iabc)-1;
fun = Iabc(1:N)';
fun2 = Iabc(1:N+1);
k = 14;
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

harmonics_ERA = 100*todo2(:,2)./max(todo2(:,2));
suma = sum(todo2(2:end,2).^2);
fundamental = max(todo2(:,2));
THD_ERA = 100*sqrt(suma) ./ fundamental

%% FFT
load FFT_teo_harmonics4
harmonics_FFT = 100*spec2./max(spec2);
suma = sum(spec2(2:end).^2);
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental
%
figure;
algo = 1:2:k-1;
harmonics = [harmonics_FFT(algo) harmonics_ERA];
bar(harmonics)
ylabel('% of fundamental')
xlabel('Harmonic number, #')
% ylim([0 5])
% xlim([-0.5 30])
legend('FFT', 'ERA')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

I_era = 0;
I_fft = 0;
for i = 1:k/2
    I_era = todo2(i,2).*cos(todo2(i,1)*2*pi*t(1:N) + todo2(i,3) ) + I_era;
    I_fft = spec2(algo(i)).*cos( fre_fft(algo(i))*2*pi*t(1:N) + aspec2(algo(i)) ) + I_fft;
    I_eram(i,:) = todo2(i,2).*cos(todo2(i,1)*2*pi*t(1:N) + todo2(i,3) );
    I_fftm(i,:) = spec2(algo(i)).*cos( fre_fft(algo(i))*2*pi*t(1:N) + aspec2(algo(i)) );
    amp_fft(i,1) = spec2(algo(i));
    amp_era(i,1) = todo2(i,2);
    fre_fft2(i,1) = fre_fft(algo(i));
    fre_era(i,1) = todo2(i,1);
end
amp_fft_era = [amp_fft amp_era]
fre_fft_era = [fre_fft2 fre_era]
figure;
plot(t(1:N), fun', t(1:N), I_fft, '-.', t(1:N), I_era, '--')
ylabel('s(t)')
xlabel('Time (s)')
legend('Actual', 'FFT', 'ERA')
ax1 = gca; % Store handle to axes 1.
% Create smaller axes in top right, and plot on it
% Store handle to axes 2 in ax2.
ax2 = axes('Position',[0.2 0.65 0.25 0.15]);
box on;
plot(t(1:N), fun', t(1:N), I_fft, '-.', t(1:N), I_era, '--')
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
xlim([t(1012) t(1036)])
ylim([0.99 0.9961])
grid on
set(gca, 'LooseInset', [0,0,0,0]);

figure;
subaxis(4, 1, 1, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.1)
plot(t(1:N), I_fftm(1,:), t(1:N), I_eram(1,:), '--')
ylabel('s_1(t)')
grid on
subaxis(4, 1, 2, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.1)
plot(t(1:N), I_fftm(2,:), t(1:N), I_eram(2,:), '--')
ylabel('s_3(t)')
grid on
subaxis(4, 1, 3, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.1)
plot(t(1:N), I_fftm(3,:), t(1:N), I_eram(3,:), '--')
ylabel('s_5(t)')
grid on
subaxis(4, 1, 4, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.1)
plot(t(1:N), I_fftm(4,:), t(1:N), I_eram(4,:), '--')
ylabel('s_7(t)')
xlabel('Time (s)')
grid on
% % % subaxis(5, 1, 5, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.12, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.1)
% % % plot(t(1:N), fun', t(1:N), I_fft, '-.', t(1:N), I_era, '--')
% % % ylabel('s(t), $\hat{s}(t)$')
% % % xlabel('Time (s)')
% % % grid on
% % % ax1 = gca; % Store handle to axes 1.
% % % % Create smaller axes in top right, and plot on it
% % % % Store handle to axes 2 in ax2.
% % % ax2 = axes('Position',[0.7 0.15 0.15 0.15]);
% % % box on;
% % % plot(t(1:N), fun', t(1:N), I_fft, '-.', t(1:N), I_era, '--')
% % % set(gca,'XTickLabel',[])
% % % set(gca,'YTickLabel',[])
% % % xlim([t(1012) t(1036)])
% % % ylim([0.99 0.9961])
% % % grid on
% % % set(gca, 'LooseInset', [0,0,0,0]);


RMSE_FFT = sqrt(mean((fun' - I_fft).^2)) 
RMSE_ERA = sqrt(mean((fun' - I_era).^2))

a = norm(fun',2)^2;
b_fft = norm(fun,2)^2 - norm(I_fft,2)^2;
b_era = norm(fun,2)^2 - norm(I_era,2)^2;
SNR_FFT = 10*log10(a/b_fft);
SNR_ERA = 10*log10(a/b_era);
SSE_FFT = sum((fun' - I_fft).^2);
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
amp_era = [todo2(7,2) ; todo2(7,2) ; todo2(6,2) ; todo2(6,2) ; todo2(5,2) ; todo2(5,2) ; todo2(1,2) ; todo2(1,2) ; todo2(4,2) ; todo2(4,2) ; todo2(2,2) ; todo2(2,2) ; todo2(3,2) ; todo2(3,2)];
figure;
Zp = pinv(ZZ);
% figure;
subaxis(3, 1, 1, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.12, 'MarginTop', 0.04, 'MarginBottom', 0.15)
plot(t(1:N), fun')
ylabel('s(t)')
xlabel('Time (s)')
grid on
subaxis(3, 1, 2, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.12, 'MarginTop', 0.04, 'MarginBottom', 0.15)
for kk = 1:k
    ir = Zp(kk,:);
    IR(kk,:) = fftshift(fft(ir));
    plot( fn, amp_era(kk)*abs( IR(kk,:)/max(abs(max(IR(kk,:)))) ) ); %
    hold on
end
xlim([-11 11])
set(gca,'XTickLabel',[])
% title('Frequency responses of filters')
ylabel('Magnitude')
% xlabel('Normalised Frequency u=f_1T_1')
grid on
set(gca, 'LooseInset', [0,0,0,0]);
f_sf = Fs*(-(N/2):(N/2-1))/N;
subaxis(3, 1, 3, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.1, 'SpacingVert', 0.12, 'MarginTop', 0.04, 'MarginBottom', 0.15)
plot(f_sf/F1, spec);
xlim([-11 11])
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

