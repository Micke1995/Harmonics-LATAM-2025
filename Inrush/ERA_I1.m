% function [Frec, Amp,theta,thd] = ERA(N,k,Ts,fun)
clear all
clc
load 'I_harmonics_YN_D_V02'
% % % [z,p] = butter(2,360/Fs,'low');
% % % ia = filter(z, p, st(:,1));
% % % ib = filter(z, p, st(:,2));
% % % ic = filter(z, p, st(:,3));
% % % st = [ia ib ic];
F1 = 60;
% Fs = Fs; % Sample frequency
% No = Fs/F1; %samples per period
% Ts = 1/Fs;
%load CB25kV2
%load AEMC
%load Vc
% load VabcB25kW
signals = staa;%st(:,1) - mean(st(:,1));
% signals = downsample(VabcB25kW,2);
t = t11;%signals(0.8*Fs:end,1); %
Iabc = signals(:,1);%signals(0.8*Fs:end,2);%
N = length(Iabc)-1;%12*No;%
fun = Iabc(1:N)';
fun2 = Iabc(1:N+1)';
k = 104;
r = round(N/2) - 1;
%L = ceil((N) / 2)
%for i=1:(N-L)
 %   H0(i,:)=fun(i:(i+L));
%end
%for i=2:(N-L)
 %   H1(i,:)=fun(i:i+L);
%end
H0 = hankel(fun(1:r),fun(r:N-2));
H1 = hankel(fun(2:r+1),fun(r+1:N-1));
[U,S,V] = svds(H0,k);
A = (S^-(1/2))*U'*H1*V*(S^-(1/2));
z = eig(A);
% % V = zeros();
% % for i = 1:N
% %     for j = 1:k
% %         V(i,j) = z(j)^(i-1);
% %     end
% % end
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
% rows3 = [3,6,9,12,14,16,17,20,22,24,26,28,29,31,32,35,36,37,38];%
rows3 = [2 4 5 7 8 9 11 12 13 14 15 17 19 20 21 23];
%find(ampli2>0);
% % % ampli3 = ampli2(rows3);
% % % freq3 = ceil(freq2(rows3));
% % % phase3 = phase2(rows3);
% % % dampr3 = dampr2(rows3);
% % % mfreq = 60.*(1:length(freq3))';
% % % ce = 0;
% % % for k = 1:length(ampli3)
% % %     rowss = find(freq3 == mfreq(k));
% % %     for m = 1:length(rowss)
% % %     rows4(k+ce) = rowss(m);
% % %     ce = ce + 1;
% % %     end
% % % end
% % % rows5 = find(rows4>0);
% % % rows6 = rows4(rows5);
% % % ampli4 = ampli3(rows6);
% % % freq4 = freq3(rows6);
% % % phase4 = phase3(rows6);
% % % dampr4 = dampr3(rows6);


 disp('------------------ERA------------------')
disp('    Freq       Amp       Phase    DampR')
% apli5 = ampli3;%100*ampli4./max(ampli4);
todo22 = todo_orden(rows3,:);
todo3 = zeros(17,4);
todo2(1:13,:) = todo22(1:13,:);
todo2(14,:) = [14*60 0 0 0];
todo2(15:17,:) = todo22(14:16,:);

disp(sortrows(todo2))

% IHD(:,k) = 100*(amp_h(2:end,k))./(amp_h(1,k));
harmonics_ERA = 100*todo2(:,2)./max(todo2(:,2));
% harmonics2_ERA = todo2(:,2);
suma = sum(todo2(2:end,2).^2);
fundamental = max(todo2(:,2));
THD_ERA = 100*sqrt(suma) ./ fundamental

%% FFT
load harmonics_FFT_1_2
% f_sf = Fs*(-(N/2):(N/2-1))/N;
aspec = angle(Esp_sf);
ph_fft = aspec(N/2+12:12:end);
fre_fft = f_sf(N/2+12:12:end)';
amp_fft = spec2;

harmonics_FFT = 100*spec2./max(spec2);
% harmonics2_FFT = spec2(:,2);
suma = sum(spec2(2:end).^2);
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental

figure;
harmonics = [harmonics_FFT(1:17)' harmonics_ERA];
bar(harmonics)
ylabel('% of fundamental')
xlabel('Harmonic number, #')
% ylim([0 5])
xlim([1.5 17.5])
legend('FFT', 'ERA')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

I_era = 0;
I_fft = 0;
for i = 1:17
    I_era = todo2(i,2).*cos(todo2(i,1)*2*pi*t(1:N) + todo2(i,3) - pi/4 ) + I_era;
    I_fft = amp_fft(i).*cos(fre_fft(i)*2*pi*t(1:N) - ph_fft(i) - pi/4 - pi/2 ) + I_fft;
end

error_fft = fun' - I_fft;
error_era = fun' - I_era;

figure;
plot(t(1:N), fun, t(1:N), I_fft, '-.', t(1:N), I_era, '--')
ylabel('Current (A)')
xlabel('Time (s)')
xlim([5 5.2])
legend('Actual', 'FFT', 'ERA')
grid on
set(gca, 'LooseInset', [0,0,0,0]);


RMSE_FFT = sqrt(mean((fun' - I_fft).^2)) 
RMSE_ERA = sqrt(mean((fun' - I_era).^2))

a = norm(fun',2)^2;
b_fft = norm(fun,2)^2 - norm(I_fft,2)^2;
b_era = norm(fun,2)^2 - norm(I_era,2)^2;
SNR_FFT = 10*log10(a/b_fft);
SNR_ERA = 10*log10(-a/b_era);
SSE_FFT = norm((fun' - I_fft),2)^2;
% SSE_FFT = sum((fun' - I_fft).^2);
SSE_ERA = sum((fun' - I_era).^2);
MSE_FFT = (1/N)*SSE_FFT;
MSE_ERA = (1/N)*SSE_ERA;
RMSE_fft = sqrt(MSE_FFT);
RMSE_era = sqrt(MSE_ERA);
Errors_FFT = [SNR_FFT ; SSE_FFT ; MSE_FFT ; RMSE_FFT]
Errors_ERA = [SNR_ERA ; SSE_ERA ; MSE_ERA ; RMSE_ERA]


%% Frequency Response
tam = N;%2^15;
fn = Fs*(-tam/2:tam/2)/(tam*60);%linspace(-Fs/2,Fs/2,tam);
figure;
Zp = pinv(ZZ);
for k = 1:16
    algo1(k) = find(Frec == todo22(k,1));
end
algo = sortrows(algo1',1);
algo(17:32) = algo + 1;
amp_era = [todo2(end:-1:7,2) ; todo2(5:6,2) ; todo2(4:-1:3,2) ; todo2(1:2,2) ; todo2(end:-1:7,2) ; todo2(5:6,2) ; todo2(4:-1:3,2) ; todo2(1:2,2)];
amp_era([4 21]) = [];
subaxis(2, 1, 1, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.12, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.2)
for kk = 1:32
    ir = Zp(algo(kk),:);
%     ir(tam) = 0;
    IR(kk,:) = fftshift(fft(ir));
    plot( fn, amp_era(kk)*abs( IR(kk,:) ));%/max(abs(max(IR(kk,:)))) ) );%amp_era(kk)*
    hold on
end
% plot(fn,1*abs(IR))%/max(abs(max(IR))));
% % % for kk = 3:length(algo)
% % %     ir = Zp(algo(kk),:);
% % % %     ir(tam) = 0;
% % %     IR = fftshift(fft(ir));
% % %     plot(fn,1*abs(IR));
% % %     hold on
% % % end
% ylim([0 1.1])
xlim([-18 18])
% title('Frequency responses of filters')
set(gca,'XTickLabel',[])
ylabel('Magnitude')
% xlabel('Normalised Frequency u=f_1T_1')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

subaxis(2, 1, 2, 'Spacing', 0.1, 'MarginRight', 0.02, 'MarginLeft', 0.12, 'SpacingVert', 0.08, 'MarginTop', 0.04, 'MarginBottom', 0.2)
plot(f_sf/60, spec);
xlim([-18 18])
% title('Frequency responses of filters')
% ylabel('Magnitude')
xlabel('Normalised Frequency u=f_1T_1')
grid on

figure;
plot(t(1:N), fun)
ylabel('Current (A)')
xlabel('Time (s)')
xlim([5 5.2])
grid on
set(gca, 'LooseInset', [0,0,0,0]);





