clear all
clc
close all
load 'I_harmonics_YN_D_V02'
F1 = 60;
signals = staa;
t = t11;
Iabc = signals(:,1);
N = length(Iabc)-1;%12*No;%
fun = Iabc(1:N)';
[ns, Nm] = size(fun);

sf = fun - mean(fun);
% sf(:,2^15) = 0;
%
f_sf = linspace(-Fs/2,Fs/2,length(sf(1,:)));
%
Esp_sf = fftshift(fft(sf));
% spec = abs( Esp_sf(1:length(f_sf)) / max(max (Esp_sf')));
spec = 2*abs( Esp_sf )/N;
figure;
plot(f_sf, spec);
xlim([-20*60 20*60])
title('Frequency responses of filters')
ylabel('Magnitude')
xlabel('Frequency (Hz)')
grid on

spec2 = spec(N/2+12:12:end);
harmonics_FFT = spec2;%100*spec2./max(spec2);
% harmonics2_FFT = spec(:,2);
suma = sum(spec2(2:end).^2);
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental

figure;
% bar(harmonics(:,3))
% hold on
bar(harmonics_FFT)
ylabel('% of fundamental')
xlabel('Harmonic number, #')
% ylim([0 5])
% xlim([1 20])
% legend('Osplines', 'ERA')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

% figure;
% plot(f_sf, spec)
% grid on


