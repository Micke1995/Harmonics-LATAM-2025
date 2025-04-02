clear all
clc
close all
load 'Irectifier'
F1 = 60;
Ts = 1e-4;
Fs = 1/Ts; % Sample frequency
No = round(Fs/F1); %samples per period
signals = Irs(1:12*No,1);
Iabc = signals;
N = length(Iabc);%12*No;%
fun = Iabc(1:N)';
[ns, Nm] = size(fun);
Ts = t(2)-t(1);
Fs = 1/Ts;
sf = fun - mean(fun);
% sf(:,2^15) = 0;
%
% f_sf = linspace(-Fs/2,Fs/2,length(sf));
f_sf = Fs*(-(N/2):(N/2-1))/N;
%
Esp_sf = fftshift(fft(sf));
% spec_norm = abs(fftshift(fft(sf)/N));
% spec_norm(1:12:end) = spec_norm(1:12:end)./spec_norm(1:12:end);
spec = 2*abs( Esp_sf/N );
aspec = angle( Esp_sf );
figure;
plot(f_sf, spec);
xlim([-30*60 30*60])
title('Frequency responses of filters')
ylabel('Magnitude')
xlabel('Frequency (Hz)')
grid on

spec2 = spec(N/2+13:12:end);
aspec2 = aspec(N/2+13:12:end);
fre_fft = f_sf(N/2+13:12:end);
harmonics_FFT = spec2;%./max(spec2);
% harmonics2_FFT = spec(:,2);
algo = [1 5 7 11 13 17 19 23 25 29];
suma = sum(spec2(2:end).^2);
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental

figure;
bar(harmonics_FFT,0.5)
ylabel('% of fundamental')
xlabel('Harmonic number, #')
% ylim([0 5])
% xlim([1 20])
legend('FFT')
grid on
set(gca, 'LooseInset', [0,0,0,0]);

I_fft = 0;
for i = 1:length(algo)
    I_fft = spec2(algo(i)).*cos( fre_fft(algo(i))*2*pi*t(1:N) + aspec2(algo(i)) ) + I_fft;
end
figure;
plot(t(1:N), fun, t(1:N), I_fft)
grid on
