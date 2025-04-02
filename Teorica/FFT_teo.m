clear all
clc
close all
Fs = 4*1920;
F1 = 60;
Ts = 1/Fs;
t = 0:Ts:1;
phi = 1 - exp(-t/1).*cos(2*pi*t/5); % pi/4;%
xt = cos(2*pi*F1*t + phi);
st = tanh(3*xt);
No = round(Fs/F1); %samples per period
signals = st(1:12*No);
Iabc = signals;
N = length(Iabc);%12*No;%
fun = Iabc(1:N)';
[ns, Nm] = size(fun);
f_sf = Fs*(-(N/2):(N/2-1))/N;
sf = fun;
%
Esp_sf = fftshift(fft(sf));
spec = 2*abs( Esp_sf/N );
aspec = angle( Esp_sf );
figure;
plot(f_sf, spec);
% xlim([-30*60 30*60])
title('Frequency responses of filters')
ylabel('Magnitude')
xlabel('Frequency (Hz)')
grid on

spec2 = spec(N/2+13:12:end);
aspec2 = aspec(N/2+13:12:end);
fre_fft = f_sf(N/2+13:12:end);
harmonics_FFT = spec2;%./max(spec2);
suma = sum(spec2(2:end).^2);
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental
save otro aspec2 spec2 spec fre_fft 
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
for i = 1:13 %length(spec2)
    I_fft = spec2(i).*cos( fre_fft(i)*2*pi*t(1:N) + aspec2(i) ) + I_fft;
end
figure;
plot(t(1:N), fun, t(1:N), I_fft)
grid on
RMSE_FFT = sqrt(mean((fun' - I_fft).^2)) 


