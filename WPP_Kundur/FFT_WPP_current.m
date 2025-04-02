clear all
clc
close all
% load 'voltajes armonicos'
load caso_Zamora_01 %caso_19may %
t = simout.Time;
signals = simout.Data;
Ts = t(2) - t(1);
Fs = 1/Ts; % Sample frequency
F1 = 60;
No = round(Fs/F1); %samples per period
signals = signals(end-12*No:end,10);%(No+1:13*No+1,7);
Iabc = signals;
N = length(Iabc);%12*No;%
fun = Iabc(1:N)';
[ns, Nm] = size(fun);
Ts = t(2)-t(1);
Fs = 1/Ts;
sf = fun - mean(fun);
% sf(:,2^15) = 0;
%
f_sf = linspace(-Fs/2,Fs/2,length(sf));%/(F1);
% f_sf = Fs*(-(N/2):(N/2-1))/(N*F1);
%
Esp_sf = fftshift(fft(sf));
% spec_norm = abs(fftshift(fft(sf)/N));
% spec_norm(1:12:end) = spec_norm(1:12:end)./spec_norm(1:12:end);
spec = 2*abs( Esp_sf/N );
aspec = angle( Esp_sf );
figure;
plot(f_sf, spec);
xlim([-12*60 12*60])
title('Frequency responses of filters')
ylabel('Magnitude')
xlabel('Frequency (Hz)')
grid on

specT = spec((N+1)/2+12:12:end);
specT2 = spec((N+1)/2-12:-12:1);
spec2 = [spec((N+1)/2-5) specT2(1) specT2(5) spec((N+1)/2-12*5-6) specT2(7) specT2(11) spec((N+1)/2+5) specT(1) specT(5) spec((N+1)/2+12*5+6) specT(7) specT(11)];
aspecT = aspec((N+1)/2+12:12:end);
aspecT2 = aspec((N+1)/2-12:-12:1);
aspec2 = [aspec((N+1)/2-5) aspecT2(1) aspecT2(5) aspec((N+1)/2-12*5-6) aspecT2(7) aspecT2(11) aspec((N+1)/2+5) aspecT(1) aspecT(5) aspec((N+1)/2+12*5+6) aspecT(7) aspecT(11)];
fre_fft1 = f_sf((N+1)/2+12:12:end);
fre_fft2 = f_sf((N+1)/2-12:-12:1);
fre_fft = [f_sf((N+1)/2-5) fre_fft2(1)  fre_fft2(5) f_sf((N+1)/2-12*5-6) fre_fft2(7) fre_fft2(11) f_sf((N+1)/2+5) fre_fft1(1)  fre_fft1(5) f_sf((N+1)/2+12*5+6) fre_fft1(7) fre_fft1(11)];
% % % harmonics_FFT = spec2;%./max(spec2);
% % % % harmonics2_FFT = spec(:,2);
% % % algo = [1 5 7 11 13 17 19 23 25 29];
suma = sum(spec2(9:end).^2) + spec2(7).^2;
fundamental = max(spec2);
THD_FFT = 100*sqrt(suma) ./ fundamental

% % % figure;
% % % bar(harmonics_FFT,0.5)
% % % ylabel('% of fundamental')
% % % xlabel('Harmonic number, #')
% % % % ylim([0 5])
% % % xlim([1 12])
% % % legend('FFT')
% % % grid on
% % % set(gca, 'LooseInset', [0,0,0,0]);

I_fft = 0;
for i = 1:6
    I_fft = spec2(i+6).*cos( fre_fft(i+6)*2*pi*t(1:N) + aspec2(i+6) ) + I_fft;
end
figure;
plot(t(1:N), fun, t(1:N), I_fft)
grid on
