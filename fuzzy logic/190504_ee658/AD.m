clear all;
close all;
clc;
[Y FS]=audioread('speech1_10k.wav');
% Windowing the signal %
N=250;
Ham_Window=hamming(N);
Ham_Window_speech=Ham_Window.*Y;
% Computing the autocorrelation %
[Autocorr lags]=xcorr(Ham_Window_speech);
figure(1);
plot(lags,Autocorr);
title('Autocorrelation of windowed speech');
xlabel('lags');
ylabel('Autocorrelation magnitude');
% Building a 4-pole model %
Rp=Autocorr(250:254);
Hermition_toeplitz=toeplitz(Rp(1:1:4));
% Determining the coefficients %
A_levin=levinson(Rp,4);
Coeff=A_levin(2:1:5);
Energy=sqrt(Rp(1)-Coeff*Rp(2:1:5));
% Evaluating inverse filter A(z) %
f=[1 -1*Coeff];
M=256; % 256 point dft
A_levin_dft=fft(f,M);
% Evaluating the all-pole model H(z) %
H_dft=zeros(1,M);
for i=1:1:M,
H_dft(i)=Energy/A_levin_dft(i);
end
Log_Mag_H_dft=log(abs(H_dft));
figure(2);
%subplot(2,1,1);
plot(Log_Mag_H_dft);
title('Estimated of vocal tract response LPC');
xlabel('Samples');
ylabel('Log-Magnitude');
figure(3);
%subplot(2,1,2);
plot(log(abs(fft(Ham_Window_speech,M))));
title('Magnitude of fourier transform of windowed speech');
xlabel('Samples');
ylabel('Log-Magnitude');
% Evaluating the prediction error %
y_dft=fft(Ham_Window_speech,M);
Prediction_error_dft=y_dft'.*A_levin_dft;
Prediction_error=ifft(Prediction_error_dft);
figure(4)
plot(Prediction_error);
title('Prediction error LPC');
xlabel('Samples');
ylabel('Magnitude');
