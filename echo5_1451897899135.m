clc;
clear all;
%warning off;
close all;
 
% Input signals 

load farspeech
ss=x;
%ss= audioread('C:\Users\OM\Desktop\clean\sp12.wav');
subplot 411;
plot(ss);

n=audioread('echoedsignal.wav');
subplot 412;
plot(n);% subplot (1,2,2), disp (n);
load nearspeech;

subplot 413;
plot(v);
ss1=n-ss;
for i=1:1:20115
    a(i)=v(i)+ss1(i);
    s(i)=ss1(i);
end
%wavwrite(a,'C:\Users\OM\Desktop\d')
subplot 414;
plot(a);
% plot(a);
%Initialization
 N=20115;
 %Hpsd=dspdata.psd(N);
p=1024;
 mu=0.001;
w=zeros(p,1);
% x=zeros(N,1);
% d=zeros(N,1);
 
%  figure;
%  plot(s);
%Algorithm
for i=p:N
    xvec=n(i:-1:i-p+1);
    y(i)=w'*xvec;
    e(i)=a(i)-y(i);
    w=w+mu*e(i)*xvec;
end

% plot((s'-e));
% xlabel('time index'); 
% ylabel('signal value');
%Calculating MSE
[R C]=size(s); 
err = (((s-e).^2)/(R*C)); % size of signals must be same length 
MSE=sqrt(err); 

% for i=1:N
%      err(i)=(s(i)-e(i)).^2;
%     % nn(i)=n(i).^2;
% %     err(i)=(s(i)-e(i));
% %      nn(i)=n(i);
% %     ss(i)=s(i);
% end
% MSE=sqrt(err)
 plot(MSE);
%  Calculating SNR INPUT
% SNR=snr(e,(ss-e))
 rms_signal=sqrt(mean(a.^2));
  rms_echo=sqrt(mean((s-a).^2));
  Lsig=10*log10(rms_signal);
  Lech=10*log10(rms_echo);
  ERLEo=Lsig-Lech

%  
% Calculating SNR OUTPUT
% SNR=snr(e,(ss-e))
 rms_signal=sqrt(mean(e.^2));
  rms_echo=sqrt(mean((s-e).^2));
  Lsig=10*log10(rms_signal);
  Lech=10*log10(rms_echo);
  ERLEi=Lsig-Lech