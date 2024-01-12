 % Coded by Vaaiibhav
clc;
clear all;
close all;
%%  room impulse response

%% Select ROom Impulse Response Wave file
[Fname3,Pname3] = uigetfile('RIR.wav','Select RIR FIle'); 
[h,Fs3] = audioread(strcat(Pname3,Fname3));    %Far-end signal 

H = length(h)/62;
%% 

fs = 8000;
M = fs/2 + 1;
%% nearspeech
H=32;
% load nearspeech;
% v=wavread('rashnearend.wav');
[Fname1,Pname1] = uigetfile('*.wav','Select nearspeech FIle'); 


[v,Fs1] = audioread(strcat(Pname1,Fname1));

near = v;
figure(1);
plot(near);
title('nearspeech');
nearlength = length(near);

%% farspeech

%  load farspeech;
% x=wavread('rashfarend.wav');
[Fname2,Pname2] = uigetfile('*.wav','Select farspeech FIle'); 


[x,Fs2] = audioread(strcat(Pname2,Fname2));

far = x;
figure(2);
plot(far);
title('farspeech');
farlength = length(far);
farfilterd = filter(H,1,x);

%% creating nearplusfar 

nearplusfar = v + x;
nearplusfar=2*nearplusfar/(max(nearplusfar)-min(nearplusfar));
nearplusfar=nearplusfar-min(nearplusfar)-1;
figure (3);
plot (nearplusfar);
% sound(nearplusfar,8000);
title('nearplusfar');

%% far and echoed speech
x = x(1:length(x));
% sound(x, 8000);

dhat = filter(H,1,x);
d = dhat + v+0.001*randn(length(v),1);
figure (8);
plot (d);
figure (11);
plot (dhat);
title('farandechoed'); 
%  sound(dhat,8000);
%  sound(d,8000);

%% 
pause(30);


micSignal = d + v +0.001*randn(length(v),1);
               
figure (5);
plot (micSignal);
title('micSignal'); 
% p8 = audioplayer(micSignal,fs);
%  playblocking(p8);



% playthis
%  sound(micSignal,8000);


%% LMS Algo for removing echo
p = 2;

N = farlength;
%modu = abs(mod(x(i)));
%mu = beta(x(i)/modu^2)*e(i);
mu=0.12;
alfa=0.002;
c=0.001;
w=zeros(1,N);

% x = x(1:length(W0)*floor(length(x)/length(W0)));
% d = d(1:length(W0)*floor(length(d)/length(W0)));

% lmsout =  lms(micSignal,nearplusfar,mu,1,12);\%Algorithm
% for i=p:N
%     xvec=nearplusfar(i:-1:i-p+1);
%     y(i)=w'*xvec;
%     e(i)=micSignal(i)-y(i);
%     w=w+mu*e(i)*xvec;
% end



for i=1:N
%     p= logical(p);
%     p = round(p);
%    x(i)=d(i:-1:i-p+1);

tic
    %e(i)=1;
  
mu=alfa/(c+(x(i)'*x(i)));
 
   e(i) = d(i) - w(i)' * x(i);
   %ce(i)=abs(e(i));
   ce(i)=e(i);
 %  w(i+1) = w(i) + (mu * ce(i)) * x(i);
   w(i+1) = w(i) + mu * ce(i) * x(i);
   
  
 toc
 timeForEachIteration(i) = toc;
end
for i=1:N
yd(i) = sum(w(i)' * x(i));  

end

figure(6);
plot(e);
title('NLMS Out');
% p8 = audioplayer(e,fs);
%  playblocking(p8);
  sound(e,8000);
[R C]=size(x);

%% MSE 
error = e;
se = (e.^2);
mse= se;

figure(8);
plot(abs(mse));
title('MSE Graph');



%% Plotting erle

Hd2 = dfilt.dffir(ones(1,1000));

e=transpose(e);
erle = filter(Hd2,(e-v(1:length(e))).^2)./ ...
    (filter(Hd2,micSignal(1:length(e)).^2));
erledB = 10*log10(erle);
erledB = abs(erledB);
figure(7);
plot(erledB);
xlabel('Samlpes]');
ylabel('ERLE [dB]');
title('Echo Return Loss Enhancement');
set(gcf, 'Color', [1 1 1])
maxer=max(erledB)
disp(maxer);


plot(timeForEachIteration);
title('Convergence Graph');
                