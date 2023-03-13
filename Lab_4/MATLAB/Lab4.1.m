clear all;
close all;

%% ---------- 1 -------------
%%1.A%%

Vdc = 100;
R = 10;
L = 0.025;
f = 50;
a1 = 30*pi/180;
a2 = 18*pi/180;
T = 1/f;
w = 2*pi*f;
dt = 2*(10^(-5));

%Initialization of Currents
IQ1 = zeros(1,3001);
IQ2 = zeros(1,3001);
IQ3 = zeros(1,3001);
IQ4 = zeros(1,3001);


%Systems parameters
A_1 = -R/L;
B_1 = 1/L;
C_1 = 1;
D_1 = 0;

%System
sys_1 = ss(A_1,B_1,C_1,D_1);

%Discrete system
sys_1_dis= c2d(sys_1,dt);

k=1;
%Square pulse construction when a=30 degrees
for wt_1 = 0:dt*2*pi*f:6*pi
    wt = mod(wt_1,2*pi);
    if(wt <= a1)
        v1(k) = 0;
    elseif (wt <= pi - a1)
        v1(k) = Vdc;
    elseif (wt <= pi + a1)
        v1(k) = 0;
    elseif (wt <= 2*pi -a1)
        v1(k) =  -Vdc;
    else 
        v(k) = 0;
    end;
    k=k+1;
end;

k=1;
%Square pulse construction when a=18 degrees
for wt_1 = 0:dt*2*pi*f:6*pi
    wt = mod(wt_1,2*pi);
    if(wt <= a2)
        v2(k) = 0;
    elseif (wt <= pi - a2)
        v2(k) = Vdc;
    elseif (wt <= pi + a2)
        v2(k) = 0;
    elseif (wt <= 2*pi -a2)
        v2(k) =  -Vdc;
    else 
        v2(k) = 0;
    end;
    k=k+1;
end;



%Current when a = 30
I_L1(1) = 0;
%Current when a = 18
I_L2(1) = 0;


for k = 1:1:3000
I_L1(k+1) = sys_1_dis.A*I_L1(k) + sys_1_dis.B*v1(k);
I_L2(k+1) = sys_1_dis.A*I_L2(k) + sys_1_dis.B*v2(k);
end;


wt_1 = 0:dt*2*pi*f:6*pi;
figure();
plot(wt_1,v1);
axis([0 6*pi -120 120]);
title('Voltage and Current output when a=30°')
grid on;
hold on; 
plot(wt_1,I_L1);
legend('Vo','Io');
hold off;


wt_1 = 0:dt*2*pi*f:6*pi;
figure();
plot(wt_1,v2);
axis([0 6*pi -120 120]);
title('Voltage and Current output when a=18°')
grid on;
hold on; 
plot(wt_1,I_L2);
legend('Vo','Io');
hold off;

k=1;

%========= For a=30 ===========%

%%Current calculation
for wt_1 = 0:dt*2*pi*f:6*pi
wt = mod(wt_1,2*pi);

if (wt<=a1)
    IQ2(k) = I_L1(k);
elseif (wt<=pi)
    IQ1(k) = I_L1(k);
elseif (wt<=2*pi-a1)
    IQ4(k) = I_L1(k);
else 
    IQ3(k) = I_L1(k);
end;
k = k+1;
end;

figure;

wt_1 = 0:dt*2*pi*f:6*pi;
subplot(2,2,1);
plot(wt_1,IQ1);
axis([0 6*pi -20 20]);
title('IQ1 when a=30°');
xlabel('wt(rad/sec)');
ylabel('IQ1(A)');
grid on;

subplot(2,2,2);
plot(wt_1,IQ2);
axis([0 6*pi -20 20]);
title('IQ2 when a=30°');
xlabel('wt(rad/sec)');
ylabel('IQ2(A)');
grid on;

subplot(2,2,3);
plot(wt_1,IQ3);
axis([0 6*pi -20 20]);
title('IQ3 when a=30°');
xlabel('wt(rad/sec)');
ylabel('IQ3(A)');
grid on;

subplot(2,2,4);
plot(wt_1,IQ4);
axis([0 6*pi -20 20]);
title('IQ4 when a=30°');
xlabel('wt(rad/sec)');
ylabel('IQ4(A)');
grid on;


%========================%
%Initialization of Currents
IQ1 = zeros(1,3001);
IQ2 = zeros(1,3001);
IQ3 = zeros(1,3001);
IQ4 = zeros(1,3001);


%Current when a = 30
I_L1(1) = 0;
%Current when a = 18
I_L2(1) = 0;

k=1;
for k = 1:1:3000
I_L1(k+1) = sys_1_dis.A*I_L1(k) + sys_1_dis.B*v1(k);
I_L2(k+1) = sys_1_dis.A*I_L2(k) + sys_1_dis.B*v2(k);
end;

%=========== For  a=18 ============%

k=1;
%%Current calculation
for wt_1 = 0:dt*2*pi*f:6*pi
wt = mod(wt_1,2*pi);

if (wt<=a2)
    IQ2(k) = I_L2(k);
elseif (wt<=pi)
    IQ1(k) = I_L2(k);
elseif (wt<=2*pi-a2)
    IQ4(k) = I_L2(k);
else
    IQ3(k) = I_L2(k);
end;
k = k+1;
end;

figure;

wt_1 = 0:dt*2*pi*f:6*pi;
subplot(2,2,1);
plot(wt_1,IQ1);
axis([0 6*pi -20 20]);
title('IQ1 when a=18°');
xlabel('wt(rad/sec)');
ylabel('IQ1(A)');
grid on;

subplot(2,2,2);
plot(wt_1,IQ2);
axis([0 6*pi -20 20]);
title('IQ2 when a=18°');
xlabel('wt(rad/sec)');
ylabel('IQ2(A)');
grid on;

subplot(2,2,3);
plot(wt_1,IQ3);
axis([0 6*pi -20 20]);
title('IQ3 when a=18°');
xlabel('wt(rad/sec)');
ylabel('IQ3(A)');
grid on;

subplot(2,2,4);
plot(wt_1,IQ4);
axis([0 6*pi -20 20]);
title('IQ4 when a=18°');
xlabel('wt(rad/sec)');
ylabel('IQ4(A)');
grid on;


%================= 1.B Fourier Transform ==============%



 NFFT = 3000;
 Fs = 50000;
 V1 = fft(v1,NFFT)/3000; 
 I1 = fft(I_L1,NFFT)/3000;
 V2 = fft(v2,NFFT)/3000; 
 I2 = fft(I_L2,NFFT)/3000; 
 f = Fs/2*linspace(0,1,NFFT/2+1);

figure;
stem(f,2*abs(V1(1:NFFT/2+1)))
title('Voltage Harmonics (V1)')
xlabel('Frequency (Hz)')
ylabel('Voltage(V)')

figure;
stem(f,2*abs(V2(1:NFFT/2+1)))
title('Voltage Harmonics (V2)')
xlabel('Frequency (Hz)')
ylabel('Voltage(V)')


%================ Active Power Calculation ===============%

p1=0;
p2=0;

%Instant Active Power
k=0;
for k = 1:1:3001
    P1(k)= v1(k)*I_L1(k);  
    P2(k)= v2(k)*I_L2(k);
end

figure;
plot(wt_1,P1)
title('P1(t)');
xlabel('wt(rad/sec)')
ylabel('instant power (W)')

figure;
plot(wt_1,P2)
title('P2(t)');
xlabel('wt(rad/sec)')
ylabel('instant power (W)')

%Active Power 

for k = 1:1:3001
    p1= P1(k)+ p1;
    p2= P2(k)+ p2;
    
end

p1 = p1/3001
p2 = p2/3001

%===== RMS values calculation for Power Factor Calculation =======%

%Initialization

k=0;
V1_num=0;
V2_num=0;
I1_num=0;
I2_num=0;

for k= 1:1:3000
    V1_num = ((rms(V1(k))^2))+V1_num;
    V2_num = ((rms(V2(k))^2))+V2_num;
    I1_num = ((rms(I1(k))^2))+I1_num;
    I2_num = ((rms(I2(k))^2))+I2_num;
end

V1_rms=sqrt(V1_num);
V2_rms=sqrt(V2_num);
I1_rms=sqrt(I1_num);
I2_rms=sqrt(I2_num);

s1=V1_rms*I1_rms
s2=V2_rms*I2_rms

%Power Factor
PF_1= p1/s1
PF_2= p2/s2




