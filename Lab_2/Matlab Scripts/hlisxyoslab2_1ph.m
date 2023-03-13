clear all;
close all;
clc;

%% 1 phase

Vrms=230;       %Vrms
f=50;           %frequency
dt=0.00002;     %step
T=10;           %time duration
t=0:dt:T-dt;    %sampling time
t1=0:dt:0.02;   %first period

w=0:360/1000:(T/0.02)*360-(360/1000);   %angle in rad 2*pi*f*t*180/pi

V1ph = 230*sqrt(2)*sin(2*pi*f*t);

figure
plot(w,V1ph)
hold on;
plot(w,-V1ph,'r')
hold off;
axis([0 360 -400 400])
legend('V','-V')
xlabel('Angle (degrees)')
title('1-ph input Voltage (V)')
 
 
% load voltage
i=1;
for t=0:dt:T-dt
    if (mod(2*pi*f*t,2*pi) < pi) 
        V_L(i) = 230*sqrt(2)*sin(2*pi*f*t);
    else
        V_L(i) = 230*sqrt(2)*sin(2*pi*f*t+pi);
    end
    i=i+1;
end
    
 
figure
plot(w,V_L)
axis([0 720 -350  350])
xlabel('Angle (degrees)')
ylabel('Voltage (v)')
title('1-phase load voltage ,a=0')


a = 0;     %firing angle
R = 2.5;   %resistance of loads 
L1 = 0.04; %load 1
L2 = 0.08; %load 2

%L=0.04
A1 = -R/L1;
B1 = 1/L1;
C1 = 1;

%L=0.08 
A2 = -R/L2;
B2 = 1/L2;
C2 = 1;

sys1=ss(A1,B1,C1,0);
sys2=ss(A2,B2,C2,0);

sys1_d=c2d(sys1,dt);
sys2_d=c2d(sys2,dt);

%creating matrices A,B,C (D is zero)
I_L1 = 0:1:500000-1;
I_L2 = 0:1:500000-1;

Ad1 = sys1_d.A;
Bd1 = sys1_d.B;
Cd1 = sys1_d.C;

Ad2 = sys2_d.A;
Bd2 = sys2_d.B;
Cd2 = sys2_d.C;

% load current
for k = 1:1:500000-1
    I_L1(k+1) = Ad1*I_L1(k) + Bd1*V_L(k);   % current for L= 0.04 H
    I_L2(k+1) = Ad2*I_L2(k) + Bd2*V_L(k); % current for L= 0.08 H
end;


figure
plot(w,I_L1)
axis([18000 18360 0  90]);
title('1-ph Load Current a=0, L=0.04')
xlabel('Angle (degrees)')
ylabel('Current (A)')

figure
plot(w,I_L2)
axis([18000 18360 0  90]);
title('1-ph Load Current a=0, L=0.08')
xlabel('Angle (degrees)')
ylabel('Current (A)')

 
%Thyristors currents
i=1;
 for t=0:dt:T-dt
    if (mod(2*pi*f*t,2*pi) < pi) 
        I12a(i) = I_L1(i);   % L =0.04 H
        I12b(i)= I_L2(i);   % L =0.08 H
        I34a(i) = 0;
        I34b(i)= 0;
    else
        I12a(i) = 0;
        I12b(i)= 0;
        I34a(i) = I_L1(i);  % L =0.04 H
        I34b(i)= I_L2(i);  % L =0.08 H
    end
    i=i+1;
 end
 
 figure
 subplot(1,2,1);
 plot(w,I12a)
 title('Current of thyristors T1,T2 for L = 0.04 for a=0');
 axis([18000 18360 -5  90]);
 
 subplot(1,2,2);
 plot(w,I34a)
 title('Current of thyristors T3,T4 for L = 0.04 for a=0');
 axis([18000 18360 -5  90]);
 
 figure
 subplot(1,2,1);
 plot(w,I12b)
 title('Current of thyristors T1,T2 for L = 0.08 for a=0');
 axis([18000 18360 -5  90]);
 
 subplot(1,2,2);
 plot(w,I34b)
 title('Current of thyristors T3,T4 for L = 0.08 for a=0');
 axis([18000 18360 -5  90]);

 
 % load voltage
  i=1;
 for t=0:dt:T-dt  
    if (mod(2*pi*f*t,2*pi) < pi) 
        if(mod(2*pi*f*t,2*pi) < pi/2)
           V_L(i)=230*sqrt(2)*sin(2*pi*f*t+pi) ;
        else
           V_L(i)=230*sqrt(2)*sin(2*pi*f*t);
        end
    else
        if(mod(2*pi*f*t,2*pi) < 3*pi/2)
           V_L(i)=230*sqrt(2)*sin(2*pi*f*t) ;
        else
           V_L(i)=230*sqrt(2)*sin(2*pi*f*t+pi);
        end
    end
    i=i+1;
 end
 
 figure
 plot(w,V_L)
 axis([0 360 -350  350])
 title('1-ph Load Voltage,a=90')
 xlabel('Angle (degrees)')
 ylabel('Voltage (v)')
 
 % load current
 t=0:dt:T-dt;


for k = 1:1:500000-1    
    I_L1(k+1) = Ad1*I_L1(k) + Bd1*V_L(k);
    if I_L1(k+1)<0
       I_L1(k+1) =0;
    end
    I_L2(k+1) = Ad2*I_L2(k) + Bd2*V_L(k);
    if I_L2(k+1)<0
       I_L2(k+1) =0;
    end
end;

figure
plot(w,I_L1)
axis([18000 18720 -10  100]);
title('1-ph Load Current,a=90,L=0.04')
xlabel('Angle (degrees)')
ylabel('Current (A)')

figure
plot(w,I_L2)
axis([18000 18720 -10  100]);
title('1-ph Load Current,a=90,L=0.08')
xlabel('Angle (degrees)')
ylabel('Current (A)')

%Thyristors current
i=1;
 for t=0:dt:T-dt
    if (mod(2*pi*f*t,2*pi) > pi/2 && mod(2*pi*f*t,2*pi) < 3*pi/2) 
        I12a(i) = I_L1(i);  % L =0.04 H
        I12b(i)= I_L2(i); % L =0.08 H
        I34a(i) = 0;
        I34b(i)= 0;
    else
        I12a(i) = 0;
        I12b(i)= 0;
        I34a(i) = I_L1(i);  % L =0.04 H
        I34b(i)= I_L2(i); % L =0.08 H
    end
    i=i+1;
 end
 
 figure
 subplot(1,2,1);
 plot(w,I12a)
 title('Current of thyristors T1,T2 for L = 0.04 a=90');
 axis([18000 18360 -5 100]);
 
 subplot(1,2,2);
 plot(w,I34a)
 title('Current of thyristors T3,T4 for L = 0.04 a=90');
 axis([18000 18360 -5 100]);
 
 figure
 subplot(1,2,1);
 plot(w,I12b)
 title('Current of thyristors T1,T2 for L = 0.08 a=90');
 axis([18000 18360 -5 100]);
 
 subplot(1,2,2);
 plot(w,I34b)
 title('Current of thyristors T3,T4 for L = 0.08 a=90');
 axis([18000 18360 -5 100]);
