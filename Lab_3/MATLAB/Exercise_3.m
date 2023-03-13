close all;
clear all;
clc;

%% Square pulse generation
Fs = 10000;         
T = 6*1/Fs;
t = 0:1/(100*Fs):T-1/(100*Fs);

pulse_tr = 150*sawtooth(2*pi*Fs*t,0.5)+150;
pulse_1 = 150*ones(600,1);
pulse_2 = 240*ones(600,1);


for i = 1:1:600
    %DT=0.5
    if pulse_tr(i)>= pulse_1(i)
        Vdc_1(i)= 0;
    else
        Vdc_1(i)= pulse_1(i);
    end
    
    %DT=0.8
    if pulse_tr(i)>= pulse_2(i)
        Vdc_2(i)= 0;
    else
        Vdc_2(i)= pulse_2(i);
    end    
end

figure;
subplot(2,1,1);
plot(t,pulse_tr)
hold on;
plot(t,Vdc_1,'r')
title('Triangular Pulse DT=0.5');
xlabel('Time(sec)');
ylabel('Volt(V)');

subplot(2,1,2);
plot(t,pulse_tr)
hold on;
plot(t,Vdc_2,'r')
title('Triangular Pulse DT=0.8');
xlabel('Time(sec)');
ylabel('Volt(V)');

%% DT = 0.5, L = 0.001

%Parameters 
R = 3;          %Resistance
L = 0.0005;     %Inductor
C = 10^(-5);    %Capacitor
L1 = 10^(-3);   %Inductor A 
L2 = 10^(-2);   %Inductor B 
L3 = 10^(-3);   %Inductor C 

%State space model
%Matrixes
A1=[0 0 -1/L1 ; 0 -R/L 1/L ; 1/C -1/C 0 ];
B1=[1/L1 0 ; 0 -1/L ; 0 0];

A2=[0 0 -1/L2 ; 0 -R/L 1/L ; 1/C -1/C 0 ];
B2=[1/L2 0 ; 0 -1/L ; 0 0];

A3=[0 0 -1/L3 ; 0 -R/L 1/L ; 1/C -1/C 0 ];
B3=[1/L3 0 ; 0 -1/L ; 0 0];

%Matrix B for phase 2
B = [0 0 ; 0 -1/L ; 0 0];
%C and D stay the same
C=[1 0 0; 0 1 0; 0 0 1];
D=[0 0 ; 0 0 ; 0 0];


%Phase 1 discrete system
sys1_1 = ss(A1,B1,C,D);
sys1_1_dis=c2d(sys1_1,10^(-6));

A11_dis = sys1_1_dis.A;
B11_dis = sys1_1_dis.B;
C11_dis = sys1_1_dis.C;
D11_dis = sys1_1_dis.D;

%Phase 2 discrete system
sys1_2 = ss(A1,B,C,D);
sys1_2_dis=c2d(sys1_2,10^(-6));

A12_dis = sys1_2_dis.A;
B12_dis = sys1_2_dis.B;
C12_dis = sys1_2_dis.C;
D12_dis = sys1_2_dis.D;

%Initializing everything
V = 90*ones(10000,1);
Vs = 300*ones(10000,1);
Vd = -300*ones(10000,1);
Id = 0*ones(10000,1);
VL = 0*ones(10000,1);
Vc(1)=0;
If(1)=0;
Il(1)=0;
t = 0;

%for-loop
for p = 1:1:10000
   
    if((mod(t,1/Fs))<=1/(Fs*2)) 
        %Phase 1 calculation
        Il(p+1) = A11_dis(1)*Il(p)+ A11_dis(4)*If(p) + A11_dis(7)*Vc(p) + B11_dis(1)*Vs(p) + B11_dis(4)*V(p);
        If(p+1) = A11_dis(2)*Il(p)+ A11_dis(5)*If(p) + A11_dis(8)*Vc(p) + B11_dis(2)*Vs(p) + B11_dis(5)*V(p);
        Vc(p+1) = A11_dis(3)*Il(p)+ A11_dis(6)*If(p) + A11_dis(9)*Vc(p) + B11_dis(3)*Vs(p) + B11_dis(6)*V(p);
        Vd(p+1) = -300;
        Id(p+1) = 0;
        VL(p+1) = Vs(p+1) - Vc(p+1);
    else
        %Phase 2 calculation
        Il(p+1) = A12_dis(1)*Il(p)+ A12_dis(4)*If(p) + A12_dis(7)*Vc(p) + B12_dis(1)*Vs(p) + B12_dis(4)*V(p);
        If(p+1) = A12_dis(2)*Il(p)+ A12_dis(5)*If(p) + A12_dis(8)*Vc(p) + B12_dis(2)*Vs(p) + B12_dis(5)*V(p);
        Vc(p+1) = A12_dis(3)*Il(p)+ A12_dis(6)*If(p) + A12_dis(9)*Vc(p) + B12_dis(3)*Vs(p) + B12_dis(6)*V(p);
        Vd(p+1) = 0;
        Id(p+1) =Il(p+1);
        VL(p+1) = -Vc(p+1); 
    end;
    t= t+ 1/(100*Fs);
 end;

t = 0:1/(100*Fs):10000/(100*Fs);

%Displaying requested plots
figure;
subplot(2,1,1);
plot(t,Il);
title('Coil current for DT = 0.5, L = 0.001');
axis([0 0.01 0 25]);
xlabel('Time(sec)');
ylabel('IL(A)');

subplot(2,1,2);
plot(t,VL);
title('Coil voltage for DT = 0.5, L = 0.001');
axis([0 0.01 -200 300]);
xlabel('Time(sec)');
ylabel('VL(V)');

figure;
subplot(2,1,1);
plot(t,Il-If);
title('Capacitor current for DT = 0.5, L = 0.001');
axis([0 0.01 -15 20]);
xlabel('Time(sec)');
ylabel('Ic(A)');

subplot(2,1,2);
plot(t,Vc);
title('Capacitor voltage for DT = 0.5, L = 0.001');
axis([0 0.01 0 200]);
xlabel('Time(sec)');
ylabel('Vc(V)');

figure;
plot(t,If);
title('Load current for DT=0.5, L = 0.001');
axis([0 0.01 -10 25]);
xlabel('Time(sec)');
ylabel('ILoad(A)');

%% DT = 0.5, L = 0.01 

%State space model
%Phase 1 discrete system
sys2_1 = ss(A2,B2,C,D);
sys2_1_dis=c2d(sys2_1,10^(-6));

A21_dis = sys2_1_dis.A;
B21_dis = sys2_1_dis.B;
C21_dis = sys2_1_dis.C;
D21_dis = sys2_1_dis.D;

%Phase 2 discrete system
sys2_2 = ss(A2,B,C,D);
sys2_2_dis=c2d(sys2_2,10^(-6));

A22_dis = sys2_2_dis.A;
B22_dis = sys2_2_dis.B;
C22_dis = sys2_2_dis.C;
D22_dis = sys2_2_dis.D;

%Initializing
V = 90*ones(20000,1);
Vs = 300*ones(20000,1);
Vd = -300*ones(20000,1);
Id = 0*ones(20000,1);
VL = 0*ones(20000,1);
Vc(1)=0;
If(1)=0;
Il(1)=0;
t = 0;

for p = 1:1:20000   
    if((mod(t,1/Fs))<=1/(Fs*2))
        %Phase 1 calculations
        Il(p+1) = A21_dis(1)*Il(p)+ A21_dis(4)*If(p) + A21_dis(7)*Vc(p) + B21_dis(1)*Vs(p) + B21_dis(4)*V(p);
        If(p+1) = A21_dis(2)*Il(p)+ A21_dis(5)*If(p) + A21_dis(8)*Vc(p) + B21_dis(2)*Vs(p) + B21_dis(5)*V(p);
        Vc(p+1) = A21_dis(3)*Il(p)+ A21_dis(6)*If(p) + A21_dis(9)*Vc(p) + B21_dis(3)*Vs(p) + B21_dis(6)*V(p);
        Vd(p+1) = -300;
        Id(p+1) = 0;
        VL(p+1) = Vs(p+1) - Vc(p+1); 
    else
        %Phase 2 calculations
        Il(p+1) = A22_dis(1)*Il(p)+ A22_dis(4)*If(p) + A22_dis(7)*Vc(p) + B22_dis(1)*Vs(p) + B22_dis(4)*V(p);
        If(p+1) = A22_dis(2)*Il(p)+ A22_dis(5)*If(p) + A22_dis(8)*Vc(p) + B22_dis(2)*Vs(p) + B22_dis(5)*V(p);
        Vc(p+1) = A22_dis(3)*Il(p)+ A22_dis(6)*If(p) + A22_dis(9)*Vc(p) + B22_dis(3)*Vs(p) + B22_dis(6)*V(p);
        Vd(p+1) = 0;
        Id(p+1) =Il(p+1);
        VL(p+1) = -Vc(p+1);
    end;
   t= t+ 1/(100*Fs);
 end;

t = 0:1/(100*Fs):20000/(100*Fs);

%Displaying requested plots
figure;
subplot(2,1,1);
plot(t,Il);
title('Coil current for DT = 0.5, L = 0.01');
axis([0 0.01 0 20]);
xlabel('Time(sec)');
ylabel('IL(A)');

subplot(2,1,2);
plot(t,VL);
title('Coil voltage for DT = 0.5, L = 0.01');
axis([0 0.01 -200 300]);
xlabel('Time(sec)');
ylabel('VL(V)');

figure;
subplot(2,1,1);
plot(t,Il-If);
title('Capacitor current for DT = 0.5, L = 0.01');
axis([0 0.01 -5 10]);
xlabel('Time(sec)');
ylabel('Ic(A)');

subplot(2,1,2);
plot(t,Vc);
title('Capacitor voltage for DT=0.5, L = 0.01');
axis([0 0.01 0 150]);
xlabel('Time(sec)');
ylabel('Vc(V)');

figure;
plot(t,If);
title('Load current for DT = 0.5, L = 0.01');
axis([0 0.01 -10 20]);
xlabel('Time(sec)');
ylabel('ILoad(A)');

%% DT = 0.8, L = 0.001

%State space model
%Phase 1 discrete system
sys3_1 = ss(A3,B3,C,D);
sys3_1_dis=c2d(sys3_1,10^(-6));

A31_dis = sys3_1_dis.A;
B31_dis = sys3_1_dis.B;
C31_dis = sys3_1_dis.C;
D31_dis = sys3_1_dis.D;

%Phase 2 discrete system
sys3_2 = ss(A3,B,C,D);
sys3_2_dis=c2d(sys3_2,10^(-6));

A32_dis = sys3_2_dis.A;
B32_dis = sys3_2_dis.B;
C32_dis = sys3_2_dis.C;
D32_dis = sys3_2_dis.D;

%Initializing
V = 90*ones(20000,1);
Vs = 300*ones(20000,1);
Vd = -300*ones(20000,1);
Id = 0*ones(20000,1);
VL = 0*ones(20000,1);
Vc(1)=0;
If(1)=0;
Il(1)=0;
t = 0;

for p = 1:1:20000
   
    if((mod(t,1/Fs))<=(1/Fs)*0.8)
        %Phase 1 calculations
        Il(p+1) = A31_dis(1)*Il(p)+ A31_dis(4)*If(p) + A31_dis(7)*Vc(p) + B31_dis(1)*Vs(p) + B31_dis(4)*V(p);
        If(p+1) = A31_dis(2)*Il(p)+ A31_dis(5)*If(p) + A31_dis(8)*Vc(p) + B31_dis(2)*Vs(p) + B31_dis(5)*V(p);
        Vc(p+1) = A31_dis(3)*Il(p)+ A31_dis(6)*If(p) + A31_dis(9)*Vc(p) + B31_dis(3)*Vs(p) + B31_dis(6)*V(p);
        Vd(p+1) = -300;
        Id(p+1) = 0;
        VL(p+1) = Vs(p+1) - Vc(p+1);
        
    else
        %Phase 2 calculations
        Il(p+1) = A32_dis(1)*Il(p)+ A32_dis(4)*If(p) + A32_dis(7)*Vc(p) + B32_dis(1)*Vs(p) + B32_dis(4)*V(p);
        If(p+1) = A32_dis(2)*Il(p)+ A32_dis(5)*If(p) + A32_dis(8)*Vc(p) + B32_dis(2)*Vs(p) + B32_dis(5)*V(p);
        Vc(p+1) = A32_dis(3)*Il(p)+ A32_dis(6)*If(p) + A32_dis(9)*Vc(p) + B32_dis(3)*Vs(p) + B32_dis(6)*V(p);
        Vd(p+1) = 0;
        Id(p+1) = Il(p+1);
        VL(p+1) = -Vc(p+1);
    end;
   t= t+ 1/(100*Fs);
 end;

t = 0:1/(100*Fs):20000/(100*Fs);

%Displaying requested plots
figure;
subplot(2,1,1);
plot(t,Il);
title('Coil current for DT = 0.8, L = 0.001');
axis([0 0.01 0 60]);
xlabel('Time(sec)');
ylabel('IL(A)');

subplot(2,1,2);
plot(t,VL);
title('Coil voltage for DT = 0.8, L = 0.001');
axis([0 0.01 -300 350]);
xlabel('Time(sec)');
ylabel('VL(V)');

figure;
subplot(2,1,1);
plot(t,Il-If);
title('Capacitor current for DT = 0.8, L = 0.001');
axis([0 0.01 -20 25]);
xlabel('Time(sec)');
ylabel('Ic(A)');

subplot(2,1,2);
plot(t,Vc);
title('Capacitor voltage for DT = 0.8, L = 0.001');
axis([0 0.01 0 300]);
xlabel('Time(sec)');
ylabel('Vc(V)');

figure;
plot(t,If);
title('Load current for DT = 0.8, L = 0.001');
axis([0 0.01 -10 65]);
xlabel('Time(sec)');
ylabel('ILoad(A)');

