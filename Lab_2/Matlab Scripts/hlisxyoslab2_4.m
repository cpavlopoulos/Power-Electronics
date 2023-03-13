% 4
clear all;
close all;

Vrms = 230;         % initial Vrms
f = 50;             % freq
dt = 0.0002;         % step
t = 0:dt:10;        % time vector
t1 = 0:dt:5-dt;
t2 = 5:dt:10-dt;        

Vm = Vrms*sqrt(2)*sqrt(3);  

%creating Va,Vb,Vc
Vab =Vm*sin(2*pi*f*t);
Vbc =Vm*sin(2*pi*f*t-(2*pi/3));
Vca =Vm*sin(2*pi*f*t-(4*pi/3));

Vba = -Vm*sin(2*pi*f*t);
Vcb = -Vm*sin(2*pi*f*t-(2*pi/3));
Vac = -Vm*sin(2*pi*f*t-(4*pi/3));


Iref = 75;  % reference current
R= 2.5;     % ohm
L= 0.04;    % henry

% state space parametres
A = -R/L;
B = 1/L;
C = 1;
D = 0;

sys_3 = ss(A,B,C,D);
sys_d3 = c2d(sys_3, dt);

%Ki=7 Kp=0.8
Ac = 0;
Bc = 1;
Cc = 7;
Dc = 0.8;

sys_c = ss(Ac,Bc,Cc,Dc);
sys_dc = c2d(sys_c, dt);


I_L1 = zeros(1,size(t,2)+1);
V_L = zeros(1,size(t,2));
x = zeros(1,size(t,2));
ang = zeros(1,size(t,2));
i=1;

for j = t1-dt
     
    if i==1
        if (mod(2*pi*f*t1(i)*180/pi, 2*180) < ang(i)) && (mod(2*pi*f*t1(i)*180/pi_t, 2*180) > 0)  
            V_L(i) = Vca(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= ang(i)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<60+ ang(i))  
            V_L(i) = Vcb(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 60 + ang(i)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<120 + ang(i))   
            V_L(i) = Vab(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 120 + ang(i)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<180 + ang(i))  
            V_L(i) = Vac(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 180 + ang(i)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<240 + ang(i))  
            V_L(i) = Vbc(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 240 +ang(i)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<300 + ang(i))   
            V_L(i) = Vba(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 300 + ang(i)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<=360)  
            V_L(i) = Vca(i);
        end
        
    else
        
        if (mod(2*pi*f*t1(i)*180/pi, 2*180) < ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180) > 0)    
            V_L(i) = Vca(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<60+ ang(i-1))  
            V_L(i) = Vcb(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 60 + ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<120 + ang(i-1))   
            V_L(i) = Vab(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 120 + ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<180 + ang(i-1))  
            V_L(i) = Vac(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 180 + ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<240 + ang(i-1))  
            V_L(i) = Vbc(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 240 +ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<300 + ang(i-1))   
            V_L(i) = Vba(i);
        elseif (mod(2*pi*f*t1(i)*180/pi, 2*180) >= 300 + ang(i-1)) && (mod(2*pi*f*t1(i)*180/pi, 2*180)<=360)    
            V_L(i) = Vca(i);
        end
        
    end
    
    %Load current
    I_L1(i+1) = sys_d3.A*I_L1(i)+sys_d3.B*V_L(i);
     
    %Model x,y
    x(i+1) = sys_dc.A*x(i) + sys_dc.B*(I_L1(i) - Iref);
    ang(i) = sys_dc.C*x(i) + sys_dc.D*(I_L1(i) - Iref);
   
    i=i+1;
    
end


% After 5 seconds we decrease the voltage to 200V
Vrms_new = 200;         % new Vrms
Vm = Vrms_new*sqrt(2);  % Final voltage amplitude

%creating Va,Vb,Vc
Vab =Vm*sin(2*pi*f*t);
Vbc =Vm*sin(2*pi*f*t-(2*pi/3));
Vca =Vm*sin(2*pi*f*t-(4*pi/3));

Vba = -Vm*sin(2*pi*f*t);
Vcb = -Vm*sin(2*pi*f*t-(2*pi/3));
Vac = -Vm*sin(2*pi*f*t-(4*pi/3));

for j=t2-dt
    
    if (mod(2*pi*f*t(i)*180/pi, 2*180) < ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180) > 0)    
        V_L(i) = Vca(i);
    elseif (mod(2*pi*f*t(i)*180/pi, 2*180) >= ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180)<60+ ang(i-1))  
        V_L(i) = Vcb(i);
    elseif (mod(2*pi*f*t(i)*180/pi, 2*180) >= 60 + ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180)<120 + ang(i-1))   
        V_L(i) = Vab(i);
    elseif (mod(2*pi*f*t(i)*180/pi, 2*180) >= 120 + ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180)<180 + ang(i-1))  
        V_L(i) = Vac(i);
    elseif (mod(2*pi*f*t(i)*180/pi, 2*180) >= 180 + ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180)<240 + ang(i-1))  
        V_L(i) = Vbc(i);
    elseif (mod(2*pi*f*t(i)*180/pi, 2*180) >= 240 +ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180)<300 + ang(i-1))   
        V_L(i) = Vba(i);
    elseif (mod(2*pi*f*t(i)*180/pi, 2*180) >= 300 + ang(i-1)) && (mod(2*pi*f*t(i)*180/pi, 2*180)<=360)    
        V_L(i) = Vca(i);
    end
    
    % current in load
    I_L1(i+1) = sys_d3.A*I_L1(i)+sys_d3.B*V_L(i);
    
    %angle calculation
    x(i+1) = sys_dc.A*x(i) + sys_dc.B*(I_L1(i) - Iref);
    ang(i) = sys_dc.C*x(i) + sys_dc.D*(I_L1(i) - Iref);
    
    i=i+1;
    
end


figure;
plot(t,I_L1(1:end-1));
hold on
plot(t,75*ones(1,size(t,2)),'LineWidth',1);
hold off
xlabel('time(sec)');
ylabel('I_L(A)');
title('Load Current');
legend('I_Load','75 A');

figure;
plot(t,ang);
xlabel('time(sec)');
title('Angle');

figure;
plot(t,V_L);
grid on;
xlabel('time(sec)');
title('Load Voltage');

