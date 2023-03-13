%% 4.1
clear all;
close all;
clc;

Vs = 300;                       % Source voltage
C = 0.00001;                    % Capacitance
R = 3;                          % Resistance
L = 0.001;                      % Induction
fs = 10000;                     % Frequency
Ts = 1/fs;                      % Period


dt = 0.000001;                  % Step
n = 1000;                       % Amount of periods
t = 0:dt:n*Ts-dt;               % Time

Kp = 0.001;
Ki = 0.5;

%Motor parameters
Ra = 3;             % Resistance
La = 0.0005;        % Induction
J = 0.005;          % kgm^2
Kt = 0.3;           % Nm/A
Ke = 0.3;           % v/rad/s
TL = 15;            % Torque
omega = 80;         % Reference value

pulse = (sawtooth(2*pi*fs*t,(1/2))+1)/2;

% Phase 1 state space
A1 = [0 0 -1/L; 0 -Ra/La 1/La; 1/C -1/C 0];
B1 = [1/L 0; 0 -1/La; 0 0];
C1 = eye(3);
D1 = zeros(3,2);

sys1 = ss(A1,B1,C1,D1);
sys_dis1 = c2d(sys1,dt);

sysa = sys_dis1.A;
sysb = sys_dis1.B;
sysc = sys_dis1.C;
sysd = sys_dis1.D;

% Phase 2 state space
A2 = [0 0 -1/L; 0 -Ra/La 1/La; 1/C -1/C 0];
B2 = [0 0; 0 -1/La; 0 0];
C2 = eye(3);
D2 = zeros(3,2);


sys2 = ss(A2,B2,C2,D2);
sys_dis2 = c2d(sys2,dt);

sys2a = sys_dis2.A;
sys2b = sys_dis2.B;
sys2c = sys_dis2.C;
sys2d = sys_dis2.D;

% Motor state space
A_motor = Kt/J;
B_motor = -1/J;
C_motor = Ke;
D_motor = 0;

sysmotor = ss(A_motor,B_motor,C_motor,D_motor);
sys_dismotor = c2d(sysmotor,dt);

motora = sys_dismotor.A;
motorb = sys_dismotor.B;
motorc = sys_dismotor.C;
motord = sys_dismotor.D;

% Controller state space
Ac = 0;
Bc = 1;
Cc = Ki;   
Dc = Kp;

sysc = ss(Ac,Bc,Cc,Dc);
sysc_dis = c2d(sysc,dt);

cona=sysc_dis.A;
conb=sysc_dis.B;
conc=sysc_dis.C;
cond=sysc_dis.D;

x1 = 0;  
x2 = 0; 
g = 100; 
v = 0;   
j = 1;

% Helping vector
IL = zeros(size(t));                        
VL = zeros(size(t));  
Va = zeros(size(t));
Ia = zeros(size(t));                        
Vc = zeros(size(t));                        
Ic = zeros(size(t));                                                
Te = zeros(size(t));                        
wm = zeros(size(t));                          
error = zeros(size(t));                     
x = zeros(size(t));                          
duty = zeros(size(n));                      
step_p = zeros(size(g));                    


for k=0:1:n-1
    if (j==1)
        % Generating pulse
        for i=1:g
            if(pulse(i)>= duty(j))
                step_p(i) = min(pulse);
            else
                step_p(i) = max(pulse);
            end
        end
    else
        for i=1:g
            if(pulse(i)>= duty(j-1))
                step_p(i) = min(pulse);
            else
                step_p(i) = max(pulse);
            end
        end
    end
    
    % Calculating current and voltages
    for u = 1:Ts/dt
        if(step_p(u) == max(step_p))
            % phase 1
            IL(u+1+v) = sysa(1,:)*[IL(u+v); Ia(u+v); Vc(u+v)] + sysb(1,:)*[Vs; Va(u+v)];
            Ia(u+1+v) = sysa(2,:)*[IL(u+v); Ia(u+v); Vc(u+v)] + sysb(2,:)*[Vs; Va(u+v)];
            Ic(u+v) = IL(u+v) - Ia(u+v);
            Vc(u+1+v) = sysa(3,:)*[IL(u+v); Ia(u+v); Vc(u+v)] + sysb(3,:)*[Vs; Va(u+v)];
            VL(u+v) = Vs - Vc(u+v);
            
        else
            % phase 2
            IL(u+1+v) = sys2a(1,:)*[IL(u+v); Ia(u+v); Vc(u+v)] + sys2b(1,:)*[Vs; Va(u+v)];
            Ia(u+1+v) = sys2a(2,:)*[IL(u+v); Ia(u+v); Vc(u+v)] + sys2b(2,:)*[Vs; Va(u+v)];
            Ic(u+v) = IL(u+v) - Ia(u+v);
            Vc(u+1+v) = sys2a(3,:)*[IL(u+v); Ia(u+v); Vc(u+v)] + sys2b(3,:)*[Vs; Va(u+v)];
            VL(u+v) = - Vc(u+v);
            
        end
        
        %omega computing
        Te(u+v) = Kt*Ia(u+v);
        wm(u+1+v) = motora*Ia(u+v) + motorb*TL;
        Va(u+1+v) = motorc*wm(u+1+v);
        
        % error computing
        error(u+v) =  wm(u+v) - omega;
        x(u+1+v) = cona*x(u+v) + conb*error(u+v);
        
        if(u == g)
            x1 = x(u+v);
            x2 = error(u+v);
        end
        
    end
    
    duty(j) = abs(conc*x1 + cond*x2);
    
    %incrementing constants
    v = v+g;
    j = j+1;
end


% Plots
figure;
plot(t,IL(1:end-1));
grid on;
xlabel('Sec');
ylabel('A');
title('Inductor current T_L = 15 Nm');

figure;
plot(t,VL);
grid on;
xlabel('Sec');
ylabel('V');
title('Inductor voltage T_L = 15 Nm');

figure;
plot(t,Ia(1:end-1));
grid on;
xlabel('Time(sec)');
ylabel('I_{Load}(A)');
title('Load current T_L = 15 Nm');

figure;
plot(t,Vc(1:end-1));
grid on;
xlabel('Sec)');
ylabel('V');
title('Capacitor voltage T_L = 15 Nm');

figure;
plot(t,Ic);
grid on;
xlabel('Sec');
ylabel('A');
title('Capacitor current T_L = 15 Nm');

figure;
plot(t,wm(1:end-1));
grid on;
xlabel('Sec)');
ylabel('rad/s');
title('Angular frequency T_L = 15 Nm');

figure;
plot(t,Te);
grid on;
xlabel('Sec)');
ylabel('Nm');
title('Electric torque T_L = 15 Nm');
