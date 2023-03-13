clear all;
close all;

% ---------------- 3 -----------------

V = 100;                 % Input voltage
amp = 1;                 % Amplitude
R = 10;                  % Resistance
L = 0.025;               % Induction
f = 50;                  % Frequency
om = 2*pi*f;             % Omega
dt = 2*(10^-5);          % Step
n = 2;              
t = 0:dt:n*(1/f);        % Time vector
om_t = om*t;             % Calculation of w*t
w_t = mod(om*t,2*pi);    % on [0,2pi]

% Initialization 
p1 = zeros(size(t));
p2 = zeros(size(t));
p3 = zeros(size(t));
p4 = zeros(size(t));
p5 = zeros(size(t));
p6 = zeros(size(t));



% Pulse creation
for k=1:size(t,2)
    
    if (w_t(k)>=0 && w_t(k)<pi)  
        p1(k) = amp;
        p4(k) = 0;
    elseif (w_t(k)>=pi && w_t(k)<=2*pi) 
        p4(k) = amp;
        p1(k) = 0;
    end
    
    if (w_t(k)>=0 && w_t(k)<2*pi/3) 
        p3(k) = 0;
        p6(k) = amp;
    elseif (w_t(k)>=2*pi/3 && w_t(k)<5*pi/3) 
        p3(k) = amp;
        p6(k) = 0;
    elseif (w_t(k)>=5*pi/3 && w_t(k)<=2*pi) 
        p3(k) = 0;
        p6(k) = amp;
    end
    
    if (w_t(k)>=0 && w_t(k)<pi/3)     
        p5(k) = amp;
        p2(k) = 0;
    elseif (w_t(k)>=pi/3 && w_t(k)<4*pi/3)    
        p5(k) = 0;
        p2(k) = amp;
    elseif (w_t(k)>=4*pi/3 && w_t(k)<=2*pi)   
        p5(k) = amp;
        p2(k) = 0;
    end
    
end


% Initialization
Vao = zeros(size(t));
Vbo = zeros(size(t));
Vco = zeros(size(t));
Vab = zeros(size(t));
Vbc = zeros(size(t));
Vca = zeros(size(t));
Van = zeros(size(t));
Vbn = zeros(size(t));
Vcn = zeros(size(t));

% Voltages calculations
for k=1:size(t,2)
    
    if (p1(k) == amp)
        Vao(k) = V;
    elseif ( p4(k) == amp)
        Vao(k) = 0;
    end
    
    if (p3(k) == amp)
        Vbo(k) = V;
    elseif ( p6(k) == amp)
        Vbo(k) = 0;
    end
    
    if (p5(k) == amp)
        Vco(k) = V;
    elseif (p2(k) == amp)
        Vco(k) = 0;
    end
    
end

figure()
plot(om_t,Vao)
title('VA0')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -20 120])

figure()
plot(om_t,Vbo)
title('VB0')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -20 120])

figure()
plot(om_t,Vco)
title('VC0')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -20 120])

% Polar Voltages
Vab=Vao-Vbo;
Vca=Vco-Vao;
Vbc=Vbo-Vco;

% Phase Voltages
Van=(Vab-Vca)/3;
Vbn=(Vbc-Vab)/3;
Vcn=(Vca-Vbc)/3;

% Polar Voltages
figure()
plot(om_t,Vab)
title('Polar VAB')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -120 120])

figure()
plot(om_t,Vbc)
title('Polar VBC')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -120 120])

figure()
plot(om_t,Vca)
title('Polar VCA')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -120 120])

%Phase Voltages
figure()
plot(om_t,Van)
title('Phase VAN')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -100 100])

figure()
plot(om_t,Vbn)
title('Phase VBN')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -100 100])

figure()
plot(om_t,Vcn)
title('Phase VCN')
xlabel('ω(rad/s)')
ylabel('Voltage(V)');
axis([0 4*pi -100 100])


% Initialization
IAn = zeros(size(t));
IBn = zeros(size(t));
ICn = zeros(size(t));

% State space parametres
A_3 = -R/L;
B_3 = 1/L;
C_3 = 1;
D_3 = 0;

% State space model
sys_3 = ss(A_3,B_3,C_3,D_3);
% Discrete system
sys_3_discrete = c2d(sys_3, dt);
sys3A = sys_3_discrete.A;
sys3B = sys_3_discrete.B;
sys3C = sys_3_discrete.C;
sys3D = sys_3_discrete.D;

%Load's current
for k=1:size(t,2)
    IAn(k+1) = sys3A*IAn(k) + sys3B*Van(k);
    IBn(k+1) = sys3A*IBn(k) + sys3B*Vbn(k);
    ICn(k+1) = sys3A*ICn(k) + sys3B*Vcn(k);
end

figure()
plot(om_t,IAn(1:end-1))
title('Phase IAN')
xlabel('ω(rad/s)')
ylabel('I(A)');
axis([0 4*pi -15 15])

figure()
plot(om_t,IBn(1:end-1))
title('Phase IBN')
xlabel('ω(rad/s)')
ylabel('I(A)');
axis([0 4*pi -15 15])

figure()
plot(om_t,ICn(1:end-1))
xlabel('ω(rad/s)')
title('Phase ICN')
ylabel('I(A)');
axis([0 4*pi -15 15])


% Power calculation
PoutA = (IAn(1:end-1).^2)*R;
PoutB = (IBn(1:end-1).^2)*R;
PoutC = (ICn(1:end-1).^2)*R;
figure()
plot(t,PoutC)
title('Output Power');
xlabel('t(sec)');
ylabel('P(W)');


NFFT = 2000;
p3ph = (Vcn.*ICn(1:end-1));                   % Power
P3ph = mean(p3ph((n-1)*1000+1:end));          % Active power
V3rms = sqrt(sum(Vcn(1:NFFT).^2)/(NFFT/2));   % Vrms phase A
I3rms = sqrt(sum(ICn(1:NFFT).^2)/(NFFT/2));   % Irms phase A
S3ph = V3rms * I3rms;                         % Apparent power
Pf3 = P3ph/S3ph;                              % Power factor
disp('Power Factor in each phase:');
disp(Pf3);

% Fourier transformation
F = (1/dt)/2*linspace(0,1,(NFFT/2)+1);   
VCnF = fft(Vcn,NFFT);               
figure();
stem(F,abs(VCnF(1:(NFFT/2)+1)));
title('Voltage Harmonic')
ylabel('Voltage(V)')
xlabel('Frequency(Hz)')
