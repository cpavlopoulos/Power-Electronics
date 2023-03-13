close all;
clear all;

%% ---------- 2 -------------

V = 100;            % Input voltage
R = 10;             % Resistance
L = 0.025;          % Induction
f = 50;             % Frequency
om = 2*pi*f;        % Omega
dt = 10^-6;         % Step
n =1;               % Number of periods
t = 0:dt:n*1/f;     % Vector for time
ma = 0.9;           % Ma is 0.9
mf1 = 40;           % Mf is 40 for 1st case
mf2 = 200;          % Mf is 200 for 2nd case
sine = sin(om*t)*V; % Reference signal
amp = 1/2;          % Triangle wave's width
Vm_c = V/ma;        % Carrier Vm


%Initialization
Va  = zeros(size(t));
Vb  = zeros(size(t));
Vab = zeros(size(t));
IQ1 = zeros(size(t));
IQ2 = zeros(size(t));
IQ3 = zeros(size(t));
IQ4 = zeros(size(t));
IL = zeros(size(t));

% System's parameters
A1 = -R/L;
B1 = 1/L;
C1 = 1;
D1 = 0;

% State space model
sys1 = ss(A1,B1,C1,D1);
% Discrete system
sys1_dis = c2d(sys1, dt);
dis1 = sys1_dis.A;
dis2 = sys1_dis.B;
dis3 = sys1_dis.C;
dis4 = sys1_dis.D;

%Calculations for Mf = 40
f_carrier = f*mf1;  %Carrier signal
triangle = ((sawtooth(2*pi*f_carrier*t, amp)))*Vm_c;
k = 1;  
for i = t
    if(sine(k) > triangle(k))
        Va(k) = V;
        IQ1(k) = IL(k);
    elseif(sine(k) < triangle(k))
        IQ4(k) = IL(k);
    end
        
    if(-sine(k) > triangle(k))
        Vb(k) = V;
        IQ3(k) = IL(k);
    elseif(-sine(k) < triangle(k))
        IQ2(k) = IL(k);
    end
    
    Vab(k) = Va(k) - Vb(k);
    IL(k+1) = dis1*IL(k) + dis2*Vab(k);
        
    k=k+1;
end
    
figure;
plot(t,triangle);
hold on;
plot(t,sine);
plot(t,-sine);
legend('Carrier','Sine','-Sine');
title('PWM for Mf = 40');
xlabel('t(sec)');
hold off;
   
figure;
subplot(1,2,1);
plot(t,Va);
title('Va for Mf = 40');
xlabel('t(sec)');
ylabel('Va(V)');
axis([0 0.02 -20 120]);
 

subplot(1,2,2);
plot(t,Vb);
title('Vb for Mf = 40');
xlabel('t(sec)');
ylabel('Vb(V)');
axis([0 0.02 -20 120]);
      
figure;
plot(t,Vab);
title('Output voltage for Mf = 40');
xlabel('t(sec)');
ylabel('Vab(V)');
axis([0 0.02 -120 120]);
  
figure;
plot(t,IL(1:end-1));
title('Load Current for Mf = 40');
xlabel('t(sec)');
ylabel('IL(A)');
   
figure;
subplot(2,2,1);
plot(t,IQ1);
title('Q1 for Mf = 40');
xlabel('t(sec)');
ylabel('IQ1(A)');
  
subplot(2,2,2);
plot(t,IQ2);
title('Q2 for Mf = 40');
xlabel('t(sec)');
ylabel('IQ2(A)');
  
subplot(2,2,3);
plot(t,IQ3);
title('Q3 for Mf = 40');
xlabel('t(sec)');
ylabel('IQ3(A)');
  
subplot(2,2,4);
plot(t,IQ4);
title('Q4 for Mf = 40');
xlabel('t(sec)');
ylabel('IQ4(A)');
  
% Calculation of output power
Pout = (IL(1:end-1).^2)*R;    
figure()
plot(t,Pout)
title('Output Power for Mf = 40 using PWM');
grid on;
xlabel('t(sec)');
ylabel('Pout(W)');
    

NFFT = 10000;
p1ph = (Vab.*IL(1:end-1));                          % Power
P1ph = mean(p1ph((n-1)*1000+1:end));                % Active power
Vab_rms = sqrt(sum(Vab(1:NFFT).^2)/(NFFT/2));       % Vrms of output    
Iab_rms = sqrt(sum(IL(1:NFFT).^2)/(NFFT/2));        % Irms 
S = Vab_rms * Iab_rms;                              % Apparent power
Pf = P1ph/S;                                        % Power factor
fprintf('Power Factor for Mf = 40 is ');            % Displaying Pf 
disp(Pf);  
        
% Fourier transformation
F = (1/dt)/2*linspace(0,1,(NFFT/2)+1);
Vab_F = fft(Vab,NFFT);
figure();
stem(F,(abs(Vab_F(1:(NFFT/2)+1)/Vab_F(2))));
title('Harmonics for Mf = 40');
xlabel('Frequency(Hz)')
ylabel('Vab')
       

% Calculations for Mf = 200
f_carrier = f*mf2;  %Carrier signal
triangle = ((sawtooth(2*pi*f_carrier*t, amp)))*Vm_c;
k = 1;  
for i = t
    if(sine(k) > triangle(k))
        Va(k) = V;
        IQ1(k) = IL(k);
    elseif(sine(k) < triangle(k))
        IQ4(k) = IL(k);
    end
        
    if(-sine(k) > triangle(k))
        Vb(k) = V;
        IQ3(k) = IL(k);
    elseif(-sine(k) < triangle(k))
        IQ2(k) = IL(k);
    end
    
    Vab(k) = Va(k) - Vb(k);
    IL(k+1) = dis1*IL(k) + dis2*Vab(k);
        
    k=k+1;
end
    
figure;
plot(t,triangle);
hold on;
plot(t,sine);
plot(t,-sine);
legend('Carrier','Sine','-Sine');
title('PWM for Mf = 200');
xlabel('t(sec)');
hold off;
    
figure;
subplot(1,2,1);
plot(t,Va);
title('Va for Mf = 200');
xlabel('t(sec)');
ylabel('Va(V)');
axis([0 0.02 -20 120]);

subplot(1,2,2);
plot(t,Vb);
title('Vb for Mf = 200');
xlabel('t(sec)');
ylabel('Vb(V)');
axis([0 0.02 -20 120]);
      
figure;
plot(t,Vab);
title('Output voltage for Mf = 200');
xlabel('t(sec)');
ylabel('Vab(V)');
axis([0 0.02 -120 120]);
  
figure;
plot(t,IL(1:end-1));
title('Load Current for Mf = 200');
xlabel('t(sec)');
ylabel('IL(A)');
   
figure;
subplot(2,2,1);
plot(t,IQ1);
title('Q1 for Mf = 200');
xlabel('t(sec)');
ylabel('IQ1(A)');
  
subplot(2,2,2);
plot(t,IQ2);
title('Q2 for Mf = 200');
xlabel('t(sec)');
ylabel('IQ2(A)');

subplot(2,2,3);
plot(t,IQ3);
title('Q3 for Mf = 200');
xlabel('t(sec)');
ylabel('IQ3(A)');
  
subplot(2,2,4);
plot(t,IQ4);
title('Q4 for Mf = 200');
xlabel('t(sec)');
ylabel('IQ4(A)');
  
% Calculation of output power
Pout = (IL(1:end-1).^2)*R;    
figure()
plot(t,Pout)
title('Output Power for Mf = 200 using PWM');
grid on;
xlabel('t(sec)');
ylabel('Pout(W)');
    
    
NFFT = 10000;
p1ph = (Vab.*IL(1:end-1));                          % Power
P1ph = mean(p1ph((n-1)*1000+1:end));                % Active power
Vab_rms = sqrt(sum(Vab(1:NFFT).^2)/(NFFT/2));       % Vrms of output    
Iab_rms = sqrt(sum(IL(1:NFFT).^2)/(NFFT/2));        % Irms 
S = Vab_rms * Iab_rms;                              % Apparent power
Pf = P1ph/S;                                        % Power factor
fprintf('Power Factor for Mf = 200 is ');           % Displaying Pf 
disp(Pf);  
        
% Fourier transformation
F = (1/dt)/2*linspace(0,1,(NFFT/2)+1);
Vab_F = fft(Vab,NFFT);
figure();
stem(F,(abs(Vab_F(1:(NFFT/2)+1)/Vab_F(2))));
title('Harmonics for Mf = 200');
xlabel('Frequency(Hz)')
ylabel('Vab')
