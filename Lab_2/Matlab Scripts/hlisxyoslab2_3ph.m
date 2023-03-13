clear all;
close all;

%variables used for creating Va,Vb,Vc
f = 50;
Vm = 230*sqrt(2)*sqrt(3);
i=1;
dt=2*10^(-5); %sec
T=10;   %sec
t=0:dt:T-dt;  %sampling time

%creating Va,Vb,Vc
Vab =Vm*sin(2*pi*f*t);
Vbc =Vm*sin(2*pi*f*t-(2*pi/3));
Vca =Vm*sin(2*pi*f*t-(4*pi/3));

Vba = -Vm*sin(2*pi*f*t);
Vcb = -Vm*sin(2*pi*f*t-(2*pi/3));
Vac = -Vm*sin(2*pi*f*t-(4*pi/3));

t1 = 0:dt:0.02; %first period
w = 0:360/1000:(T/0.02)*360-(360/1000); % degrees 

figure
plot(w,Vab,w,Vbc,'r',w,Vca,'g',w,Vba,'b--',w,Vcb,'r--',w,Vac,'g--')
title('3-ph input voltage');
axis( [0 360 -700  700])
legend('Vab','Vbc','Vca','Vba','Vcb','Vac')
xlabel('Angle (degrees)')
ylabel('Voltage (v)')

%for a=0
% load voltage
i=1;
for t=0:dt:T-dt
    if (mod(2*pi*f*t,2*pi) < pi/3) 
        V_L(i) = Vcb(i);
    elseif (mod(2*pi*f*t,2*pi) < 2*pi/3)
        V_L(i) = Vab(i);
    elseif (mod(2*pi*f*t,2*pi) < pi)
        V_L(i) = Vac(i);
    elseif (mod(2*pi*f*t,2*pi) < 4*pi/3)
        V_L(i) = Vbc(i);    
    elseif (mod(2*pi*f*t,2*pi) < 5*pi/3)
        V_L(i) = Vba(i);     
    else
        V_L(i) = Vca(i);      
    end
    i=i+1;
end
 
figure
hold on;
plot(w,V_L)
axis([0 360 -700  700])
title('3-ph Load Voltage a=0')
xlabel('Angle (degrees)')
ylabel('Voltage (v)')
 
%%%%% 3-phase Currents %%%%%
 
 
a = 0;    %gwnia enarxhs
R = 2.5   %Ohm
L1 = 0.04 %henry
L2 = 0.08 %henry

%matrixes for L=0.04 henry
A1 = -R/L1;
B1 = 1/L1;
C1 = 1;
D1 = 0;

%matrixes for L=0.08 henry
A2 = -R/L2;
B2 = 1/L2;
C2 = 1;
D2 = 0;

sys1=ss(A1,B1,C1,D1);

sys2 =ss(A2,B2,C2,D2);

sys1_d=c2d(sys1,dt);
sys2_d=c2d(sys2,dt);

%creating matrix
I_L1 = 0:1:500000-1;
I_L2 = 0:1:500000-1;

Ad1 = sys1_d.A;
Bd1 = sys1_d.B;
Cd1 = sys1_d.C;
Dd1 = sys1_d.D;

Ad2 = sys2_d.A;
Bd2 = sys2_d.B;
Cd2 = sys2_d.C;
Dd2 = sys2_d.D;


% load current
t=0:dt:T-dt;
for k = 1:1:500000-1
    I_L1(k+1) = Ad1*I_L1(k) + Bd1*V_L(k);   % current for L= 0.04 H
    I_L2(k+1) = Ad2*I_L2(k) + Bd2*V_L(k); % current for L= 0.08 H

end;

figure
plot(w,I_L1)
axis([16000 16360 0  250])
title('3-ph Load Current a=0 L=0.04')
xlabel('Angle (degrees)')
ylabel('Current (A)')

figure
plot(w,I_L2)
axis([16000 16360 0  250])
title('3-ph Load Current a=0 L=0.08')
xlabel('Angle (degrees)')
ylabel('Current (A)')



% thyristors currents for a=0
i=1;
for t=0:dt:T-dt
    
     % L =0.04 H
     Id1a(i) = 0;
     Id2a(i) = 0;
     Id3a(i) = 0;
     Id4a(i) = 0;
     Id5a(i) = 0;
     Id6a(i) = 0;
     
     % L =0.08 H
     Id1b(i) = 0;
     Id2b(i) = 0;
     Id3b(i) = 0;
     Id4b(i) = 0;
     Id5b(i) = 0;
     Id6b(i) = 0;
     
    if (mod(2*pi*f*t,pi) < pi/6) 
        Id1a(i) = I_L1(i);  % L =0.04 H
        Id6a(i) = I_L1(i);
        Id1b(i) = I_L2(i);  % L =0.08 H
        Id6b(i) = I_L2(i);   
    elseif (mod(2*pi*f*t,pi) < pi/3)    
        Id1a(i) = I_L1(i);
        Id2a(i) = I_L1(i);
        Id1b(i) = I_L2(i);
        Id2b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < pi/2)
        Id2a(i) = I_L1(i);
        Id3a(i) = I_L1(i);
        Id2b(i) = I_L2(i);
        Id3b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < 2*pi/3)
        Id3a(i) = I_L1(i);
        Id4a(i) = I_L1(i);
        Id3b(i) = I_L2(i);
        Id4b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < 5*pi/6)
        Id4a(i) = I_L1(i);
        Id5a(i) = I_L1(i);
        Id4b(i) = I_L2(i);
        Id5b(i) = I_L2(i);
   elseif (mod(2*pi*f*t,pi) < pi)
        Id5a(i) = I_L1(i);
        Id6a(i) = I_L1(i);
        Id5b(i) = I_L2(i);
        Id6b(i) = I_L2(i);
        
    end  

    i=i+1;
end
 
%ploting the currents of each thyristor for L = 0.04
figure

subplot(2,3,1);
plot(w,Id1a)
title('Current of thyristor T1 L = 0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,2);
plot(w,Id2a)
title('Current of thyristor T2 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,3);
plot(w,Id3a)
title('Current of thyristor T3 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,4);
plot(w,Id4a)
title('Current of thyristor T4 L=0.04');
axis([18000 18360 -5  250]);

subplot(2,3,5);
plot(w,Id5a)
title('Current of thyristor T5 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,6);
plot(w,Id6a)
title('Current of thyristor T6 L=0.04');
axis([18000 18360 -5  250]);

%ploting the currents of each thyristor for L = 0.08
figure

subplot(2,3,1);
plot(w,Id1b)
title('Current of thyristor T1 L=0.08');
axis([18000 18360 -5  250]);

subplot(2,3,2);
plot(w,Id2b)
title('Current of thyristor T2 L=0.08');
axis([18000 18360 -5  250]);

subplot(2,3,3);
plot(w,Id3b)
title('Current of thyristor T3 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,4);
plot(w,Id4b)
title('Current of thyristor T4 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,5);
plot(w,Id5b)
title('Current of thyristor T5 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,6);
plot(w,Id6b)
title('Current of thyristor T6 L=0.08');
axis([18000 18360 -5  250]);
 
%for a=67
 
%load voltage
i=1;
for t=0:dt:T-dt
    
   if (mod(2*pi*f*t,2*pi) < (7*(pi/180))) 
       V_L(i) = Vba(i);
   elseif (mod(2*pi*f*t,2*pi) < (67*(pi/180)))
       V_L(i) = Vca(i);
   elseif (mod(2*pi*f*t,2*pi) < (127*(pi/180)))
       V_L(i) = Vcb(i);
   elseif (mod(2*pi*f*t,2*pi) < (187*(pi/180)))
       V_L(i) = Vab(i);    
   elseif (mod(2*pi*f*t,2*pi) < (247*(pi/180)))
       V_L(i) = Vac(i);  
   elseif (mod(2*pi*f*t,2*pi) < (307*(pi/180)))
       V_L(i) = Vbc(i);
   else
       V_L(i) = Vba(i);      
   end
   i=i+1;
end
 
figure
plot(w,V_L)
axis([0 360 -100 700]);
title('3-ph Load Voltage a=67')
xlabel('Angle (degrees)')
ylabel('Voltage (v)')

% load current
t=0:dt:T-dt;
for k = 1:1:500000-1
    I_L1(k+1) = Ad1*I_L1(k) + Bd1*V_L(k);   % current for L= 0.04 H
    I_L2(k+1) = Ad2*I_L2(k) + Bd2*V_L(k); % current for L= 0.08 H

end;

figure
plot(w,I_L1)
axis([16000 16360 -100  250])
title('3-ph Load Current a=67 L=0.04')

figure
plot(w,I_L2)
axis([16000 16360 -100  250])
title('3-ph Load Current a=67 L=0.08')
xlabel('Angle (degrees)')
ylabel('Current (A)')
 
 %------------------------------------------------------------------------
 
% thyristors currents for a=67
i=1;
for t=0:dt:T-dt
     % L =0.04 H
     Id1a(i) = 0;
     Id2a(i) = 0;
     Id3a(i) = 0;
     Id4a(i) = 0;
     Id5a(i) = 0;
     Id6a(i) = 0;
     
     % L =0.08 H
     Id1b(i) = 0;
     Id2b(i) = 0;
     Id3b(i) = 0;
     Id4b(i) = 0;
     Id5b(i) = 0;
     Id6b(i) = 0;
     
    if (mod(2*pi*f*t,pi) < (7*(pi/180))) 
        Id5a(i) = I_L1(i); 
        Id6a(i) = I_L1(i);
        Id5b(i) = I_L2(i);  
        Id6b(i) = I_L2(i); 
        if(I_L1(i)<0)
            Id5a(i)=0;
            Id6a(i)=0;
        end;
        if(I_L2(i)<0)
            Id5b(i)=0;
            Id6b(i)=0;
        end;
    elseif (mod(2*pi*f*t,pi) < (37*(pi/180))) 
        Id1a(i) = I_L1(i);  
        Id6a(i) = I_L1(i);
        Id1b(i) = I_L2(i);  
        Id6b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < (67*(pi/180)))    
        Id1a(i) = I_L1(i);
        Id2a(i) = I_L1(i);
        Id1b(i) = I_L2(i);
        Id2b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < (97*(pi/180)))
        Id2a(i) = I_L1(i);
        Id3a(i) = I_L1(i);
        Id2b(i) = I_L2(i);
        Id3b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < (127*(pi/180)))
        Id3a(i) = I_L1(i);
        Id4a(i) = I_L1(i);
        Id3b(i) = I_L2(i);
        Id4b(i) = I_L2(i);
    elseif (mod(2*pi*f*t,pi) < (157*(pi/180)))
        Id4a(i) = I_L1(i);
        Id5a(i) = I_L1(i);
        Id4b(i) = I_L2(i);
        Id5b(i) = I_L2(i);
   elseif (mod(2*pi*f*t,pi) < (187*(pi/180)))
        Id5a(i) = I_L1(i);
        Id6a(i) = I_L1(i);
        Id5b(i) = I_L2(i);
        Id6b(i) = I_L2(i);
        
    end  

    i=i+1;
end
 
%ploting the currents of each thyristor for L = 0.04
figure
 
subplot(2,3,1);
plot(w,Id1a)
title('Current of thyristor T1 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,2);
plot(w,Id2a)
title('Current of thyristor T2 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,3);
plot(w,Id3a)
title('Current of thyristor T3 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,4);
plot(w,Id4a)
title('Current of thyristor T4 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,5);
plot(w,Id5a)
title('Current of thyristor T5 L=0.04');
axis([18000 18360 -5  250]);
 
subplot(2,3,6);
plot(w,Id6a)
title('Current of thyristor T6 L=0.04');
axis([18000 18360 -5  250]);

%ploting the currents of each thyristor for L = 0.08
figure
 
subplot(2,3,1);
plot(w,Id1b)
title('Current of thyristor T1 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,2);
plot(w,Id2b)
title('Current of thyristor T2 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,3);
plot(w,Id3b)
title('Current of thyristor T3 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,4);
plot(w,Id4b)
title('Current of thyristor T4 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,5);
plot(w,Id5b)
title('Current of thyristor T5 L=0.08');
axis([18000 18360 -5  250]);
 
subplot(2,3,6);
plot(w,Id6b)
title('Current of thyristor T6 L=0.08');
axis([18000 18360 -5  250]);
