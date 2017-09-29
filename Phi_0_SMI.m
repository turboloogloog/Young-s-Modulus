clc;clear all;close all;
%% set parameter values for genarating y, that is g--------------------------------------------------
alfa = 3;
c = 2.3;
Lamda = 785e-9;
A0 = Lamda/2*2;            % Vibration amplitude (m)
ft = 40;                % Vibration frequency ( Hz)
fs = 20000;             % Sampling frequency (Hz)
n = fs/ft*8;            % Sampling length

%% -----------
t=(0:n-1)/fs;           % time index

phi_0=4*pi/Lamda*A0*sin(2*pi*ft*t);
% att=exp(-0.00006*(t1-0)); 
% att=exp(-60*t);
% phi_0=att.*phi_0;

%%
x=phi_0+atan(alfa);     % introduce x to simplify calculation
y=zeros(1,n);           % momery locating for solutions, where y = ph_f+atan(alfa)

% Now phase equation becomes y = x-c*sin(y)

precision = 10000;       % number of steps for searching within y_range
e = zeros(1,precision+1);
yy = e;
yyy = zeros(1,30);       % assume there are maximally 30 solutions for the phase equation

% this part determine the initial y(1) by solving the phase equation when x(1) is given
%
yy=x(1)-c+2*c*[1:precision]/precision;  % make the possible solutions within [x(1)-c, x(1)+c]
   e=yy-x(1)+c*sin(yy);  % yy=y_max*[1:precision]/precision; 
   e=sign(e);            %set p and N, determing the jump position  
   e=diff(e);  
   e=sign(e);  
   e=abs(e);             % e is either 1 or 0,jump part is 1, no jumpt part is zero.by 1 determine the jumping position
   iii=1;
   for ii=1:precision-1  % chech the values of e(y)
       if (1+c*cos(yy(ii)))>0;   % disgard the solutions with negative gradient
       if e(ii)==1; yyy(iii)=yy(ii); iii=iii+1; end   % solutions are these yy with zero e(ii),stable solutions r in yyy by increasing trend
       end
   end
   
   %   number_of_roots=iii-1;
   y(1)=yyy(1); % choose the smallest solution as the initial y(1). this will not affect the steady state behavior
 
% the following part for solving for y for a given x(i)
y_r=acos(-1/c);                         % the point where y jumps when x increases and reaches x_r
y_d=2*pi-acos(-1/c);                    % the point where x drops when x decreases and recahed x_d
x_r=y_r+sqrt(c*c-1);
x_d=y_d+sqrt(c*c-1);

for i=2:n      
   yy=x(i)-c+2*c*[1:precision]/precision;  % yy -- [x(i)-c, x(i)+c]
   e=yy-x(i)+c*sin(yy);  % yy - [-y_rang/2,y_rang/2] 
   e=sign(e);  
   e=diff(e);  
   e=sign(e);  
   e=abs(e);   % e is either 1 or 0
   iii=1;
   
   for ii=1:precision-1  % check the values of e(y)
       if (1+c*cos(yy(ii)))>0;   % disgard the solutions with negative gradient
       if e(ii)==1; yyy(iii)=yy(ii); iii=iii+1; end   % solutions are saved in yyy
       end
   end
   number_of_roots=iii-1;

   
   % the following part determine which solution is the true one based on linewidth mode selection principle
   % the mode is selected based on the following rules:
   %  1. y varies nomotically with x;
   %  2. y tries to keep continuty whenever possible;
   %  3. In case of jumpping and drooping, the mode with least linewidth will be selected

   ddx=x(i)-x(i-1);         
   ddy=abs(yyy(1)-y(i-1));  % distance of the smallest solution to previous y(i-1), yyy(1) is the smallest solution
   y0=yyy(1);   
   for iii=2:number_of_roots; 
       ddyy=abs(yyy(iii)-y(i-1));  % distance of yyy(i) to y(i)
       if ddyy<ddy; y0=yyy(iii); ddy=ddyy; end;  
   end 
   % Now y0 is the solution cloest to y(i-1) 

   % The following part determine the solution when jump or drop happens. 
   % The one with smallest distance to x will be chosen

   if abs(y0-y(i-1))>pi/4;   % condition for detecting the jumping points,jump threshold.
       dyx=abs(yyy(1)-x(i));
       for iii=1:number_of_roots;
           ddyx=abs(yyy(iii)-x(i));
           if ddyx<dyx;y0=yyy(iii);dyx=ddyx;end
       end
   end
   y(i)=y0;
end; 
y=y-atan(alfa);
g=cos(y);

%%
figure
NFFT = 2^nextpow2(n); % Next power of 2 from length of y
Y = fft(g,NFFT*20)/n;
f = fs/2*linspace(0,1,NFFT*20/2);

try_01 = ceil(NFFT/2*0.3); 
subplot(2,2,4)

x = f(1:try_01);y = 2*abs(Y(1:try_01));

plot(x,y)

title('Spectrum of g')
xlabel('Frequency (Hz)')
ylabel('Amplitude of Fourier transformation')


NFFTy = 2^nextpow2(n); % Next power of 2 from length of y
Yy = fft(phi_0,NFFT*20)/n;
fy = fs/2*linspace(0,1,NFFT*20/2);
subplot(2,2,2)

x=fy(1:try_01);y=2*abs(Yy(1:try_01));

plot(x,y)

title('Spectrum of phi_0')
xlabel('Frequency (Hz)')
ylabel('Amplitude of Fourier transformation')


subplot(2,2,1)
plot(t(1:n),phi_0);
title('Phase of laser without feedback')
xlabel('Times (s)')
ylabel('Amplitude of Phi_0 (rad)')
set(gca,'xlim',[0 max(t)]);
axis([0 max(t) -max(phi_0)*1.5 max(phi_0)*1.5]);
% grid on;


subplot(2,2,3)
% figure
plot(t,g);
title('Normalized self-mixing signal')
xlabel('Times (s)')
ylabel('Amplitude of G(t)')
set(gca,'xlim',[0 max(t)]);
axis([0 max(t) -1.5 1.5]);
% grid on;
