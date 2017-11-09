
% Plot the simulated LFP with frequency power spectrum of the thalamic network
% The m-file generates subplots of Figure 2 of the Li et al. paper
% The varialbe "FLAG_OSC" needs to be set to the corresponding simulated
% oscillation state so the figure is generated properly
% Written by Guoshi Li (guoshi_li@med.unc.edu)

clc;
clear all;
close all;

% Select which oscillation state to plot based on simulation
FLAG_OSC = 1; % 1: Delta; 2: Spindle; 3: Alpha: 4: Gamma

load tc1_all;
load tc2_all
load in_all
load re_all;
 
if (FLAG_OSC == 1)
   T0 = 2000;    
   T1 = T0+1000;
elseif (FLAG_OSC == 2)
   T0 = 1100;    
   T1 = T0+1000;
elseif (FLAG_OSC == 3)
   T0 = 1000;    
   T1 = T0+1000;
else
   T0 = 1000;    
   T1 = T0+1000;
end


C1 = tc1_all(:, 2:end);
C2 = tc2_all(:, 2:end);

C = [C1 C2];

FILORDER = 1000;

[row, col]=size(C);

TC = C;
lfp = sum(TC,2)/(col);  % simulated local field potential


DT = 0.2;               % sampling time: ms
Fs = 1/DT*1000;         % sampling frequency: Hz

Fmax = 50;              % maximal frequency to plot
Fc   = [0.5 80];        % Cut-off frequency
Wc   = Fc/(Fs/2);   

n1 = T0/DT+1;
n2 = T1/DT+1;

t = tc1_all(:,1);

if (FLAG_OSC == 1)
   t = t(n1:n2)-T0;
elseif (FLAG_OSC == 2)
   t = t(n1:n2)-500;
elseif (FLAG_OSC == 3)
   t = t(n1:n2)-T0;
else
   t = t(n1:n2)-T0;
end


y = lfp(n1:n2);
y = y-mean(y);

L = length(y);
NFFT = 2^nextpow2(L);     % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2));

f = Fs/2*linspace(0,1,NFFT/2);

m = Fmax/(0.5*Fs)*(0.5*NFFT);
m = ceil(m)+1;

h = fir1(FILORDER, Wc);
x = filtfilt(h,1, y);

X = fft(x,NFFT)/L;
XX = 2*abs(X(1:NFFT/2));


[Peak, I] = max(XX);
fo=f(I);
disp('The oscillation frequency is:');
fo
disp('The oscillation power is:');
Peak


if (FLAG_OSC == 1)
   xmin  = 0;
   xmax  = 1000;
   ymin  = -20;
   ymax  = 50;
   ypmax = 11;
elseif (FLAG_OSC == 2)
   xmin = 600;
   xmax = 1600;
   ymin = -20;
   ymax = 50;
   ypmax = 20;
elseif (FLAG_OSC == 3)
   xmin = 0;
   xmax = 1000;
   ymin = -25;
   ymax =  25;
   ypmax = 11;
else
   xmin = 0;
   xmax = 1000;
   ymin = -25;
   ymax =  25;  
   ypmax = 11;
end


figure;
subplot(2,1,1);
plot(t, x, 'k-', 'LineWidth',1);
xlabel('ms', 'FontSize',14);
ylabel('sLFP (mV)', 'FontSize',14);
set(gca, 'FontSize',12);
set(gca, 'YTick', [-20:20:40]);
axis([xmin, xmax, ymin, ymax]);
box('off');

if (FLAG_OSC == 1)
  title('Delta OSC', 'FontSize',16);
elseif (FLAG_OSC == 2)
  title('Spindle OSC', 'FontSize',16);
elseif (FLAG_OSC == 3)
  title('Alpha OSC', 'FontSize',16); 
else
  title('Gamma OSC', 'FontSize',16);  
end


% Plot single-sided amplitude spectrum.
subplot(2,1,2);
plot(f(1:m),XX(1:m),'k-');
xlabel('Frequency (Hz)', 'FontSize',14)
ylabel('Power', 'FontSize',14)
set(gca, 'FontSize',12);
axis([0, Fmax, 0, ypmax]);
box('off');







