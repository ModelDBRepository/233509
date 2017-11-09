
% Plot the spiking activiites of the whole thalamic network
% The m-file generates subplots of Figure 2 of the Li et al. paper
% The varialbe "FLAG_OSC" needs to be set to the corresponding simulated
% oscillation state so the figure is generated properly
% Written by Guoshi Li (guoshi_li@med.unc.edu)


clc;
clear all;
close all;

% Select which oscillation state to plot based on simulation
FLAG_OSC = 1; % 1: Delta; 2: Spindle; 3: Alpha: 4: Gamma


min_T = 0;

if (FLAG_OSC == 1)
  T0 = 1000;
  T1 = T0+1000;
  max_T = 1000;
elseif (FLAG_OSC == 2)
  T0 = 500;
  T1 = T0+3000;
  max_T = 3000;
elseif (FLAG_OSC == 3)
  T0 = 950;
  T1 = T0+1000; 
  max_T = 1000;
else
  T0 = 1000;
  T1 = T0+1000;  
  max_T = 1000;
end


Ntc1x = 7;
Ntc1y = 7;

Ntc2x = 12;
Ntc2y = 12;

Nin1 = 8;
Nin2 = 8;

Nre1 = 10;
Nre2 = 10;

Ntc1  = Ntc1x*Ntc1y;     
Ntc2  = Ntc2x*Ntc2y;    
Nin   = Nin1*Nin2;
Nre   = Nre1*Nre2;   

dh = 0.3;
di = 0.3;
dr = 0.3;
de = 0.3;


T_Start = T0;            % Start time of the rasterplot     
T_End   = T1;            % End time of the rasterplot

T  = T_End - T_Start;    % Total duration in ms


TO1 = T_Start;    
TO2 = T_End;       

TO  = (TO2-TO1)/1000;


%============================================
%        Generate raster plot 
%============================================
% For TC1 cells
figure;

subplot(4,1,1);

for i = 0:1:(Ntc1x-1)
   for j = 0:1:(Ntc1y-1) 
       
   n = i*Ntc1y+j+1;
   
   s = ['load TC1' '_' int2str(i) '_'  int2str(j) ';'];    
   eval(s);
    
   ss = ['SpkT = TC1' '_'  int2str(i) '_'  int2str(j) ';'];    
   eval(ss);  
   
   % Firing rate 
    A = find (SpkT>=TO1 & SpkT<=TO2); 
    SHTC(n,1) = length(A)/TO;
   
   L = length(SpkT);
   if (L~=0)  
    for k = 1:L
     if (SpkT(k)>=T_Start & SpkT(k)<=T_End)
      
      x = [SpkT(k)-T0  SpkT(k)-T0];
      y = [n-dh        n+dh ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end
    end
   end
   
 end
end

% xlabel('ms', 'FontSize',14);
ylabel('HTC #', 'FontSize',14);
set(gca, 'FontSize',12);
set(gca,'XTickLabel',[]);
set(gca, 'YTick', [0:25:50]);
axis([min_T,max_T,0,Ntc1+1]);
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



% Calculate Average firing rate
FHTC = SHTC/T*1000;
fHTC = mean(FHTC);
disp('The average HTC firing rate is:');
fHTC

clear SpkT;


% For Interneurons
subplot(4,1,2);

for i = 0:1:(Nin1-1)
   for j = 0:1:(Nin2-1) 
       
    n = i*Nin2+j+1;
    
    s = ['load IN' '_' int2str(i) '_'  int2str(j) ';'];    
    eval(s);
   
    ss = ['SpkT = IN' '_'  int2str(i) '_'  int2str(j) ';'];    
    eval(ss);  
   
   % Odor rate 
    A = find (SpkT>=TO1 & SpkT<TO2); 
    SIN(n,1) = length(A)/TO;
   
   L = length(SpkT);
   if (L~=0)  
    for k = 1:L
     if (SpkT(k)>=T_Start & SpkT(k)<=T_End)
      
      x = [SpkT(k)-T0  SpkT(k)-T0];
      y = [n-di        n+di ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end
    end
   end
   
 end
end

% xlabel('ms', 'FontSize',14);
ylabel('IN #', 'FontSize',14);
set(gca, 'FontSize',12);
axis([min_T,max_T,0,Nin+1]);
set(gca,'XTickLabel',[]);
set(gca, 'YTick', [0:30:60]);
box('off');


% Calculate Average firing rate
FIN = SIN/T*1000;
fIN = mean(FIN);
disp('The average IN firing rate is:');
fIN

clear SpkT;


% For TC2 cells
subplot(4,1,3);

for i = 0:1:(Ntc2x-1)
   for j = 0:1:(Ntc2y-1) 
       
    n = i*Ntc2y+j+1;
    
    s = ['load TC2' '_' int2str(i) '_'  int2str(j) ';'];    
    eval(s);
    
    ss = ['SpkT = TC2' '_'  int2str(i) '_'  int2str(j) ';'];    
    eval(ss);  
     
   % Odor rate 
    A = find (SpkT>=TO1 & SpkT<TO2); 
    SRTC(n,1) = length(A)/TO;
   
   L = length(SpkT);
   if (L~=0)  
    for k = 1:L
     if (SpkT(k)>=T_Start & SpkT(k)<=T_End)
      
      x = [SpkT(k)-T0 SpkT(k)-T0];
      y = [n-dr       n+dr ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end
    end
   end
   
 end
end

% xlabel('ms', 'FontSize',14);
ylabel('RTC #', 'FontSize',14);
set(gca, 'FontSize',12);
axis([min_T,max_T,0,150]);
set(gca,'XTickLabel',[]);
box('off');

% Calculate Average firing rate
FRTC = SRTC/T*1000;
fRTC = mean(FRTC);
disp('The average RTC firing rate is:');
fRTC


clear SpkT;


%=========================================================

% For RE cells
subplot(4,1,4);

for i = 0:1:(Nre1-1)
   for j = 0:1:(Nre2-1) 
       
    n = i*Nre2+j+1;
    
    s = ['load RE' '_' int2str(i) '_'  int2str(j) ';'];    
    eval(s);
    
    ss = ['SpkT = RE' '_'  int2str(i) '_'  int2str(j) ';'];    
    eval(ss);  
   
   % Odor rate 
    A = find (SpkT>=TO1 & SpkT<TO2); 
    SRE(n,1) = length(A)/TO;
   
   L = length(SpkT);
   if (L~=0)  
    for k = 1:L
     if (SpkT(k)>=T_Start & SpkT(k)<=T_End)
      
      x = [SpkT(k)-T0  SpkT(k)-T0];
      y = [n-de        n+de ];
   
       plot(x,y,'k','LineWidth',1);
       hold on;
     end
    end
   end
   
 end
end

xlabel('ms', 'FontSize',14);
ylabel('RE #', 'FontSize',14);
set(gca, 'FontSize',12);
set(gca, 'YTick', [0:50:100]);
axis([min_T,max_T,0,Nre+1]);
box('off');

% Calculate average RE firing rate
FRE = SRE/T*1000;
fRE = mean(FRE);
disp('The average RTC firing rate is:');
fRE




