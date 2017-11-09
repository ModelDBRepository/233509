
% Plot the membrane voltages of representative neurons in the thalamic
% network (generate subplots of Figure 2 of the Li et al. paper)
% The varialbe "FLAG_OSC" needs to be set to the corresponding simulated
% oscillation state so the figure is generated properly
% Written by Guoshi Li (guoshi_li@med.unc.edu)


clc;
clear all;
close all;

% Select which oscillation state to plot based on simulation
FLAG_OSC = 1; % 1: Delta; 2: Spindle; 3: Alpha: 4: Gamma

load tc1_all;
load tc2_all;
load in_all;
load re_all;


c1 = tc1_all;
c2 = tc2_all;
c3 = in_all;
c4 = re_all;

if (FLAG_OSC == 1)
  t1 = c1(:,1)-2000;
elseif (FLAG_OSC == 2)
  t1 = c1(:,1)-500;
elseif (FLAG_OSC == 3)
  t1 = c1(:,1)-950;  
else
  t1 = c1(:,1)-1000;   
end
  
v1 = c1(:,2);
v2 = c1(:,3);

v3 = c2(:,4);
v4 = c2(:,5);

if (FLAG_OSC == 4)
   v5 = c3(:,10);
   v6 = c3(:,11);
else
   v5 = c3(:,2);
   v6 = c3(:,3);    
end

v7 = c4(:,2);  
v8 = c4(:,3); 


if (FLAG_OSC == 1)
   xmin = 0;
   xmax = 1000;
   ymin = -100;
   ymax =  60;
elseif (FLAG_OSC == 2)
   xmin = 0;
   xmax = 3000;
   ymin = -100;
   ymax =  60;
elseif (FLAG_OSC == 3)
   xmin = 0;
   xmax = 1000;
   ymin = -80;
   ymax =  60;
else
   xmin = 0;
   xmax = 1000;
   ymin = -80;
   ymax = 60;
end


%==================================================================

figure;
subplot(4,1,1);
plot(t1,v1,'b',t1,v2,'r');
set(gca,'FontSize',12);
set(gca,'XTickLabel',[]);
axis([xmin, xmax, ymin, ymax]);
legend('HTC1', 'HTC2');
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
  

subplot(4,1,2);
plot(t1,v5,'b',t1,v6,'r');
set(gca,'FontSize',12);
set(gca,'XTickLabel',[]);
ylabel('mV', 'FontSize',14);
axis([xmin, xmax, ymin, ymax]);
legend('IN1', 'IN2');
box('off');

subplot(4,1,3);
plot(t1,v3,'b',t1,v4,'r');
set(gca,'FontSize',12);
set(gca,'XTickLabel',[]);
axis([xmin, xmax, ymin, ymax]);
legend('RTC1', 'RTC2');
box('off');

subplot(4,1,4);
plot(t1,v7,'b',t1,v8,'r');
set(gca,'FontSize',12);
xlabel('ms','FontSize',14);
axis([xmin, xmax, ymin, ymax]);
legend('RE1', 'RE2');
box('off');



