% to plot the spectrum of the tight binding model 
clear all; close all; clc; 

x1 = (-1.4:0.01:1.4)*pi;
y1 = zeros(1, length(x1));

y2 = -2.5:0.01:2.5;
x2 = zeros(1, length(y2));

x3 = pi*(-1:0.01:1);
y3 = -2 * cos(x3);

x4 = pi*(-1:0.05:1);
y4 = -2*cos(x4);

y5 = 0:0.01:2;
x5 = -pi*ones(1, length(y5));

y6 = 0:0.01:2;
x6 = pi*ones(1, length(y5));

x7 = pi/2*(0.3:0.01:1.7);
y7 = 2*(x7-pi/2);

x8 = -pi/2*(0.3:0.01:1.7);
y8 = -2*(x8+pi/2);

h1 = figure;
plot (x1, y1, 'linewidth', 2)
hold on 
plot (x2, y2, 'linewidth', 2)
plot (x3, y3, 'linewidth', 2, 'color', 'r')
plot(x4, y4, 'o')

plot (x5, y5, ':', 'linewidth', 2, 'color', 'r')
plot (x6, y6,':', 'linewidth', 2, 'color', 'r')

plot (x7, y7,':', 'linewidth', 2, 'color', 'b')
plot (x8, y8,':', 'linewidth', 2, 'color', 'b')

r= 0.9;
plot(pi/2+r*cos(pi*(0:0.01:2)) , r*sin(pi*(0:0.01:2)), ':', 'linewidth', 2)
plot(-pi/2+r*cos(pi*(0:0.01:2)) , r*sin(pi*(0:0.01:2)), ':', 'linewidth', 2)

arrow([-1.4*pi 0],[pi*1.4 0],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')
arrow([0,-2.5],[ 0,2.5 ],'width',1.5,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')

text(-1.1*pi ,-0.35 ,'$- \pi $','fontsize',22,'Interpreter','latex')
text(0.95 *pi ,-0.35 ,'$  \pi $','fontsize',22,'Interpreter','latex')
text(1.3*pi ,-0.5 ,'$ q $','fontsize',22,'Interpreter','latex')
text(0.1*pi ,2.3, '$ \varepsilon $','fontsize',24,'Interpreter','latex')
text(0.5*pi ,-0.35 ,'$ q_i $','fontsize',22,'Interpreter','latex')
text(-0.70*pi ,-0.35 ,'$ -q_i $','fontsize',22,'Interpreter','latex')
% text(0.1*pi ,2 ,'$ 2 $','fontsize',22,'Interpreter','latex')

axis off 
% axis equal 

str = strcat('plot.eps');
print(h1,'-depsc',str)