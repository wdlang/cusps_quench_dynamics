% 2016.jan.14 
clear all; close all; clc; myfont = 22;

r = 4;
s = 3.2;
theta = 1.1*2*pi/5;

x1 = -(r+1):0.01:(r+s);
y1 = zeros(1, length(x1));

x2 = y1;
y2 = x1;

list = 0:0.01:1; 

h1 = figure;

plot(x1, y1, 'linewidth', 2)
hold on 
plot(x2, y2, 'linewidth', 2)
plot(r*cos(2*pi*(0:0.01:1)), r*sin(2*pi*(0:0.01:1)), 'linewidth', 2 )

arrow(r*[1 0],r*[cos(-theta) -sin(theta)],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','r','EdgeColor','r')
arrow(r*[cos(theta) -sin(theta)],r*[cos(-2*theta) -sin(2*theta)],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','r','EdgeColor','r')
arrow(r*[cos(2*theta) -sin(2*theta)],r*[cos(-3*theta) -sin(3*theta)],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','r','EdgeColor','r')
arrow(r*[cos(3*theta) -sin(3*theta)],r*[cos(-4*theta) -sin(4*theta)],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','r','EdgeColor','r')
arrow(r*[cos(4*theta) -sin(4*theta)],r*[cos(-5*theta) -sin(5*theta)],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','r','EdgeColor','r')

plot( r*(1:-0.01:0), r*(0:-0.01:-1),'--','linewidth', 2,'color', 'g' )
plot( -r*(1:-0.01:0), -r*(0:-0.01:-1),'--','linewidth', 2,'color', 'g' )
plot( -r*(1:-0.01:0), r*(0:-0.01:-1),'--','linewidth', 2,'color', 'g' )
plot( r*(1:-0.01:0), -r*(0:-0.01:-1),'--','linewidth', 2,'color', 'g' )

arrow(r*[0.05 -0.95],r*[0 -1],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','g','EdgeColor','g')
arrow(r*[-0.95 -0.05],r*[-1 0],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','g','EdgeColor','g')
arrow(r*[-0.05 0.95],r*[0 1],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','g','EdgeColor','g')
arrow(r*[0.95 0.05],r*[1 0 ],'width',1,'TipAngle',15,'BaseAngle',30,'FaceColor','g','EdgeColor','g')

arrow([r 0],[r+s 0],'width',1,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')
arrow([0,r],[ 0,r+s ],'width',1,'TipAngle',24,'BaseAngle',30,'FaceColor','b','EdgeColor','b')

plot(r*cos(theta)*list, r*sin(-theta)*list,':' ,'linewidth', 2,'color', 'k')
r2 = 0.5;

plot(r2*cos(theta*list), r2*sin(-theta*list), 'linewidth', 2,'color', 'k')

text(0.5 ,-0.5 ,'$ \theta $','fontsize',22,'Interpreter','latex')

% text(r ,r ,'$ \psi_0 $','fontsize',22,'Interpreter','latex')
text(r+1 ,-1 ,'$ Re(\psi_0 )$','fontsize',22,'Interpreter','latex')
text(0.5 ,r+3 ,'$ Im(\psi_0 )$','fontsize',22,'Interpreter','latex')
text(r ,0.5 ,'$ (1,0)$','fontsize',22,'Interpreter','latex')
text(0 ,r+0.5 ,'$ (0,1)$','fontsize',22,'Interpreter','latex')

axis off 
axis equal 

str = strcat('trajectory.eps');
print(h1,'-depsc',str)