% dynamics of the toy model
clear all; close all; clc; myfont = 22;

delta = 1;
g = 0.5;
T = 2*pi/delta;
dt = 0.001*T;
Tmax = 3*T;
tlist = 0:dt:Tmax;
Llist = [5, 10, 20]; 
plist = zeros(length(Llist), length(tlist));
plist2 = zeros(1, length(tlist));
plist3 = zeros(1, length(tlist));
plist4 = zeros(1, length(tlist));

for ss = 1: length(Llist)
    L = Llist(ss);
    N = 2*L +1;
    
    H = zeros(N, N);
    for s = 1: N
        H(s,s) = s*delta;
    end
    
    H = H + 2* g * ones(N, N);
    [VV,DD] = eig(H);
    dd = diag(DD);
    
    v0 = zeros(N, 1);
    v0 (L+1) = 1;
    
    v1 = VV'*v0;
    for s = 1: length(tlist)
        time = (s-1)*dt;
        v = VV* (exp(-i*time*dd).*v1);
        plist(ss,s) = abs(v(L+1))^2;
    end
end

for s = 1: length(tlist)
    time = (s-1)*dt;
    n1 = floor(time/T);
    time2 = time - n1*T;

    amp1 = 1/2;
    amp2 = 0.5*(((1-i*g*T)/(1+i*g*T))^n1)*(1-i*2*g*(time2- T/2))/(1+ i*g*T);
    plist2(s) = abs( 2* amp2 )^2;
    plist3(s) = abs(amp1 + amp2 )^2;
    plist4(s) = abs(-amp1 + amp2 )^2;
end

h1 = figure;
subplot('position', [0.15  0.74  0.65  0.25])
plot(tlist/T, plist(1,:),tlist/T, plist2,'--','linewidth',1.5)
ylim([0 1])
ylabel('$|\psi_0|^2 $','fontsize',myfont,'Interpreter','latex');
set(gca,'xticklabel',[])
XL=xlim; YL=ylim;
text(0.027*(XL(2)-XL(1))+XL(1),0.85*(YL(2)-YL(1))+YL(1),strcat('(a) $M$=',num2str(Llist(1))),'fontsize',22 , 'Interpreter','latex')
set(gca,'fontsize',22)
box on

subplot('position', [0.15  0.45  0.65  0.25])
plot(tlist/T, plist(2,:),tlist/T, plist2,'--','linewidth',1.5)
ylim([0 1])
ylabel('$|\psi_0|^2$','fontsize',myfont,'Interpreter','latex');
set(gca,'xticklabel',[])
XL=xlim; YL=ylim;
text(0.027*(XL(2)-XL(1))+XL(1),0.85*(YL(2)-YL(1))+YL(1),strcat('(b) $M$=',num2str(Llist(2))),'fontsize',22 , 'Interpreter','latex')
set(gca,'fontsize',22)
box on

subplot('position', [0.15  0.16  0.65  0.25])
plot(tlist/T, plist(3,:),tlist/T, plist2,'--','linewidth',1.5)
ylim([0 1])
xlabel('$t/T$','fontsize',myfont,'Interpreter','latex')
ylabel('$|\psi_0|^2 $','fontsize',myfont,'Interpreter','latex');
XL=xlim; YL=ylim;
text(0.027*(XL(2)-XL(1))+XL(1),0.85*(YL(2)-YL(1))+YL(1),strcat('(c) $M$=',num2str(Llist(3))),'fontsize',22 , 'Interpreter','latex')
set(gca,'fontsize',22)
box on
% 
% str = strcat('toy.eps');
% print(h1,'-depsc',str)