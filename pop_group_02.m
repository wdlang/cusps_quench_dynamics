% population on the right-going group
clear all; close all; clc; myfont = 22;

L = 100;   N = 2*L+1;
ki = 50;
U = 12;
cutoff = 20; 
g = U/N;
Delta = 4*pi*sin(2*pi*ki/N)/N;
T = 2*pi/Delta;
theta = 2*atan(g*T);
steps = 1000;
dt = T/steps;
nloops = 20;
nsteps = nloops*steps;
plist = zeros(1, 1+nsteps);
plist_ana = zeros(1, 1+nsteps);

plist_kplus1 = zeros(1, 1+nsteps);
plist_kplus1_ana = zeros(1, 1+nsteps);

xlist = -L:L;
xlist = xlist';
psi0 = (1/sqrt(N))*exp(i*(2*pi*ki/N)*xlist);
blochgroup1 = zeros(N,2*cutoff+1);
for s1 = -cutoff : cutoff
    blochgroup1(:,s1 + cutoff + 1) = exp(i*2*pi*(ki+s1)/N*xlist)/sqrt(N);
end
blochgroup2 = blochgroup1';

H = zeros(N, N);
for s= 1:(N-1)
    H(s,s+1) = -1;     H(s+1,s) = -1;
end
H(1,N) = -1;  H(N,1) = -1;
H(L+1, L+1) = U;
[VV,DD] = eig(H);
dd = diag(DD);

plist(1) = 1;
plist_ana(1) = 1;
psi1 = VV'*psi0;
for s = 1:nsteps
    psi = VV*(exp(-i*dt*s*dd).*psi1);   
    proj = blochgroup2*psi;
    plist(s+1) = norm(proj).^2;
    plist_kplus1 ( s+1 ) = abs(proj(cutoff+2))^2;
   
    p = floor(s*dt/T);
    tdiff = s*dt - p*T;

    plist_ana(s+1) = cos(theta*p/2)^2 - g/sqrt(1+g^2*T^2)*sin(theta*(p+0.5))*tdiff;
     plist_kplus1_ana (s+1) = 4* g^2 * sin(Delta*tdiff/2)^2/ (1+ g^2 * T^2)/Delta^2;
end

h1 = figure;
ha = tight_subplot(2, 1, [0.02, 0.02], [0.32 0.01], [0.11 0.01]) ;
axes(ha(1));
plot(dt*(0:nsteps)/T, plist_kplus1, dt*(0:nsteps)/T, plist_kplus1_ana)
set(gca, 'fontsize', myfont)
set(ha(1),'XTickLabel','')
ylabel('$P $','fontsize',myfont,'Interpreter','latex');
axes(ha(2));
plot(dt*(0:nsteps)/T, plist, dt*(0:nsteps)/T, plist_ana)
hold on 
omega = theta/ T ; 
plot(dt*(0:nsteps)/T, 0.5*(1+ cos(omega*dt*(0:nsteps))), 'r:')
set(gca, 'fontsize', myfont)
xlim([0 nloops])
xlabel('$t/T$','fontsize',myfont,'Interpreter','latex');
ylabel('$P $','fontsize',myfont,'Interpreter','latex');

print(h1,'-depsc','kinks.eps')
