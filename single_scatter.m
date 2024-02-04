% 2d tight binding model; scattering of a plane wave against a single
% scatter
clear all; close all; clc; tic;

N = 60;
dim = N*N;
W = 4;
kx = 2;
ky = 0;
dt = 50;
tlist = 0: dt: 20000;

H = zeros(dim, dim);
for s2= 0:N-1
    for s1 = 0:N-1
        s3 = s1; s4 = mod(s2-1,N);
        H(s1+1+N*s2, s3+1+N*s4)=-1;
        s3 = s1; s4 = mod(s2+1,N);
        H(s1+1+N*s2, s3+1+N*s4)=-1;
        s3 = mod(s1-1,N); s4 = s2;
        H(s1+1+N*s2, s3+1+N*s4)=-1;
        s3 = mod(s1+1,N); s4 = s2;
        H(s1+1+N*s2, s3+1+N*s4)=-1;
    end
end
H(1, 1)= W;

[V,D]= eig(H);
dd = diag(D);
clear D

ini = zeros(N, N);
for s2= 0:N-1
    for s1 = 0:N-1
        ini(s1+1, s2+1)=exp(i*(kx*s1+ ky*s2));
    end
end

inif = reshape(ini, N*N, 1);
inif2 = V'*inif;
vec = inif;
wave = reshape(vec, N, N);
waveK = fftshift(fft2(wave));
denK = abs(waveK).^2;

h=figure;
% set(gcf,'doublebuffer','off');

% clims1=[0 max(max(denK))];
% himage1=imagesc([-N/2 N/2],[-N/2 N/2],denK,clims1);
surf(0:N-1, 0:N-1, denK)
% set(gca,'YDir','normal')
xlabel('$k_1$','fontsize',20,'Interpreter','latex')
ylabel('$k_2$','fontsize',20,'Interpreter','latex')
% axis image
set(gca,'fontsize',16)
str=strcat('W=',num2str(W));
title(str,'fontsize',16)

pause(0.01)

for s=1:length(tlist)
    s
    vec = V* (exp(-i*tlist(s)*dd).*inif2);
    wave = reshape(vec, N, N);
    waveK = fftshift(fft2(wave));
    denK = abs(waveK).^2;
    
%     refreshdata(himage1)
%     set(himage1,'CDATA',denK)
%     
%     pause(0.01)

% clims1=[0 max(max(denK))];
% himage1=imagesc([-N/2 N/2],[-N/2 N/2],denK,clims1);
% set(gca,'YDir','normal')
surf(0:N-1, 0:N-1, denK)
xlabel('$k_1$','fontsize',20,'Interpreter','latex')
ylabel('$k_2$','fontsize',20,'Interpreter','latex')
% axis image
set(gca,'fontsize',16)
str=strcat('W=',num2str(W));
title(str,'fontsize',16)

pause(0.02)
end