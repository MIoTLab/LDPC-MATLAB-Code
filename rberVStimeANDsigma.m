clear all
clc
% According to the relationship between RBER and variance(sigma) and time
% respectively, determine variance(sigma) versus time.

% References:
%   [1] "Adaptive Read Thresholds for NANA Flash"
%   [2] "Optimizing NAND Flash-Based SSDs via Retention Relaxation"

mu=[-3 -1 1 3];
v=[-inf -2 0 2 inf];
syms t tm sigma

% [1], Page4, formula (2)
Q=@(X)1/sqrt(2*pi)*int(exp(-t^2/2),[X,inf]);

% f1(RBER) is computed by the mean(mu) and variance(sigma) of the Flash cells.
f1=3/2*(1-Q((mu(1)-v(2))/sigma)+Q((mu(2)-v(2))/sigma)+1-Q((mu(2)-v(3))/sigma)+...
    Q((mu(3)-v(3))/sigma)+1-Q((mu(3)-v(4))/sigma)+Q((mu(4)-v(4))/sigma));
% f1=3/2*int(Gauss(1,sigma,t),[-inf,0]); 

% [2], page4
rber_tmax=10^(-2); % 4.5*10^(-4);
tmax=1;  % 5
rbw=rber_tmax/300;
rbr=(rber_tmax-rbw)/tmax^1.25;   
f2=rbw+rbr*tm^1.25;

tm=0.2:0.1:1;   % 1:0.5:5;
figure
rber=subs(f2,tm);
plot(tm,rber,'*-')
xlabel('time/year')
ylabel('RBER')
grid on
saveas(gcf,'rberVStime.fig')

for i=1:length(rber)
    f(i)=f1-rber(i);
    sigma(i)=solve(f(i));
end

figure
plot(tm,sigma,'<-')
xlabel('time/year')
ylabel('sigma')
title0={'time/year';'RBER';'sigma'};
grid on
saveas(gcf,'sigmaVStime.fig')

data0=num2cell(double([tm;rber;sigma]));
data=[title0,data0];

% save sigma&rber_time2.mat data



