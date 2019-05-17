% function d=Entropy(sigma,theta)
clear all
clc
% Quantization design for the read voltage: Different reference voltage
% placements of different distributions under a certain theta parameter.

% Reference:
%   [1] "Read and Write Voltage Signal Optimization for Multi-Level-Cell (MLC) NAND Flash Memory"

syms X;
mu=[-3 -1 1 3];

sigma_cells = [0.3269 0.3228 0.3184 0.3135 0.3082 0.3023 0.2954 0.2873 0.2768];    % 'rberVStimeANDsigma.m'
% theta_set=[0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
theta=0.2;   % Reference [1], page 5.
x=-3:0.001:-2;

for ii=1:length(sigma_cells)
    sigma=sigma_cells(ii);
    
    f_s=0;
    for i=1:4
        f(i)=Gauss(mu(i),sigma,X);
        f_s=f_s+f(i);
    end
    
    h=0;
    for i=1:4
        h=h+f(i)/f_s*log2(f_s/f(i));
    end
    
    hh(ii,:)=subs(h,x);
    aa(ii,:)=double(hh(ii,:))-theta;
    ad=(find(abs(aa(ii,:))==min(abs(aa(ii,:)))));
    d(ii)=-2-x(ad);
    
end
