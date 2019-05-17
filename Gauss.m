function F=Gauss(mu,sigma,x)
% Gauss function
F=1/(sqrt(2*pi)*sigma).*exp(-(x-mu).^2/(2*sigma^2));


