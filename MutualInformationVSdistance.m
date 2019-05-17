clc
clear all
% Draw the relationshipo between Mutual Information (MI) and the distance 
% from the intersection of two states to a certain point, and find the distance 
% value corresponding to the maximum mutual information (MMI).

% References:
%   [1] "Adaptive Read Thresholds for NANA Flash"

mu=[-3 -1 1 3];
sigma_set = [0.3269 0.3228 0.3184 0.3135 0.3082 0.3023 0.2954 0.2873 0.2768];   % 'rberVStimeANDsigma.m'

n=length(sigma_set);
newd=zeros(1,n);
b=zeros(n,4000);  
for ii=1:length(sigma_set) 
    disp(['......sigma',num2str(ii),'......'])
    sigma=sigma_set(ii);
    
    syms t d
    p=sym(zeros(7,4));   % p=sym(zeros(10,4));
    Q=@(X)1/sqrt(2*pi)*int(exp(-t^2/2),[X,inf]);    % Reference [1], page 4.
%     points=[-2-d -2 -2+d -d 0 d 2-d 2 2+d];     % 3 reads between 2 cells
    points=[-2-d -2+d -d d 2-d 2+d];     %  2 reads between 2 cells

    for i=1:2
        a_points=(mu(i)-points)/sigma;
        p(1,i)=Q(a_points(1));
        for j=2:6   %  j=2:9
            p(j,i)=Q(a_points(j))-Q(a_points(j-1));
        end
%         p(10,i)=1-Q(a_points(9));
        p(7,i)=1-Q(a_points(6));
    end
    
    p(:,3)=p(end:-1:1,2);
    p(:,4)=p(end:-1:1,1);
    
    % ——————mutual information——————
    hy_x=0;
    hy=0;
    
    for i=1:7  %% i=1:10 
        for j=1:4
            hy_x=hy_x+p(i,j)*log(p(i,j));
        end
        tp=sum(p(i,:))/4;
        hy=hy+tp*log(tp);
    end
    I=hy_x./4-hy;
    
    x=0.0005:0.0005:2;
    a=subs(I,x);
    b(ii,:)=real(double(a));
    newd(ii)=x(find(b(ii,:)==max(b(ii,:))));
    
end
figure
plot(x,b(1,:),x,b(2,:),x,b(3,:),x,b(4,:),x,b(5,:),x,b(6,:),x,b(7,:),x,b(8,:));  
xlabel('d')
ylabel('Mutual Information')
saveas(gcf,'MMIvsd.fig')




