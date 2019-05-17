function [vlsb,vmsb] = decodeMS(rx,H,iteration,LLR)  % ,sigma ,PDindex
% disp('..........MS..........')

rx=real(rx);

[m,n]=size(H);


flag0=0; 
rx0=PAM4demod(rx);
vlsb=rx0(1:n);
vmsb=rx0((n+1):end);
check1=H*vlsb';



check2=H*vmsb';

for i=1:m
    if mod(check1(i,1),2)==1||mod(check2(i,1),2)==1
        flag0=1;
        break;
    end
end

if flag0==1
        
cw=sum(H);
rw=sum(H,2);
% Lp=zeros(1,n);
Lq=zeros(m,n);
Lr=zeros(m,n);
ncol=zeros(max(cw),n);
nrow=zeros(m,max(rw));

for i=1:n
    ncol(1:cw(i),i)=find(H(:,i));
end

for j=1:m
    nrow(j,1:rw(j))=find(H(j,:));
end

% [LLR,P]=compute_P_LLR(rx);  % _S1

for ii=1:2
    if ii==1
%         disp('   LSB ');
        Lp=LLR(1,:);
%         for a=1:n
%             Lp(1,a)=log((exp(-(rx(a)+1)^2/(2*sigma^2))+exp(-(rx(a)+3)^2/(2*sigma^2)))...
%                 /(exp(-(rx(a)-1)^2/(2*sigma^2))+exp(-(rx(a)-3)^2/(2*sigma^2))));
%         end
    else
%         disp('   MSB ');
        Lp=LLR(2,:);
%         for a=1:n
%             Lp(1,a)=log((exp(-(rx(a)-3)^2/(2*sigma^2))+exp(-(rx(a)+3)^2/(2*sigma^2)))...
%                 /(exp(-(rx(a)-1)^2/(2*sigma^2))+exp(-(rx(a)+1)^2/(2*sigma^2))));
%         end
    end
    
%     if PDindex==2  
%             Lp((m+fix((n-m)/2)+1):n)=-100;
%     end

    for i=1:n
        for jo=1:cw(i)
            Lq(ncol(jo,i),i)= Lp(1,i);
        end
    end
    

    for lp=1:iteration
%         fprintf('Iteration : %d\n', int8(lp));
        for i=1:m
            sign_product=1;
            for j0=1:rw(i)
                sign_product=sign_product*sign(Lq(i,nrow(i,j0))); % 消息节点对数概率信息符号累乘
            end
            for jo=1:rw(i)
                delta=1;
                for l=1:rw(i)
                    temp=Lq(i,nrow(i,l));
                    Lq(i,nrow(i,l))=1000;
                    fg=find(nrow(i,:)==0);
                    nrow(i,fg)=nrow(i,l);
                    delta=min(abs(Lq(i,nrow(i,:))));
                    nrow(i,fg)=0;
                    Lq(i,nrow(i,l))=temp;
                end
                Lr(i,nrow(i,jo))=sign_product.*sign(Lq(i,nrow(i,jo)))*delta;
            end
        end
        

        for i=1:n
            for jo=1:cw(i)
                alpha1=0;
                for l=1:cw(i)
                    if l~=jo
                        alpha1=alpha1+Lr(ncol(l,i),i);
                    end
                end
                Lq(ncol(jo,i),i) = Lp(1,i) + alpha1;
            end
        end

        
        LQ = zeros(1,n);
        for i=1:n
            delta2=0;
            for jo=1:cw(i)
                delta2=delta2+ Lr(ncol(jo,i),i);
            end
            LQ(1,i) = Lp(1,i) + delta2;
        end
        

        vHat=zeros(1,n);
        for i=1:n
            if  LQ(1,i)<0
                vHat(1,i)=1;
            end
        end
        
%         vlsb=vHat;
%         vmsb=vHat;
        check1=ones(m,1);
        check2=ones(m,1);
        
        if ii==1
            vlsb=vHat;
            check1=H*vlsb';
        elseif ii==2
            vmsb=vHat;
            check2=H*vmsb';
        end
        
        flag=0; 
        for i=1:m
            if mod(check1(i,1),2)==1||mod(check2(i,1),2)==1
                flag=1; 
                break;
            end
        end
        if flag==0
            break;
        end
        
    end
    
end

end
