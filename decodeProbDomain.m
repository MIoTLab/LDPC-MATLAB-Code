function [vlsb,vmsb] = decodeProbDomain(rx, H,  iteration,Pr)  %  vHat ,k,ad ,PDindex,N0
% disp('..........BP..........')

rx=real(rx);

[m,n] = size(H);

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
 
% sigma=sqrt(N0/2);
vHat=zeros(1,n);

P=Pr;
for ii=1:2
    if ii==1
%         disp('   LSB ');
        P1=P(2,:);  % P1=Plsb1;
        P0=P(1,:);  % P0=Plsb0;
    else
%         disp('   MSB ');
        P1=P(4,:);  %  P1=Pmsb1;
        P0=P(3,:);  %  P0=Pmsb0;
    end
    
    % Initialization
    K0 = zeros(m, n);
    K1 = zeros(m, n);
    rji0 = zeros(m, n);
    rji1 = zeros(m, n);
    qij1 = H.*repmat(P1, m, 1);
    qij0 = H.*repmat(P0, m, 1);
    
    
    % Iteration
    for it = 1:iteration
%         fprintf('Iteration : %d\n', int8(it));
        % ----- Horizontal step -----
        for i = 1:m
            % Find non-zeros in the column
            c1 = find(H(i, :));
            for k = 1:length(c1)
                % Get column products of drji\c1(l)
                drji = 1;
                for l = 1:length(c1)
                    if l~= k
                        drji = drji*(qij0(i, c1(l)) - qij1(i, c1(l)));
                    end
                end % for l
                
                rji0(i, c1(k)) = (1 + drji)/2;
                rji1(i, c1(k)) = (1 - drji)/2;
            end % for k
        end % for i
        
        % ------ Vertical step ------
        for j = 1:n
            % Find non-zeros in the row
            r1 = find(H(:, j));
            for k = 1:length(r1)
                % Get row products of prodOfrij\ri(l)
                prodOfrij0 = 1;
                prodOfrij1 = 1;
                for l = 1:length(r1)
                    if l~= k
                        prodOfrij0 = prodOfrij0*rji0(r1(l), j);
                        prodOfrij1 = prodOfrij1*rji1(r1(l), j);
                    end
                end % for l
                
                % Update constants
                K0(r1(k), j) = P0(j)*prodOfrij0;
                K1(r1(k), j) = P1(j)*prodOfrij1;
                
                % Update qij0 and qij1
                qij0(r1(k), j) = K0(r1(k), j)./(K0(r1(k), j) + K1(r1(k), j));
                qij1(r1(k), j) = K1(r1(k), j)./(K0(r1(k), j) + K1(r1(k), j));
                
            end % for k
            
            % Update constants
            Ki0 = P0(j)*prod(rji0(r1, j));
            Ki1 = P1(j)*prod(rji1(r1, j));
            
            % Get Qj
            Qi0 = Ki0/(Ki0 + Ki1);
            Qi1 = Ki1/(Ki0 + Ki1);
            
            % Decode Qj
            if Qi1 > Qi0
                vHat(j) = 1;
            else
                vHat(j) = 0;
            end
            
        end % for j
        
        check1=ones(m,1);
        check2=ones(m,1);
        
        if ii==1
            vlsb=vHat;
            check1=H*vlsb';
        else
            vmsb=vHat;
            check2=H*vmsb';
        end
        
        % 译码尝试
%         check=H*vHat';
        flag=0;  %flag为0表示译码正确，为1表示译码错误
        for i=1:m
            if mod(check1(i,1),2)==1||mod(check2(i,1),2)==1
                flag=1;  % 若flag=1则继续迭代 lp+1
                break;
            end
        end
        if flag==0
%             disp(['BP迭代第',num2str(it),'次后译码正确']);
            break;
        end
        
    end % for n
end

end
% disp(['BP译码用时：',num2str(toc)]);
% disp('..........BP译码结束..........')

