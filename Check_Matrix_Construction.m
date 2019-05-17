% This program constructs girth 8 quasi-cyclic LDPC codes given
% row weight(k),column-weight (j) and size of each sub-matrix (m). Could 
% also be used to construct girth 6 codes.
% The program calculates the shift value of each check or varibale group.
% Check and variable node connections can also be varied to obtain regular,
% irregular and codes with zero-submatrices. The constructed code (parity-check)
%is stored in H.
% Author : Dr. Gabofetswe Alafang Malema,University of Botswana  
% Department of Computer science. e-mail: malemag@mopipi.ub.bw

%NOTE: this is a searh algorithm it may take many tries to obtain a code.
%it may also fail to find a code given some parameters even when
%theoretical it is possible to get a code. 
clear all
clc

j=4;k=36; % j is column-weight, k is row-weight (they are both variable)
m=256; % size of check and variable groups(or submatrix). Also variable.  % m=1024
num_gps = j+k; % this is the default number of groups (check and varible
                 % node groups) for regular codes.The number of check and
                 % variable groups could be changed.
M=j*m;           % M is the total number of check nodes. the formula 
                 % depends on the number check node groups used.
gps = j*k;     % this is the number of check-variable node group 
                 % connections. it depends on the number of groups and
                 % connections you want to make.
shift=zeros(1,gps); % holds shift values for each check-variable node connection.
rows = struct('connect',0,'counter',0,'con',0); % data structure used for
emptyset = [];                                  % holding the connections.
group = struct('mem',0,'grp',0,'con',0);

for y = 1:m*num_gps
    rows(y).counter = 0;
    rows(y).con = intersect(emptyset,rows(y).con);
    rows(y).connect = intersect(emptyset,rows(y).connect);
    
end

for i = 1:num_gps                   % populate the contents of each group.
    group(i).mem = [m*(i-1)+1:i*m];
end

% The following for loops connect check nodes and variable nodes.
% Regular and non-zero matrices are assumed. For irregular and zero
% matrices a custom variable to check node connection must be specified. The
% algorithm still works in the same way.
counter =0;
for i=j+1:k+j
    for y=1:j
        counter=counter+1;
        group(counter).grp=[y i];
    end
end

 % regular general construction. Example of how 3 check nodes groups
 % are connected to the 6 variable node groups (4 to 9). The default 
 % connection configuration above is as below when specified manually. The
 % configuration can be changed to suit your needs.
%  group(1).grp = [1 4]; group(10).grp=[1 7];
%  group(2).grp = [2 4]; group(11).grp=[2 7];
%  group(3).grp = [3 4]; group(12).grp=[3 7];
%  group(4).grp = [1 5]; group(13).grp=[1 8];
%  group(5).grp = [2 5]; group(14).grp=[2 8];
%  group(6).grp = [3 5]; group(15).grp=[3 8];
%  group(7).grp = [1 6]; group(16).grp=[1 9];
%  group(8).grp = [2 6]; group(17).grp=[2 9];
%  group(9).grp = [3 6]; group(18).grp=[3 9];


  g= 2; % determines the girth of the code. For girth 8 use g=4 and 
        % g=2 for girth 6 codes.
s_counter = 1;done=1;
for groups = 1:gps   % central groups.
    
    current_gp = group(groups).grp(1,1);
    r=1;
        i = group(current_gp).mem(r);
        grp2 = group(groups).grp(1,2);
                                         
            mem=[]; mem1 = rows(i).con;
            mem2 = [];
            for y = 1:g
                
                for x = 1:length(mem1)
                    x1 = mem1(x);
                    mem2 = union(mem2,rows(x1).con);
                end
                mem = union(mem1,mem2);
                mem1=mem;
            end
                   
                A = intersect(mem,group(grp2).mem);
                A = setxor(A,group(grp2).mem);
                

     if (isempty(A)~=1)
         r1 = ceil(rand*length(A));
         x1 = A(r1); 
         rows(i).counter = rows(i).counter + 1;
         rows(x1).counter = rows(x1).counter + 1;
                    
        a = rows(i).counter;
        c = rows(x1).counter;
                    
        rows(i).connect(a,1)=x1;
                    
        rows(i).con = union(rows(i).con,x1);
                    
        rows(x1).connect(c,1)=i;
        rows(x1).con = union(rows(x1).con,i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        x = find(group(grp2).mem == x1);
        shift(s_counter)=x;
        s_counter = s_counter + 1;
        
        p1 =circshift(group(grp2).mem,[1 -x+1]);
        index = rows(i).counter;
        for t1 = 2:m
            t = group(current_gp).mem(t1);
            rows(t).connect(index,1)=p1(t1);
            
            c1 =rows(p1(t1)).counter + 1;
                       
            rows(p1(t1)).connect(c1,1) = t;
            
            rows(p1(t1)).con = union(rows(p1(t1)).con,t);
            rows(t).con = union(rows(t).con,p1(t1));
            rows(t).counter = rows(t).counter+1;
            
            rows(p1(t1)).counter = c1;
            
        end
     else
         disp('Column not found. Code not obtained. Try again or increase m.');
         done=0;
    end   
    
end % groups for loop.

if done ==1
    H=zeros(M,m);
    for i=1:M            % constructs the code matrix based on check and
        for x =1:length(rows(i).connect) % variable node connections above.
            r=rows(i).connect(x)-M; 
            H(i,r)=1;
        end
    end
end
    

