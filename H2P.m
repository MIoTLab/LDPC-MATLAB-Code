function [P,rearranged_cols]=H2P(H)

% To facilitate LDPC coding, transform the check matrix (H) into P form
% by gaussian transformayion.

dim=size(H);
rows=dim(1);
cols=dim(2);

rearranged_cols=zeros(1, rows);

for k=1:rows
    vec = [k:cols];

    x = k;
    while (x<=cols & H(k,x)==0)
        ind = find(H(k+1:rows, x) ~= 0);
        if ~isempty(ind)
            break                  
        end                        
        x = x + 1;
    end

    if x>cols    
        error('Invalid H matrix.');
    end

    if (x~=k)
        rearranged_cols(k)=x;
        temp=H(:,k);
        H(:,k)=H(:,x);
        H(:,x)=temp;
    end

    if (H(k,k)==0)
        ind = find(H(k+1:rows, k) ~= 0);
        ind_major = ind(1);
        x = k + ind_major;
        H(k, vec) = rem(H(x, vec) + H(k, vec), 2);
    end

    ind = find(H(:, k) ~= 0)';
    for i=1:length(ind)
        x = ind(i);
        if x ~= k
            H(x, vec) = rem(H(x, vec) + H(k, vec), 2);  
        end
    end
end

P=H(:,rows+1:cols);
