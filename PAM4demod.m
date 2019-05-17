function d = PAM4demod(r)
%4PAM demapping

d2=zeros(length(r),2);
for i = 1:length(r)
    % 00 01 11 10---- -3 -1 1 3
    d2(i,1) = r(1,i)>0;
    d2(i,2) = abs(r(1,i)) < 2;

end
d=[d2(:,1)',d2(:,2)'];

end