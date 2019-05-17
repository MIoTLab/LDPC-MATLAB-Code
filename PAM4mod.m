function r = PAM4mod (s)
% 4-PAM modulation

r=zeros(1,length(s));
for i = 1:length(s)
    % 00 01 11 10---- -3 -1 1 3
    r(1,i) = (2*s(i,1) -1) * (2*(s(i,2)==0)+1); %4PAM mapping

end

end
