function [O,S]=MFUARSfinal(k,a,b,g,n)
%if k is less than or equal to 1, use the wtn as the proposal
%if k is greater than 1, use the vonmise as the proposal
if (k<=1)
    [O,S]=MFUARS(k,a,b,g,n);
elseif (k>1)
    [O,S]=MFUARSv(k,a,b,g,n);
end;
end