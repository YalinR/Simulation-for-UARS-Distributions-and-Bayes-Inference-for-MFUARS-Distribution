function [dcayley]=densitycayley(k,r)

dcayley=(1-cos(r))/(2*pi)*(pi^0.5)*gamma(2*k^2+2)*(1+cos(r))^(2*k^2)/(2^(2*k^2)*gamma(2*k^2+0.5));

end