function [dwmb]=densitywmb(k,r)

min=-0.5-10/(2*pi*k);
max=0.5+10/(2*pi*k);
mintvals=[floor(min):ceil(max)];
dwmb=0;
for h=1:length(mintvals)
    dwmb=dwmb+k^3/((2*pi)^0.5)*(2*pi*mintvals(1,h)-r)^2*exp(-k^2*(2*pi*mintvals(1,h)-r)^2/2);
end;
end
