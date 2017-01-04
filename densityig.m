function [dIG,sum,c]=densityig(k,r)
c=round(-(k^2)*log(0.6)-1)+100;
while (2*c+1)*exp((-c^2-c)/(2*k^2))/(1-5/3*exp((c+1)/(-k^2)))>0.000000000001
c=c+1;
end;
sum=0;

for h=0:c
sum=sum+(2*h+1)*exp(-h*(h+1)/(2*k^2))*sin((h+0.5)*r);
end;   
dIG=(sin(r/2)/pi)*sum;

end