function [h,kj,Shatm1,Shatm2,Shatm3,Shatm4,Shatm5,Shatm6,Shatm7,Shatm8,Shatm9]=getvector(k,a,b,g,delta,sigma,n,OS,m)
h=(1:m);
[~,~,Sj,kj]=mcmc(k,a,b,g,delta,sigma,n,OS,m);

Shatm1=zeros(m,1);
Shatm2=zeros(m,1);
Shatm3=zeros(m,1);
Shatm4=zeros(m,1);
Shatm5=zeros(m,1);
Shatm6=zeros(m,1);
Shatm7=zeros(m,1);
Shatm8=zeros(m,1);
Shatm9=zeros(m,1);
 
for c=1:m
Shatm1(c,1)=Sj(3*c-2,1);
Shatm2(c,1)=Sj(3*c-2,2);
Shatm3(c,1)=Sj(3*c-2,3);
Shatm4(c,1)=Sj(3*c-1,1);
Shatm5(c,1)=Sj(3*c-1,2);
Shatm6(c,1)=Sj(3*c-1,3);
Shatm7(c,1)=Sj(3*c,1);
Shatm8(c,1)=Sj(3*c,2);
Shatm9(c,1)=Sj(3*c,3);
end;
end
