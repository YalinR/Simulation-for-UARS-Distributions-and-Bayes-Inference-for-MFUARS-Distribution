function [mm,Khatm,Shatm1,Shatm2,Shatm3,Shatm4,Shatm5,Shatm6,Shatm7,Shatm8,Shatm9]=getmatrix(k,a,b,g,delta,sigma,n,OS,m0,h)
%h is the # of rows of mm
%Shatm1 is for the (1,1) element in Shat
%Shatm2 is for the (1,2) element in Shat
%Shatm3 is for the (1,3) element in Shat
%so we have 10 plots for kappa and S estimates
mm=zeros(h,1);
Khatm=zeros(h,1);
Shatm1=zeros(h,1);
Shatm2=zeros(h,1);
Shatm3=zeros(h,1);
Shatm4=zeros(h,1);
Shatm5=zeros(h,1);
Shatm6=zeros(h,1);
Shatm7=zeros(h,1);
Shatm8=zeros(h,1);
Shatm9=zeros(h,1);
m=m0;
for c=1:h
    [BayesShat,BayesKhat]=mcmc(k,a,b,g,delta,sigma,n,OS,m);
    mm(c,1)=m;
    Khatm(c,1)=BayesKhat;
    Shatm1(c,1)=BayesShat(1,1);
    Shatm2(c,1)=BayesShat(1,2);
    Shatm3(c,1)=BayesShat(1,3);
    Shatm4(c,1)=BayesShat(2,1);
    Shatm5(c,1)=BayesShat(2,2);
    Shatm6(c,1)=BayesShat(2,3);
    Shatm7(c,1)=BayesShat(3,1);
    Shatm8(c,1)=BayesShat(3,2);
    Shatm9(c,1)=BayesShat(3,3);
    m=m+50;
end;
end

    
    
    
    
    
    
    
    


