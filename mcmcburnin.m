function [Shat,S95,Khat,K025,K975]=mcmcburnin(k,a,b,g,delta,sigma,n,OS,m,burnin)
%  0<delta<pi
% -delta<r<delta
%  0=<a<=2*pi 0<=b<=pi  0<=g<=2*pi
%  m>burin
%  n is the number of 3 by 3 matrices in your data set
[sumo]=SUMO(n,OS);
Sj=zeros(3*m,3);
kj=zeros([m,1,1,1]);
%S is the starting value
S=[cos(a)*cos(g)-sin(a)*sin(g)*cos(b) sin(a)*cos(g)+cos(a)*sin(g)*cos(b) sin(g)*sin(b);-cos(a)*sin(g)-sin(a)*cos(g)*cos(b) -sin(a)*sin(g)+cos(a)*cos(g)*cos(b) cos(g)*sin(b);sin(a)*sin(b) -cos(a)*sin(b) cos(b)];
for c=1:m
%generate u1 u2 u3
z1=normrnd(0,1);
z2=normrnd(0,1);
z3=normrnd(0,1);
zd=sqrt(z1^2 +z2^2+z3^2);

u1=z1/zd;
u2=z2/zd;
u3=z3/zd;
%we can obtain the cdf, then by inverse the cdf we can get r 
%-delta<r<delta
x=rand;
if 2*x-1<0
    r= -delta*(1-2*x)^(1/3);
else r= delta*(2*x-1)^(1/3);
end;

%Sstar is generated value by the proposal
M=[u1^2+cos(r)-u1^2*cos(r) u1*u2-u1*u2*cos(r)-u3*sin(r) u1*u3-u1*u3*cos(r)+u2*sin(r);u1*u2-u1*u2*cos(r)+u3*sin(r) u2^2+cos(r)-u2^2*cos(r) u2*u3-u2*u3*cos(r)-u1*sin(r);u1*u3-u1*u3*cos(r)-u2*sin(r) u2*u3-u2*u3*cos(r)+u1*sin(r) u3^2+cos(r)-u3^2*cos(r)];
Sstar=S*M;

%obtain r1
rj1=exp(k*(trace(Sstar.'*sumo)-trace(S.'*sumo)));
%generate Bernoulli w
wj1 = rand(1,1) <= min(1,rj1);
Sj(3*c-2:3*c,:)= wj1*Sstar + (1-wj1)*S;

%generate kappa
%sigma is the variance parameter for log normal
S=Sj(3*c-2:3*c,:);
kstar=exp(normrnd(log(k),sigma));
[gkstar]=posteriordensityMFUARS(kstar,S,n,OS);
[gk]=posteriordensityMFUARS(k,S,n,OS);
rj2=(gkstar*kstar)/(gk*k);
wj2= rand(1,1) <= min(1,rj2);
kj(c,1)= wj2*kstar+(1-wj2)*k;
k=kj(c,1);
end;
%get the Bayes estiamte for S
[Sjsum]=SUMO(m-burnin,Sj(burnin*3+1:3*m,:));
Sbar=Sjsum/(m-burnin);
[Q,D]=eig(Sbar*Sbar.');
Shat=(Q*D^(0.5)*Q^(-1))^(-1)*Sbar;

%credible set for S
d=zeros([m-burnin,1,1,1]);
SS=Sj(burnin*3+1:3*m,:);
for c=1:(m-burnin)
d(c,1)=max(acos(diag(Shat.'*SS(3*c-2:3*c,:))));
end;
S95=quantile(d,0.95);

%obtain the Bayes estimate for kappa
Khat=sum(kj(burnin+1:m,:))/(m-burnin);

%credible set for kappa
K025=quantile(kj(burnin+1:m,:),0.025);
K975=quantile(kj(burnin+1:m,:),0.975);

end



