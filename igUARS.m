function [O,S]=igUARS(k,a,b,g,n)

%decide the number 'aa' first for rejection sampling
rvec=(-pi:0.01:pi);
result=zeros([length(rvec),1,1,1]);

for c=1:length(rvec)
qr=densityig(k,rvec(c));
gr=densitywmb(k,rvec(c));
result(c,1)=qr/gr;
end;
aa=max(result);

O=zeros(3*n,3);
S=[cos(a)*cos(g)-sin(a)*sin(g)*cos(b) sin(a)*cos(g)+cos(a)*sin(g)*cos(b) sin(g)*sin(b);-cos(a)*sin(g)-sin(a)*cos(g)*cos(b) -sin(a)*sin(g)+cos(a)*cos(g)*cos(b) cos(g)*sin(b);sin(a)*sin(b) -cos(a)*sin(b) cos(b)];
for i=1:n
z1=normrnd(0,1);
z2=normrnd(0,1);
z3=normrnd(0,1);
zd=sqrt(z1^2 +z2^2+z3^2);
u1=z1/zd;
u2=z2/zd;
u3=z3/zd;
%by rejection sampling to generate r from IG, the proposal is wmb
rn=rwmb(k,1);
x=rand;
while x>=densityig(k,rn)/(aa*densitywmb(k,rn))
    x=rand;
    rn=rwmb(k,1);
end;
r=rn;
M=[u1^2+cos(r)-u1^2*cos(r) u1*u2-u1*u2*cos(r)-u3*sin(r) u1*u3-u1*u3*cos(r)+u2*sin(r);u1*u2-u1*u2*cos(r)+u3*sin(r) u2^2+cos(r)-u2^2*cos(r) u2*u3-u2*u3*cos(r)-u1*sin(r);u1*u3-u1*u3*cos(r)-u2*sin(r) u2*u3-u2*u3*cos(r)+u1*sin(r) u3^2+cos(r)-u3^2*cos(r)];
O(3*i-2:3*i,:)=S*M;
end
