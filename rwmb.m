function [rval]=rwmb(k,n)
rval=zeros([n,1,1,1]);
min=-0.5-10/(2*pi*k);
max=0.5+10/(2*pi*k);
mintvals=floor(min):ceil(max);
for c=1:n
    %the range for m
    %cc=(ncx2inv(0.9999999999999999,3,0))^0.5=8.7975
    %use the number 10 which is bigger than 8.7975

    %generate a random variable from chisq df=3
    rchi=chi2rnd(3);
    if rand>0.5
     rchir=(rchi)^0.5;
    else
     rchir=-(rchi)^0.5;
    end;
    %transform the rchi to a r so that -pi<r<=pi
    rcandidate=zeros([length(mintvals),1,1,1]);
    for h=1:length(mintvals)
        rcandidate(h,1)=mintvals(1,h)*(2*pi)+ rchir/k;
        if -pi<rcandidate(h,1) && rcandidate(h,1)<=pi
        rval(c,1)=rcandidate(h,1);
        end;
    end;
end;
end
    
