function [sumo]=SUMO(n,OS)
%put all the observations in a big matrx which has dimension (3*n,3)
sumo=zeros(3);
sumo(1,:)=sum(OS(1:3:3*n-2,:));
sumo(2,:)=sum(OS(2:3:3*n-1,:));
sumo(3,:)=sum(OS(3:3:3*n,:));

end
    
