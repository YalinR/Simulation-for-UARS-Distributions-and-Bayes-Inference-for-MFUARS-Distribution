function [dlor]=densitylorentz(k,r)

la=k/2-0.5+2/((k+2)^2);
dlor = (1-cos(r))/(2*pi)*(1+la)*((1+2*la)^2+4*la*(la+1)*(cos(r/2))^2)/((1+2*la)^2-4*la*(la+1)*(cos(r/2))^2)^2;

end


