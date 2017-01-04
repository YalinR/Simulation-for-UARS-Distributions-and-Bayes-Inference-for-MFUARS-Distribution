function [g]=posteriordensityMFUARS(k,S,n,OS)
sumo=SUMO(n,OS);
g=(exp(k*(trace(S.'*sumo)-n)))*sqrt(2/k*(besseli(0,2*k))^2-2/(k^2)*besseli(0,2*k)*besseli(1,2*k)+(1/(k^2)-2/k)*(besseli(1,2*k))^2)/((besseli(0,2*k)-besseli(1,2*k))^(n+1));
end