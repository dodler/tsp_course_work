%Alpha0
function [ alpha0AR ] = alpha_0_ar( BettaAR,CF,M)
alpha0AR=CF(1);
for i=1:M
alpha0AR=alpha0AR-BettaAR(i)*CF(i+1);
end;
alpha0AR=sqrt(alpha0AR);
end