function [ BettaAR ] = beta_ar(CF,M)
% Из системы M уравнений находим коэффициенты b1,b2,...,bM
B=CF(2:M+1)'; %нумерация сдвинута на 1
matrixAR=zeros(M);
for i=1:M for j=1:M
matrixAR(i,j)=CF(abs(i-j)+1);
end; end;
BettaAR=matrixAR\B; %получаем вектор (b1,b2,...,bM)
end