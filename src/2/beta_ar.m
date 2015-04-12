function [ BettaAR ] = beta_ar(CF,M)
% �� ������� M ��������� ������� ������������ b1,b2,...,bM
B=CF(2:M+1)'; %��������� �������� �� 1
matrixAR=zeros(M);
for i=1:M for j=1:M
matrixAR(i,j)=CF(abs(i-j)+1);
end; end;
BettaAR=matrixAR\B; %�������� ������ (b1,b2,...,bM)
end