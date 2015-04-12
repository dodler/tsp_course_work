function [ systemOfARMA ] = alpha_ARMA_system(AlphaARMA,BettaARMA,CF,N,M )
% ������� ������������ ��� ���������� ���������� AlphaARMA
Rkn=getValuesRkn(AlphaARMA,BettaARMA,N,M); %��������� �������� ��(xi,eta)
systemOfARMA=-getColumnForARMA(BettaARMA,CF,N,M); %� ������� ����� ����� ������� ���������� � (������ ������ ����(3,3))
for k=1:N+1
for j=k:N+1
systemOfARMA(k)=systemOfARMA(k)+AlphaARMA(j)*Rkn(j-k+1); %������������ � ������ ������
end;
end;
end