function [ outputColumn ] = ARMA_column( BettaARMA,CF,N,M)
%��������� - ��� ����� ����� ������� ���������� �, ������ ������ ����(3,3)
outputColumn(1:N+1)=CF(1:N+1);
for k=1:N+1 %��� ������� ���������
for j=1:M %����������� �������� ������-�������
outputColumn(k)=outputColumn(k)-BettaARMA(j)*CF(abs(j-(k-1))+1);
end;
end;
end