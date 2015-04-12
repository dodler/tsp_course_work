function [ matrixRkn ] = correl_matr_xi_etha( BettaARMA, N,M )
%������� ������������� - ���������� �� �� ����� (�� ������� �).
matrixRkn=zeros(N+1);
matrixRkn(1,1)=1.0;
for k=2:N+1 %�� ���� �������, ������� �� 2-�
matrixRkn(k,k)=matrixRkn(k,k)+1.0;
for j=1:min(k-1,M) %��� ������� BattaARMA(j)
for l=1:N+1 %������ ���������� �� Alpha
matrixRkn(k,l)=matrixRkn(k,l)+BettaARMA(j)*matrixRkn(k-j,l);
end;
end;
end;
end