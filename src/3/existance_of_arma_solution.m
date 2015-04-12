function [ boolExistence, Alpha ] = existance_of_arma_solution(BettaARMA,CF,N,M,iteration)
boolExistence=true;
%����� ����� ������� ����(3,3) ���������� �
outputColumn = getColumnForARMA(BettaARMA,CF,N,M);
%��������� �����������
Alpha(1)=sqrt(CF(1)); Alpha(2:N+1)=0.0;
cc=0.0; %����� �������. (��� ����)
for l=1:iteration
%���������� R_ksi,eta �� Alpha,Betta
Rkn = getValuesRkn(Alpha, BettaARMA,N,M);
% ����� ������� �������� ���������� � ���������� �������� �
% ������������� ������:
k=N+1;
while k>=2
%������������ ������� � ��� �����������
Alpha(k)=outputColumn(k);
for j=k+1:N+1
Alpha(k)=Alpha(k)-Alpha(j)*Rkn(j-k+1);
end;
Alpha(k)=Alpha(k)/Alpha(1);
Rkn = getValuesRkn(Alpha, BettaARMA,N,M);
k=k-1;
end;
x=outputColumn(1);
for j=2:N+1 x=x-Alpha(j)*Rkn(j); end;
if (x<0)
disp('������� ���'); boolExistence=false; break;
else Alpha(1)=sqrt(x);
end;
cc=sqrt(sum((getSystemOfARMA(Alpha,BettaARMA,CF,N,M)).^2));
end;
if (boolExistence) disp('������� ����'); end;
%disp([N,M]); disp('����� �������:');
%disp(cc);
end