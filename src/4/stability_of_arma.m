function [ boolStability ] = stability_of_arma( BettaARMA, M )
boolStability=true;
if (M==0) disp('������ ���������'); end;
if (M>3) disp('������������ ������ �� ���������� � ������ ������� ������'); else
if (abs(BettaARMA(M))>=1)
disp('��� 1. ������ �����������'); boolStability=false; else
if ((M==2)&&(abs(BettaARMA(1))>=1-BettaARMA(2)))
disp('��� 2. ������ �����������'); boolStability=false;
end;
if ((M==3)&&((abs(BettaARMA(1)+BettaARMA(3))>=1-BettaARMA(2))||(abs(BettaARMA(2)+BettaARMA(1)*BettaARMA(3))>=abs(1-BettaARMA(3)^2))))
disp('��� 2. ������ �����������'); boolStability=false;
end;
end;
end;
end