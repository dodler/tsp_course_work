function [ boolStability ] = stability_of_arma( BettaARMA, M )
boolStability=true;
if (M==0) disp('Модель устойчива'); end;
if (M>3) disp('Устойчивость модели не определена в рамках данного метода'); else
if (abs(BettaARMA(M))>=1)
disp('Шаг 1. Модель неустойчива'); boolStability=false; else
if ((M==2)&&(abs(BettaARMA(1))>=1-BettaARMA(2)))
disp('Шаг 2. Модель неустойчива'); boolStability=false;
end;
if ((M==3)&&((abs(BettaARMA(1)+BettaARMA(3))>=1-BettaARMA(2))||(abs(BettaARMA(2)+BettaARMA(1)*BettaARMA(3))>=abs(1-BettaARMA(3)^2))))
disp('Шаг 2. Модель неустойчива'); boolStability=false;
end;
end;
end;
end