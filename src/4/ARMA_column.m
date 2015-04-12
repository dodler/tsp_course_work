function [ outputColumn ] = ARMA_column( BettaARMA,CF,N,M)
%результат - как левая часть системы приложения В, пример модели АРСС(3,3)
outputColumn(1:N+1)=CF(1:N+1);
for k=1:N+1 %для каждого уравнения
for j=1:M %определение элемента вектор-столбца
outputColumn(k)=outputColumn(k)-BettaARMA(j)*CF(abs(j-(k-1))+1);
end;
end;
end