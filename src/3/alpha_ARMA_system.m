function [ systemOfARMA ] = alpha_ARMA_system(AlphaARMA,BettaARMA,CF,N,M )
% система используется для нахождения параметров AlphaARMA
Rkn=getValuesRkn(AlphaARMA,BettaARMA,N,M); %получение значений КФ(xi,eta)
systemOfARMA=-getColumnForARMA(BettaARMA,CF,N,M); %и столбца левой части системы приложения В (пример модели АРСС(3,3))
for k=1:N+1
for j=k:N+1
systemOfARMA(k)=systemOfARMA(k)+AlphaARMA(j)*Rkn(j-k+1); %суммирование с правой частью
end;
end;
end