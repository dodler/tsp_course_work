function [ matrixRkn ] = correl_matr_xi_etha( BettaARMA, N,M )
%матрица коэффициентов - разложения КФ по Альфа (из система Г).
matrixRkn=zeros(N+1);
matrixRkn(1,1)=1.0;
for k=2:N+1 %по всем строкам, начиная со 2-й
matrixRkn(k,k)=matrixRkn(k,k)+1.0;
for j=1:min(k-1,M) %для каждого BattaARMA(j)
for l=1:N+1 %делать разложение по Alpha
matrixRkn(k,l)=matrixRkn(k,l)+BettaARMA(j)*matrixRkn(k-j,l);
end;
end;
end;
end