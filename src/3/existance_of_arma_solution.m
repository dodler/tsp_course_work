function [ boolExistence, Alpha ] = existance_of_arma_solution(BettaARMA,CF,N,M,iteration)
boolExistence=true;
%левая часть системы АРСС(3,3) приложения В
outputColumn = getColumnForARMA(BettaARMA,CF,N,M);
%начальное приближение
Alpha(1)=sqrt(CF(1)); Alpha(2:N+1)=0.0;
cc=0.0; %норма функций. (для себя)
for l=1:iteration
%Разложение R_ksi,eta по Alpha,Betta
Rkn = getValuesRkn(Alpha, BettaARMA,N,M);
% метод простых итераций начинается с последнего элемента и
% заканчивается первым:
k=N+1;
while k>=2
%приравниваем элемент к его приближению
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
disp('Решения нет'); boolExistence=false; break;
else Alpha(1)=sqrt(x);
end;
cc=sqrt(sum((getSystemOfARMA(Alpha,BettaARMA,CF,N,M)).^2));
end;
if (boolExistence) disp('Решение есть'); end;
%disp([N,M]); disp('Норма функций:');
%disp(cc);
end