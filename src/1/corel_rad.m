function [ correlationRadius ] = corel_rad(NCF)
%функция для расчета радиуса корреляции по известной нкф
i=length(NCF);
if (abs(NCF(i))>exp(-1))
disp('Необходимо больше значений НКФ')
end;
while abs(NCF(i))<exp(-1)
i=i-1;
end;
correlationRadius=i; %т.к. еще сущ. сдвиг по нумерации
disp('Радиус корреляции'); disp(correlationRadius);
end