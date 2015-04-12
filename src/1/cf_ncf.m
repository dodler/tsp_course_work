function [ CF, NCF ] = cf_ncf(input,NUMBER_CF)
%возвращает значения НКФ и КФ от выборки input, в том числе - нестационарной
n=length(input);
MX=sum(input)/n;
DX=sum(input.^2)/n-MX^2;
CF(1:NUMBER_CF)=0.0; NCF(1:NUMBER_CF)=0.0;
CF(1)=DX; NCF(1)=1.0;
x1=0.0; x2=0.0; %MX
y1=0.0; y2=0.0; %DX
for k=2:NUMBER_CF
x1=MX*n; x2=x1; y1=(DX+MX^2)*n; y2=y1;
for j=n-k+2:n x1=x1-input(j); end; x1=x1/(n-k+1); %MX1
for j=1:k-1 x2=x2-input(j); end; x2=x2/(n-k+1); %MX2
for j=n-k+2:n y1=y1-input(j)^2; end; y1=y1/(n-k+1)-x1^2; %DX1
for j=1:k-1 y2=y2-input(j)^2; end; y2=y2/(n-k+1)-x2^2; %DX2

%кф
for j=1:n-k+1
CF(k)=CF(k)+(input(j)-x1)*(input(j+k-1)-x2);
end;
CF(k)=CF(k)/(n-k+1);
% нкф
NCF(k)=CF(k)/sqrt(y1*y2);
end;
end