function [ Zeta ] = zeta_sp(initialData, Betta,MX)
%этот процесс нужен только для нахождения его КФ
M=length(Betta);
n=length(initialData);
ksi=rand(n,1);
Zeta=initialData-MX;
for i=2:M
% for j=1:i-1
% Zeta(i)=Zeta(i)-Betta(j)*initialData(i-j);
% end;
for j=1:M
if (j<=i-1)
Zeta(i)=Zeta(i)-Betta(j)*initialData(i-j);
else
Zeta(i)=Zeta(i)-Betta(j)*initialData(n+i-j);
end;
end;
end;
for i=M+1:n
for j=1:M
Zeta(i)=Zeta(i)-Betta(j)*initialData(i-j);
end;
end;
end