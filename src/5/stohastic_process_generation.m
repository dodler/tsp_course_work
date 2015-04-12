function [ Eta ] = stohastic_process_generation(Z,Alpha,Betta,N,M,MX)
%Z - случайная последовательность, распределенная по закону N(0,1).
n=length(Z); p=M; q=N;
%Z=rand(n,1);
%Z=norminv(Z,0,1);
Eta(1:n)=0;
for t=1:n-100
for j=1:p
if (t-j>0)
Eta(t)=Eta(t)+Betta(j)*Eta(t-j);
end;
end;
for j=1:q+1
%t=1. j=1 => Z(t-0), j=q+1 => Z(t-q+1)
Eta(t)=Eta(t)+Alpha(j)*Z(t-j+1+q);
end;
Eta(t)=Eta(t)+Z(t+q);
end;
Eta=Eta+MX;
end