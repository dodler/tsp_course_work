function [ tNCF ] = arma_tcf(BettaARMA,AlphaARMA,N,M )
% матрица разложения КФ через коэффициенты BettaARMA
matrixARMA=zeros(max(M,N+1));
for i=1:max(N+1,M+1)
for j=1:max(N+1,M+1)
% -1 b1 b2 b3 .. bM
% b1 b2-1 b3 ... bM-1 0
% b2 b1+b3 b4-1 ... 0 0
% b3 b4+b2 b5+b1 b6-1 0 0 0
if ((i==j)&&(i~=1))
if (i+j-2<=M)
matrixARMA(i,j)=BettaARMA(i+j-2);
end;
end;
if (j>i)
if (j+i-2<=M)
matrixARMA(i,j)=BettaARMA(j+i-2);
end;
end;
if ((j<i)&&(j~=1))
if (i-j>0)
matrixARMA(i,j)=BettaARMA(i-j);
end;
if (j+i-2<=M)
matrixARMA(i,j)=matrixARMA(i,j)+BettaARMA(i+j-2);
end;
end;
if ((j<i)&&(j==1))
matrixARMA(i,j)=BettaARMA(i-1);
end;
end;
matrixARMA(i,i)=matrixARMA(i,i)-1;
end;
%нахождение разложения "правой"(с коэффициентами Alpha) части системы А
Rkn = getValuesRkn(AlphaARMA,BettaARMA,N,M);
rigthVectorARMA(1:N+1)=0.0;
for k=1:N+1
for j=1:N-k+2
rigthVectorARMA(k)=rigthVectorARMA(k)+AlphaARMA(k+j-1)*Rkn(j);
end;
end;
if (N+1<M+1) rigthVectorARMA(N+2:M+1)=0.0; end;
%решение матричного уравнения
teorNCF=matrixARMA\(-rigthVectorARMA)';
%disp(teorNCF);
tNCF(1:max(N+1,M+1))=teorNCF(1:max(N+1,M+1));
for i=max(N+2,M+2):20
tNCF(i)=0.0;
for j=1:M
if (i>j)
tNCF(i)=tNCF(i)+BettaARMA(j)*tNCF(i-j);
end;
end;
end;
tNCF=tNCF/tNCF(1);
end