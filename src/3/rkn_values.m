function [ Rkn ] = rkn_values(Alpha,BettaARMA,N,M)
%Rkn(k)=matrixRkn(k,:)*(Alpha');
matrixRkn=getMatrixCorrelationKsiEta(BettaARMA,N,M);
Rkn(1:N+1)=0.0;
for k=1:N+1
for j=1:N+1
Rkn(k)=Rkn(k)+matrixRkn(k,j)*Alpha(j);
end;
end;
end