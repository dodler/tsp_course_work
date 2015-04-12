function [ tNCF ] = theor_ncf_ar(BettaAR,alpha0AR,M )
% NCF(i+1)=NCF(i)
% |-1 b1 b2 b3 .. bM | | CF(1)|= -alpha0^2
% |b1 b2-1 b3 ... bM-1 0 | | CF(2)|= 0
% |b2 b1+b3 b4-1 ... 0 0 | | CF(3)|= 0
% |b3 b4+b2 b5+b1 b6-1 0 0 0 | | CF(4)|= 0
matrixAR=zeros(M+1);
for i=1:M+1
for j=1:M+1
% if ((j>i)&&(i+j-2<=M)) matrixAR(i,j)=BettaAR(i+j-2); end;
% if ((i==j)&&(i>1)&&(i+j-2<=M)) matrixAR(i,j)=BettaAR(i+j-2); end;
% if ((j==1)&&(i>j)) matrixAR(i,j)=BettaAR(i+j-2); end;
% if ((j>1)&&(i>j)) matrixAR(i,j)=BettaAR(i+j-2)+BettaAR(i-j); end;
% =>
if (((j>i)&&(i+j-2<=M))||((j==1)&&(i>j))||((j>1)&&(i>j)&&(i+j-2<=M))||((i==j)&&(i>1)&&(i+j-2<=M)))
matrixAR(i,j)=BettaAR(i+j-2); end;
if ((j>1)&&(i>j)) matrixAR(i,j)=matrixAR(i,j)+BettaAR(i-j); end;
end;
matrixAR(i,i)=matrixAR(i,i)-1.0;
end;
% disp(matrixAR);
rigthVectorAR(1)=-alpha0AR^2; rigthVectorAR(2:M+1)=0.0;
teorNCF=matrixAR\rigthVectorAR';
tNCF(1:(M+1))=teorNCF(1:(M+1));
for i=(M+2):20
tNCF(i)=0.0;
for j=1:M
tNCF(i)=tNCF(i)+BettaAR(j)*tNCF(i-j);
end;
end;
%disp('Корреляционные функции ');
% disp(tNCF);
tNCF=tNCF/tNCF(1);
end