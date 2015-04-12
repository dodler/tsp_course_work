function [ vectorNCF ] = ncf_ma( Alpha )
% Finding the values of the NCF
sizeAlpha=length(Alpha);
for k=1:sizeAlpha
vectorNCF(k)=0.0;
for j=1:sizeAlpha-k+1
vectorNCF(k)=vectorNCF(k)+Alpha(j)*Alpha(j+k-1);
end;
end;
vectorNCF=vectorNCF/vectorNCF(1);
end