u1=rand(ceil(n/2),1);
u2=rand(ceil(n/2),1);

for i=1:fix(n/2)	
	X(i*2-1)=sqrt(-2*log(u2(i)))*cos(2*pi*u1(i));
	X(i*2)=sqrt(-2*log(u2(i)))*sin(2*pi*u1(i));
	X(i*2-1)=X(i*2-1)*sigma+a;
	X(i*2)=X(i*2)*sigma + a;
end;

if (fix(n/2)<ceil(n/2))
	X(n)=sqrt(-2*log(u2(ceil(n/2))))*cos(2*pi*u1(ceil(n/2)));
	X(n)=X(n)*sigma+a;
end

save lyan X;