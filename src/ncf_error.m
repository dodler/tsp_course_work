function [x] ncf_error(theor_ncf, ncf)
	for i=2:11
		x(i-1)=theor_ncf(i)-ncf(i);
	end;
	e=sum(x.^2);
end;