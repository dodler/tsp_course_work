function [ flag ] = solution_exist( CF,N )
%�������� �� ������������� �������
flag=true; %������� ���
for j=1:N+1 koef(j)=0.0; end;
approxAlpha(1)=sqrt(CF(1)); for i=2:N+1 approxAlpha(i)=0.0; end;
while (flag)
k=N+1;
while (k>=2)
koef(k)=0.0;
for j=2:N+1-(k-1)
koef(k)=koef(k)+approxAlpha(j)*approxAlpha((k-1)+j);
end;
approxAlpha(k)=(CF(k)-koef(k))/approxAlpha(1);
k=k-1;
end;
ppp=CF(1)-sum(approxAlpha.^2)+approxAlpha(1)^2;
if (ppp<0) disp('������� ���'); break; else
approxAlpha(1)=sqrt(ppp);
end;
flag=(sqrt(sum(getSystemOfMA(approxAlpha,CF,N).^2))>0.0001); %����� ������� ������ 0.0001
end;
flag=~flag;
if (~flag) disp('����������� (��������� ���������� ��������) ');
disp(sqrt(sum(getSystemOfMA(approxAlpha,CF,N).^2)));
end;
end