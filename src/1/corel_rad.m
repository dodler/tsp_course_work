function [ correlationRadius ] = corel_rad(NCF)
%������� ��� ������� ������� ���������� �� ��������� ���
i=length(NCF);
if (abs(NCF(i))>exp(-1))
disp('���������� ������ �������� ���')
end;
while abs(NCF(i))<exp(-1)
i=i-1;
end;
correlationRadius=i; %�.�. ��� ���. ����� �� ���������
disp('������ ����������'); disp(correlationRadius);
end