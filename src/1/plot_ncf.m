function [] = plot_ncf( ncf )
%PLOT_NCF ������� ������ ������ ����������� ��� �� �������
n=length(ncf);
plot(0:n-1,ncf,'r');
hold on;
grid on;
plot([0 n],[exp(-1) exp(-1)],'g')
ylabel('���, r_\eta (\tau)');
xlabel('�������� �������, ��� �����, \tau');
legend('���������� ���','exp(-1)');
hold off;

end

