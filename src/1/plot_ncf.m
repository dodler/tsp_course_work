function [] = plot_ncf( ncf )
%PLOT_NCF функция строит график зависимости нкф от времени
n=length(ncf);
plot(0:n-1,ncf,'r');
hold on;
grid on;
plot([0 n],[exp(-1) exp(-1)],'g')
ylabel('НКФ, r_\eta (\tau)');
xlabel('Разность времени, или номер, \tau');
legend('Выборочные НКФ','exp(-1)');
hold off;

end

