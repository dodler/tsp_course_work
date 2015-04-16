output = 'output.xlsx';

tsp_data = load('data.txt'); % загрузка данных
print_data(tsp_data); % вывод последовательности
n= length(tsp_data);

mx = sum(tsp_data)/n; % матожидание
dx = sum(tsp_data.^2)/n-mx^2; % дисперсия

[cf, ncf] = cf_ncf(tsp_data,20); % корелляционная функция и норм. кф

rad=corel_rad(ncf); % радиус корреляции

plot_ncf(ncf); % график нкф

%вывод 
xlswrite(output,{'Матожидание'}, 'f1:f1');
xlswrite(output,mx, 'g1:g1');
xlswrite(output,{'Дисперсия'}, 'f2:f2');
xlswrite(output, dx, 'g2:g2');
xlswrite(output, {'Корреляционная функция'}, 'f3:f3');
xlswrite(output, cf, 'g3:g3');
xlswrite(output, {'Радиус корреляции'}, 'f4:f4');
xlswrite(output, rad, 'g4:g4');

xlswrite(output, {'Номер НКФ', 'Значение НКФ'}, 'a1:b1');
xlswrite(output,(1:20)','a2:a21');
xlswrite(output, ncf(1:20)', 'b2:b21');
% конец задания 1

beta_ar1=beta_ar(cf, 1);
beta_ar2=beta_ar(cf, 2);
beta_ar3=beta_ar(cf, 3); % нашли бета 1,2,3

xlswrite(output,{'Задание 2'},'a22:a22');
xlswrite(output,{'Коэффициенты скользящего среднего'},'a23:a23');
xlswrite(output,{'1','2','3'},'b22:d22');
xlswrite(output,{beta_ar1},'b23:b23');
xlswrite(output,beta_ar2(1:2),'c23:c24');
xlswrite(output,beta_ar3(1:3),'d23:d25');

alpha_0_ar1=alpha_0_ar(beta_ar1,cf,1);
alpha_0_ar2=alpha_0_ar(beta_ar2,cf,2);
alpha_0_ar3=alpha_0_ar(beta_ar3,cf,3);

xlswrite(output,{'Первое приближение АР'},'a26:a26');
xlswrite(output,{'1','2','3',},'b26:d26');
xlswrite(output,{alpha_0_ar1},'b27:b27');
xlswrite(output,{alpha_0_ar2},'c27:c27');
xlswrite(output,{alpha_0_ar3},'d27:d27');

theor_ncf(1,:)=theor_ncf_ar(beta_ar1,alpha_0_ar1_z,1);
theor_ncf(2,:)=theor_ncf_ar(beta_ar2,alpha_0_ar2_z,2);
theor_ncf(3,:)=theor_ncf_ar(beta_ar3,alpha_0_ar3_z,3);

xlswrite(output, {'теоретические нкф'}, 'h1:h1');
xlswrite(output,theor_ncf(1,:),'h2:h21');
xlswrite(output,theor_ncf(2,:),'i2:i21');
xlswrite(output,theor_ncf(3,:),'j2:j21');


% графическое отображение моделей СП
osx=0:10;
plot(osx,ncf(1:11),'k');
hold on;
grid on;
plot(osx,theor_ncf(1,1:11),'r');
plot(osx,theor_ncf(2,1:11),'g');
plot(osx,theor_ncf(3,1:11),'b');
xlabel('Номер НКФ, \tau');
ylabel('Значение НКФ, r_\eta (\tau)');
legend('Выборочные НКФ','НКФ модели АР(1)','НКФ модели АР(2)','НКФ модели АР(3)');
hold off;
%оценка погрешностей
ar_error=zeros(3,20);
for i=1:3
ar_error(i)=ncf_error(theor_ncf(i,:),ncf(:));
end;

% задание 3 из второг задания можно убрать приблжение альфа порядка 0 при
% помощи вспомогательной последовательности
% 
% N=0;
% alpha0_ma=sqrt(cf(1));
% theor_f_ncf0=ncf_ma(alpha0);
% theor_ncf0(1)=theor_f_ncf0(1);
% for i=2:11
% 	theor_ncf0(i)=0;
% end;
% existance_of_arma_solution(cf,N);
% % добавить вывод в xlsx
% 
% N=1;
% for i=1:N+1
% 	alpha01_ma(i)=sqrt(cf(1)/N+1)-0.5*i;
% end;
% [alpha1, fval1] = fsolve(@getSystemOfMA, alpha0_ma, options, cf, N);
% theor_f_ncf1=ncf_ma(alpha1);
% theor_ncf(1,1:N+1)=theor_f_ncf1(1:N+1);
% theor_ncf(1,N+2:11)=0;
% existance_of_arma_solution(cf,N);
% % лог
% 
% N=2;
% alpha_02_ma(1)=sqrt(cf(1));
% alpha_02_ma(2:N+1)=0;
% [alpha2,fval2]=fsolve(@getSystemOfMA, alpha_02_ma,options, cf,N);
% theor_f_ncf2=ncf_ma(alpha2);
% theor_ncf(2, 1:N+1)=theor_f_ncf2(1:N+1);
% theor_ncf(2,2+N:20)=0;
% existance_of_arma_solution(cf,N);
% 
% N=3;
% alpha3(1)=sqrt(cf(1));
% alpha3(2:N+1)=0;
% [alpha3,fval3]=fsolve(@getSystemOfMA, alpha3, options,cf,N);
% theor_f_ncf3=ncf_ma(alpha3);
% theor_ncf(3,1:N+1)=theor_f_ncf3(1:N+1);
% theor_ncf(3,N+2:11)=0;
% 
% 
% ox_=0:10;
% plot(os_, ncf(1:11),'k');
% hold on;  
% grid on;  
% plot(newOsX,tvNCF0(1:11),'y');  
% plot(newOsX,tvNCF(1,1:11),'r');  
% plot(newOsX,tvNCF(2,1:11),'g');  
% %plot(newOsX,tvNCF(3,1:11),'b');
% xlabel('Номер НКФ, \tau');  
% ylabel('Значение НКФ, r_\eta (\tau)');  
% legend('Выборочные НКФ','НКФ модели СС(0)','НКФ модели СС(1)','НКФ модели СС(2)');  
% hold off;   
% 
% ma0_error=error(theor_f_ncf0,ncf);
% for i=1:3
% 	ma_error(i)=error(theor_ncf(i,:),ncf);
% end;
% 
% % конец задания 3 сделать вывод данных в файл
% 
% for N=1:3
% 	for M=1:3
% 		beta_arma(N,M,:)=zeros(1,3);
% 		alpha_arma(N,M,:)=zeros(1,4);
% 		theor_ncf_arma(N,M,:)=zeros(1,20);
% 		cf_zeta(N,M,:)=zeros(1,7);
% 	end;
% end;
% 
% error_arma=zeros(3);
% 
% for N=1:3
% 	for M=1:3
% 		beta(N,M,1:M)=beta_arma(cf,N,M);
% 		stability(N,M)=stability_of_arma(beta(N,M,1:M),M);
% 		solution(N,M)=existance_of_arma_solution(beta(N,M,1:M),cf,N,M,10);
% 		if (stability(N,M) && solution(N,M) == true)
% 			alpha_arma(N,M,1)=sqrt(cf(1));
% 			alpha_arma(N,M,2:N+1)=0;
% 			alpha_arma(N,M,1:N+1)=fsolve(@getSystemOfMA, alpha_arma(N,M,1:N+1),options,beta_arma(N,M,1:M),cf,N,M);
% 		end;
% 	end;
% end;
% 
% % надо ли второй метод?
% 
% for N=1:3
% 	for M=1:3
% 		if (stability(N,M) && solution(N,M)==true)
% 			theor_ncf_arma(N,M,:)=arma_tcf(beta_arma(N,M,:),alpha_arma(N,M,:),N,M);
% 			arma_error(N,M)=error_arma(theor_ncf_arma(N,M,:),ncf);
% 		end;
% 	end;
% end;
% 
% for N=1:3
% 	for M=1:3
% 		if (stability(N,M)&&solution(N,M)==true)
% 			sma2(N,M)=sum(arma_system(alpha_arma(N,M,:),beta_arma(N,M,:),cf,N,M));
% 			
% 		end;
% 	end;
% end;
% 
% % коне цзадания 6
% X = load('lyan.txt');
% 
% N=0;M=3;
% etha2=stohastic_process_generation(X, alpha0AR3, BettaAR3, N,M,mx);
% etha03=etha2(901, 5900);
% 
% N=2;M=0
% etha2=stohastic_process_generation(X, alpha2, 0,N,M,mx);
% etha20=etha2(901, 5900); % возможно тут по варикам
% 
% print_data(eta);
% n=length(eta);
% m_etha=sum(eta)/n;
% d_etha=(eta.^2)/n-m_etha^2;