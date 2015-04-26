clear();
SOURCE_FILE_NAME = 'C:\sample.txt'; // входные данные(выборка)
FLOAT_FORMAT = '%16.4f'; // формат представления чисел с плавающей точкой
INT_FORMAT = '%d'; // формат представления целых чисел
EPSILON = 1.0E-6; // точность
MAX_AR_LEVEL = 3;
MAX_MA_LEVEL = 2;
IMITATION_LENGTH = 5000;

// ввод данных

x = fscanfMat(SOURCE_FILE_NAME); // выборка
n = length(x); // объем выборки

// вывод матрицы 

function printMat(M, mformat),
  [n, m] = size(M);
  for i = 1:n,
    for j = 1:m,
      printf(mformat + " ", M(i, j));
    end;
    printf("\n");
  end;
endfunction;

// график СП 
//x=fscanfMat(SOURCE_FILE_NAME);
//f1=mean(x);
//f2=cmoment(x,2,1);
//f3=sqrt(f2);
//plot(x(1:120,:),'k-');
//plot(0:120:120,f1,'r-');
//plot(0:120:120,f1-f3,'b-');
//plot(0:120:120,f1+f3,'b-');
//xtitle('','Индекс','Значения последовательности');
//legend('Исходный процесс','Среднее','СКО',4,%f);
//xgrid(111);



// ЗАДАНИЕ 1 Оценка корреляционной функции 

function R = correlation(k, x)
  if (k < 0),
    k = -k;
  end;
  R = 0.0;
  n = length(x);
  meanx = mean(x);
  for i = 1:(n-k),
    R = R + (x(i) - meanx) * (x(i+k) - meanx);
  end;
  R = R / (n - k - 1);
endfunction;

// нормированная корреляционная функция

function r = ncorrelation(k, x)
  r = correlation(k, x) / correlation(0, x);
endfunction;

// радиус корреляции

function T = corrdist(x)
  coefficient = 0.01;
  T = coefficient * length(x) - 2;
  em1 = exp(-1);
  while (T >= 0) & (abs(ncorrelation(T, x)) < em1),
    T = T - 1;
  end;
  T = T + 1;
endfunction;

// график нормированной корреляционной функции

function corrplot()
  m = 10;
  t = [0:m];
  y = zeros(length(t), 1);
  for i = 1:length(t),
    y(i) = ncorrelation(t(i), x);
  end;
  xgrid();
xgrid();
  
plot2d(t, y, 11);
  xtitle('Корелляция', 'Номер индекса', 'Нормализованная корелляционная функция');
  a = gca();
  t = [0, 10];
  plot2d(t, [1/2.71828, 1/2.71828], style=color("green"))
  plot2d(t, [1/-2.71828, 1/-2.71828], style=color("green"))
endfunction;


// Main
meanx = mean(x); // выборочное среднее
svx = variance(x); // выборочная дисперсия

m = 10; // R = zeros(m, 1);
r = zeros(m, 1);

for i = 0:m,
  R(i+1) = correlation(i, x);
  r(i+1) = ncorrelation(i, x);
end;

Tcorr = corrdist(x);
printf("Выборочное среднее: " + FLOAT_FORMAT + "\n", meanx);
printf("Выборочная дисперсия: " + FLOAT_FORMAT + "\n", svx);
printf("Оценка корелляционной функции:\n");
printMat(R, FLOAT_FORMAT);
printf("Оценка нормальзованной кореллиционной функции:\n");
printMat(r, FLOAT_FORMAT);
printf("Радиус корелляции: " + INT_FORMAT + "\n", Tcorr);

scf(1);
//corrplot();


// Задание 2

// Поиск коэффициентов бета авторегрессии

function betas = ar(x, arLevel, maLevel)
    if arLevel=0 then arLevel=1
    end
  R = zeros(2 * arLevel, 1);
  for i = (maLevel - arLevel + 1) : (maLevel + arLevel),
    R(i - maLevel + arLevel) = correlation(i, x);
  end;
  Rmm = zeros(arLevel, arLevel);
  _n = maLevel;
  for i = 1:arLevel, // ошибка при arLevel=0
    _m = _n;
    for j = 1:arLevel,
      Rmm(i, j) = R(_m - maLevel + arLevel);
      _m = _m - 1;
    end;
    _n = _n + 1;
  end;
 Rm = -R(arLevel + 1 : 2 * arLevel);
  
  
  betas = linsolve(Rmm, Rm);
endfunction;

// Корреляционная функция

function rrr = mcorrelation(k, alph, betas)
  rrr = alph(k+1);
  len = min(k, length(betas));
  for j = 1 : len,
    rrr = rrr + betas(j) * mcorrelation(k - j, alph);
  end;
endfunction;

// Поиск коэффициентов скользящего среднего

function alphas = ma(x, arLevel, maLevel, betas)
  for i = 0 : max([arLevel, maLevel]),
    R(i+1) = correlation(i, x);
  end;
  function zr = syst(alph)
    for k = 0 : maLevel,
      zr(k+1) = -R(k+1);
      for i = k : maLevel,
        zr(k+1) = zr(k+1) + alph(i+1) * mcorrelation(i - k, alph, betas);
      end;
      for j = 1 : arLevel,
        zr(k+1) = zr(k+1) + betas(j) * R(abs(k - j) + 1);
      end;
    end;
  endfunction;
  [alphas, values, info] = fsolve([1 : (maLevel+1)], syst);
  for i = 1:length(values),
    if (abs(values(i)) > EPSILON | info == 4) then
      alphas(1) = %i;
      break;
    end;
  end;
endfunction;

// Вектор

function s = image(v)
  s = %F;
  for i = 1:length(v),
    if (imag(v(i)) <> 0) then
      s = %T;
      break;
    end;
  end;
endfunction;

// стабильность модели

function s = stable(betas)
  p = poly([pertrans(-betas) 1], "z", "coeff");
  z = roots(p);
  s = %T;
  for i = 1:length(z),
    if (abs(z(i)) >= 1) then
      s = %F;
      break;
    end;
  end;
endfunction;

// Main
alphas_list = list();
betas_list = list();
for i = 0 : MAX_AR_LEVEL,
  for j = 0 : MAX_MA_LEVEL,
    betas = ar(x, i, j);
    alphas = ma(x, i, j, betas);
    alphas_list($+1) = alphas;
    betas_list($+1) = betas;
    printf("ARMA(" + INT_FORMAT + "," + INT_FORMAT + ")\n", i, j);
    if (image(alphas)) then
      printf("Model does not exist.\n");
      continue;
    end;
    if (~stable(betas)) then
      printf("Model exists, but not stable.\n");
      continue;
    end;
    printf("alpha:\n")
    printMat(alphas, FLOAT_FORMAT);
    printf("beta:\n");
    printMat(betas, FLOAT_FORMAT);
  end;
end;


//Задание3

// Теорет. корреляция для АРМА

function R = theoretical_corr(betas, startR, k)
  nm = length(startR) - 1;
  k = abs(k);
  if (k > nm) then
     R = 0;
     M = length(betas);

     for j = 1 : M,
         R = R + betas(j) * theoretical_corr(betas, startR, k - j);
    end;
  else
    R = startR(k + 1);
  end;
endfunction;

// Норм. теор. корр. для АРМА

function r = norm_theoretical_corr(betas, startR, k)
    r = theoretical_corr(betas, startR, k) / theoretical_corr(betas, startR, 0);
endfunction;

// Квадратичное отклонение

function epsilon = quadratic_error(x, y)
   epsilon = 0;
   m = min(length(x), length(y));
   for j = 1 : m,
       epsilon = epsilon + (x(j) - y(j))^2;
   end;
endfunction;

function corrplot2(betas, R, N, M, m_)
   p.thickness = 6;
   m = m_;
   t = [-m:m];
   y = zeros(length(t), 1);
   for i = 1:length(t),
        y(i) = norm_theoretical_corr(betas, R(1 : N + M + 1), t(i));
   end;

   plot2d3(t, y, axesflag=5, style=2);

   a = gca();
   p = a.children.children;
   p.thickness = 3;
   p.mark_mode = "on";
   p.mark_size_unit = "point";
   p.mark_style = 11;
   p.mark_size = 3;
endfunction;
function eta = imitate(alphas, betas, meanx, count)
    defect = 1000;
    eta = zeros(count + defect + 1, 1);
    ksi = grand(count + defect + 1, 1, 'nor', 0, 1);
    N = length(alphas) - 1;
    M = length(betas);
    for k = 1 : count + defect + 1,
        eta(k) = 0;
        for i = 0 : N,
            if (k - i > 0) then
               eta(k) = eta(k) + alphas(i+1) * ksi(k - i);
            end;
         end;

       for j = 1 : M,
            if (k - j > 0) then
                 eta(k) = eta(k) + betas(j) * eta(k - j);
           end;
        end;
    end;
 
    eta = eta(defect + 2 : count + defect + 1) + meanx;
endfunction;

m = 10; // Analysis Depth

epsilon = zeros(MAX_AR_LEVEL, MAX_MA_LEVEL);
best_ar_eps = %inf;
best_ma_eps = %inf;
best_arma_eps = %inf;

best_ar_alpha = [];
best_ma_alpha = [];
best_arma_alpha = [];
best_ar_beta = [];
best_arma_beta = [];

R = zeros(m+1, 1);
r = zeros(m+1, 1);

for k = 0 : m,
    R(k+1) = correlation(k, x);
    r(k+1) = R(k+1) / R(1);
end;

r_model = zeros(m+1, 1);

for i = 0 : MAX_AR_LEVEL,
   for j = 0 : MAX_MA_LEVEL,
       alphas = alphas_list(i * (MAX_MA_LEVEL + 1) + j + 1);
       betas = betas_list(i * (MAX_MA_LEVEL + 1) + j + 1);
       
       for k = 0 : m,
           r_model(k+1) = norm_theoretical_corr(betas, R(1 : (i + j + 1)), k);
       end;

       epsilon(i+1, j+1) = quadratic_error(r, r_model);

       if (i == 0) & (j~=3)&(j~=2)&(j~=1)&(epsilon(1, j+1) < best_ma_eps) then
           best_ma_eps = epsilon(1, j+1);
           best_ma_alpha = alphas;
       elseif (j == 0) & (epsilon(i+1, 1) < best_ar_eps) then
               best_ar_eps = epsilon(i+1, 1);
               best_ar_alpha = alphas;
               best_ar_beta = betas;
        elseif (epsilon(i+1, j+1) < best_arma_eps) then
                best_arma_eps = epsilon(i+1, j+1);
                best_arma_alpha = alphas;
                best_arma_beta = betas; 
        end;
    end;
  end;

printf("Epsilon:\n");
printMat(epsilon, '%16.7f');
printf("Best models:\nAR(" + INT_FORMAT + "), MA(" + INT_FORMAT + "), ARMA(" + INT_FORMAT + "," + INT_FORMAT + ").\n", length(best_ar_beta), length(best_ma_alpha) - 1, length(best_arma_beta), length(best_arma_alpha) - 1);



