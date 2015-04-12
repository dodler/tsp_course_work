tsp_data = load('data.txt');
print(tsp_data);

mx = sum(tsp_data)/n;
dx = sum(tsp_data.^2)/n-mx^2;

[cf, ncf] = cf_ncf(tsp_data,20);

rad=corel_rad(ncf);

plot_ncf(ncf);