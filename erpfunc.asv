function erp = erpfunc(coef,par,grid,gh,k)

%Allocating memory
p = zeros(grid.size,1);
er = zeros(grid.size,1);
rf = zeros(grid.size,1);
temp = zeros(gh.size,1);
temp_d = zeros(gh.size,1);
temp_p = zeros(gh.size,1);
temp_f = zeros(gh.size,1);

%Expected rate of return
for i = 1:grid.size
    p(i,1) = pfunc(grid.d(i),coef);
    for j = 1:gh.size
        temp_d(j,1) = par.mud + par.rhod*grid.d(i) + sqrt(2)*par.sigma*gh.e(j);
        temp_p(j,1) = pfunc(temp_d(j,1), coef);
        temp(j,1) = gh.w(j)/sqrt(pi)*(temp_d(j,1) + temp_p(j,1));

        temp_f(j,1) = gh.w(j)/sqrt(pi)*(temp_d(j,1)/grid.d(i))^(-par.gamma);
    end
    er(i,1) = sum(temp)/p(i,1)-1;
    rf(i,1) = 1/(par.beta*sum(temp_f)) - 1;
end

%Risk-free rate
if k == 1
    erp = er-rf;
end
if k == 2
    erp = er;
end
if k == 3
    erp = rf;
end


end