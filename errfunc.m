%==========================================================================
% Description: This function calculates the RHS of the price Euler equation
% using projection on the state variables and the error between the LHS and
% the RHS.
% =========================================================================

function [ssr, error] = errfunc(coef,par,grid,gh)
% note that the second argument is ignored by optimizer

%QUESTION 1.1. COMPLETE THIS FUNCTION.

%Allocating memory
lhs = zeros(grid.size,1);
rhs = zeros(grid.size,1);
temp = zeros(gh.size,1);

%Double for loop
for i = 1:grid.size 
    %CALCULATE HERE THE LHS OF THE EULER EQUATION. HINT: YOU CAN USE PFUNC.
    lhs(i,1) = pfunc(grid.d(i), coef);
    for j = 1:gh.size
        %CALCULATE HERE THE RHS OF THE EULER EQUATION FOR EACH GAUSSIAN
        %HERMITE NODE.
        temp_d(j,1) = par.mud + par.rhod*grid.d(i) + sqrt(2)*par.sigma*gh.e(j);
        temp_c(j,1) = temp_d(j,1);
        temp_p(j,1) = pfunc(temp_d(j,1), coef);
        temp(j,1) = gh.w(j)/sqrt(pi) * ((temp_c(j,1)/grid.d(i))^(-par.gamma)*(temp_d(j,1) + temp_p(j,1)));
    end
    %CALCULATE HERE THE RHS OF THE EULER EQUATION BY COMBINING THE ELEMENTS
    %OF TEMP. (NOTE THAT THE ELEMENTS OF TEMP WILL BE OVERWRITTEN BY THE
    %CALCULATIONS FOR THE NEXT GRID POINT.)
    rhs(i,1) = par.beta*sum(temp);
end

%Sum of squared Euler equation errors
%CALCULATE HERE THE SUM OF SQUARED EULER EQUATION ERRORS

error = lhs-rhs;
ssr = sum(error.^2);
end