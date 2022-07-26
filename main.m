% =========================================================================
% LUCAS TREE ASSET PRICING MODEL: (Based on code by Wouter Den Haan)
% =========================================================================

%% Cleaning
clear all; close all; clc

%% Parameters:
par.beta  = 0.9;
par.gamma = 3.0;
par.mud   = 0.1;
par.rhod  = 0.9;
par.sigma = 0.1;

%% Grid

%Construct grid of the state variable d
d_mean = par.mud/(1-par.rhod);                      %AR(1) mean
d_std  = par.sigma/sqrt((1-par.rhod^2));            %AR(1) std. dev.
d_low  = d_mean-3*d_std;
d_high = d_mean+3*d_std;
d_numb = 6;                                         %Number of grid points
d_step = (d_high-d_low)/(d_numb-1);
d_grid = (d_low:d_step:d_high)';

%Put in struct
grid.d    = d_grid;
grid.size = d_numb;

%% Gaussian Hermite
% Set up the Gaussian hermite grid for n = 5 (using standardised code) 
gh.size = 5; [gh.e,gh.w] = hernodes(gh.size); %gh.e are the GH nodes and gh.w are the GH weights.

%% Projection for price

%Order
order = 5; %If you change the order, you might also want to change the number of grid points (line 29 above).

%Initials
initials = 'fresh'; %Type 'fresh' for fresh initial values or 'previous' to load initial values from init.mat.
switch initials
    case 'fresh'
        init = zeros(order+1,1); init(order+1,1) = par.beta/(1-par.beta); %Analytical solution for log utility.
    case 'previous'
        load init
end

%Minimization routine
routine = 'fminsearch'; %Type 'fminsearch' to use fminsearch or 'csminwel' to use csminwel.
switch routine
    case 'fminsearch'
        options = optimset('Display','Iter','MaxFunEvals',1E5,'MaxIter',...
                  1E5,'TolFun',1E-10,'TolX',1E-10);
        coef_p    = fminsearch(@(coef) errfunc(coef,par,grid,gh),init,options);
        coef_p    = fminsearch(@(coef) errfunc(coef,par,grid,gh),coef_p,options); 
        %Always run twice using the solution of the first attempt as initial conditions
        %for the second use.
        %The reason is that fminsearch is a "global" procedure that uses 
        %different step sizes initially!!
    case 'csminwel'
        addpath([cd,'\csminwel']) %Need to add directory with csminwel (subdirectory of current directory) to Matlab path.
        precision = 1E-10; maxiter = 1E5; init_invHes = 1E-10*eye(length(init));
        [fmin,coef_p,gmin,H,itct,fcount,retcodeh] = csminwel('errfunc',init,...
        init_invHes,[],precision,maxiter,par,grid,gh);
end

%The syntax to use fminsearch is very different from the syntax to use
%csminwel. With both routines the coefficients of the policy function are
%returned in coef. With fminsearch you can see the sum of squared Euler
%equation errors on the screen and with csminwel the sum of squared Euler
%equation errors is returned in fmin.

%Save new initials
init = coef_p; save init init

% get the price as a function of the state 
p = pfunc(grid.d, coef_p); 

% error evaluation
[~,error] = errfunc(coef_p,par,grid,gh);
figure(1)
scatter(grid.d,error)
title('Projection error'), box off
xlabel('\it{d_{t}}')

%% Plot: Price with state variable

%Plot price
figure(2), plot(grid.d,p)
title('Price'), box off
xlabel('\it{d_{t}}'), ylabel('\it{p_{t}}','Rotation',0)

%% Projection for risk premium and (expected) returns

%erpfunc(coef,par,grid,gh,k) function calculates
% the risk premium for k = 1
% the expected return for k = 2
% the riskfree rate for k = 3

x_rp     = erpfunc(coef_p,par,grid,gh,1);      
coef_rp = polyfit(grid.d,x_rp,order);      %Matlab's polyfit(x,y,n) function regresses/projects y on the nth-order polynomial in x.

x_er     = erpfunc(coef_p,par,grid,gh,2);      
coef_er = polyfit(grid.d,x_er,order);      

x_rf     = erpfunc(coef_p,par,grid,gh,3);      
coef_rf = polyfit(grid.d,x_rf,order);      

% %Plot these three objects
figure(3), plot(grid.d,[x_rp x_er x_rf])
title('Risk Premium, Expected Return, Risk-Free Rate'), box off
xlabel('\it{d_{t}}')
legend("risk premium", "expected return", "rf"); 

%% Simulation

%Settings
sim.size = 500;                                     %Number of periods

%Stochastics
randn('state',519)                                  %State of generator
shock = par.sigma*randn(sim.size,1);                %Shock to dividend

%Simulation
sim.d = zeros(sim.size,1); sim.d(1) = d_mean;
sim.x_rp = zeros(sim.size,1);
sim.x_er = zeros(sim.size,1);
sim.x_rf = zeros(sim.size,1);

for i = 2:sim.size
    sim.d(i) = par.mud+par.rhod*sim.d(i-1)+shock(i);
    sim.x_rp(i) = xfunc(sim.d(i), coef_rp);
    sim.x_er(i) = xfunc(sim.d(i), coef_er);
    sim.x_rf(i) = xfunc(sim.d(i), coef_rf);
end

% %Plot simulation
figure(4)
subplot(4,1,1), plot(1:sim.size,sim.d)
title('Dividend'), box off
xlabel('\it{t}'), ylabel('\it{d_{t}}','Rotation',0)
subplot(4,1,2), plot(1:sim.size,sim.x_rp)
title('Risk premium'), box off
xlabel('\it{t}'), ylabel('\it{x_{t}}','Rotation',0)
subplot(4,1,3), plot(1:sim.size,sim.x_er)
title('Expected Return'), box off
xlabel('\it{t}'), ylabel('\it{x_{t}}','Rotation',0)
subplot(4,1,4), plot(1:sim.size,sim.x_rf)
title('Risk Free Rate'), box off
xlabel('\it{t}'), ylabel('\it{x_{t}}','Rotation',0)