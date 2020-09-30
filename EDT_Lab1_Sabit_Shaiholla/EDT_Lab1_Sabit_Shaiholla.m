clc
clear all
close all
%---"1st Lab Assignment - Sabit Shaiholla 201345123"----------%
% "ALLFUNCS.m" file should be in the same folder as this file %
%----------------May take couple of minutes to run------------%

%--------------------------INPUTS-----------------------------%
obj_fun_min = @(x)(100 - x)^2; % to minimize
obj_fun_max = @(x)(100 - x)^2*(-1); %to maximize
% for methods requiring derivatives 
syms x
obj_fun_der = (100 - x)^2;
x0 = 1; % initial guess
xL = 60; % lower bound of search
xU = 150; % upper bound of search
e = 1e-04; % stopping criterion
x_act = 100; % actual min val of x
f_act = 0; % actual min val of a function
%------------------------------------------------------------%

% calling sub-functions from the ALLFUNCS class
[x_1, f_1] = ALLFUNCS.IntervalHalving(obj_fun_min,xL,xU,e);
[x_2, f_2] = ALLFUNCS.GoldenSection(obj_fun_max,xL,xU,e);
[x_3, f_3] = ALLFUNCS.QuadraticEstimation(obj_fun_min,xL,xU);
[x_4, f_4] = ALLFUNCS.NewtonRaphson(obj_fun_der,x,x0,e);
[x_5, f_5] = ALLFUNCS.Bisection(obj_fun_der,x,xL,xU,e);
[x_6, f_6] = ALLFUNCS.Secant(obj_fun_der,x,xL,xU,e);
% using built-in fmincon
[x_7, f_7] = fmincon(obj_fun_min,x0,[],[],[],[],xL,xU); 

% array of opt. func names
func_names = ["IntervalHalving","GoldenSection"...
    "QuadraticEstimation","NewtonRaphson"...
    "Bisection","Secant","Fmincon"];
x_val = [x_1,x_2,x_3,x_4,x_5,x_6,x_7]; % array of x vals of functions
f_val = [f_1,f_2,f_3,f_4,f_5,f_6,f_7]; % array of f vals of functions

% array of error vals with respect to actual min val of x
err = zeros(1,length(x_val));
for i = 1:length(err)
    err(i) = (x_val(i) - x_act)./x_act;
end
% finding the index of in err array with value closest to 0 
err_min_index = find(err == min(abs(err)));

% computing the time to execute functions 
ft_1 = @()ALLFUNCS.IntervalHalving(obj_fun_min,xL,xU,e);
t1 = timeit(ft_1);
ft_2 = @()ALLFUNCS.GoldenSection(obj_fun_max,xL,xU,e);
t2 = timeit(ft_2);
ft_3 = @()ALLFUNCS.QuadraticEstimation(obj_fun_min,xL,xU);
t3 = timeit(ft_3);
ft_4 = @()ALLFUNCS.NewtonRaphson(obj_fun_der,x,x0,e);
t4 = timeit(ft_4);
ft_5 = @()ALLFUNCS.Bisection(obj_fun_der,x,xL,xU,e);
t5 = timeit(ft_5);
ft_6 = @()ALLFUNCS.Secant(obj_fun_der,x,xL,xU,e);
t6 = timeit(ft_6);
ft_7 = @()fmincon(obj_fun_min,x0,[],[],[],[],xL,xU);
t7 = timeit(ft_7);
ft_array = [t1,t2,t3,t4,t5,t6,t7];
% finding the index of item in ft_array with the minimum time value
ft_min_index = find(ft_array == min(ft_array));

% printing the result
fprintf('\nThe function with lowest error compared to actual minimum point')
fprintf('\nof objective function was found to be <strong>%s</strong> method\n', func_names(err_min_index))
fprintf('\nThe function with lowest execution time of optimization')
fprintf('\nof objective function was found to be <strong>%s</strong> method\n', func_names(ft_min_index))

