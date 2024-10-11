%% ESE441 Case Study 1: 
%% Introduction
% * Authors:                  Will Burgess, Mack LaRosa
% * Class:                    ESE 441
% * Date:                     Created 10/10/2024, Last Edited 10/18/2024
%%
code = "finished";
v = 2; % V1, infection rate
k = [1,2].';% Sat constant for: infection, recovery
r = 0.1; % recovery rate

x = [1e6,1]; % Initial: susceptible, infected individuals
a = 1e-3; % Rate of reinfection/loss of immunity
u = [0,0]; % Control inputs

xdot = [-1*((v * x(1) * x(2))/(k(1)+x(2)))-a*x(2) + u(1);
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2)];