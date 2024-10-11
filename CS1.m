%% ESE441 Case Study 1: 
%% Introduction
% * Authors:                  Will Burgess, Mack LaRosa
% * Class:                    ESE 441
% * Date:                     Created 10/10/2024, Last Edited 10/18/2024
%% Housekeeping
close all
clear
clc
code = "finished";
%%
v = 0.05; % V1, infection rate (between 0 and 1)
k = [1,2].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 2; % recovery rate

x = [1e6,1]; % Initial: susceptible, infected individuals
a = 0.01; % Rate of reinfection/loss of immunity (hundreds place)
u = [0,0]; % Control inputs

xdot = [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(1);
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2)];

% xpplane = [-1*((v * x * y)/(k+y))+ a * y;
%     ((v * x * y)/(k + y)) - (r * y)/(y + l) - a*y];
%y = x2. k = k1, l = k2
%%

IC1 = [1e6,10]; %Initial susceptible, Initial susceptible 
t = 0:1:60;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))-a*x(2) + u(1);
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a *x(2) + u(2)];

[t, x1] = ode45(system, t, IC1);

fh1 = figure(1);
plot(t, x1(:,1).', 'linewidth', 1.5);
hold on;
plot(t, x1(:,2), 'linewidth', 1.5);
title({'Zero-Input Simulation #1',sprintf('V_{1} =%.1f, K_{1} =%.1f, K_{2} =%.1f, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)});
xlabel('Time (weeks)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on

%% Save images
% filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. Control Systems\Case Study 1\Figure export";
% exportgraphics(fh1, fullfile(filepath, 'V Distribution.jpg'), 'resolution', 300);
% exportgraphics(fh2, fullfile(filepath, 'E Magnitude.jpg'), 'resolution', 300);
% exportgraphics(fh3, fullfile(filepath, 'E Field.jpg'), 'resolution', 300);