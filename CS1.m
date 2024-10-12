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
% v = 0.05; % V1, infection rate (between 0 and 1)
% k = [100000,1].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
% r = 1000; % recovery rate
% 
% x = [1e6,1]; % Initial: susceptible, infected individuals
% a = -0.01; % Rate of reinfection/loss of immunity (hundreds place)
% u = [0,3e3]; % Control inputs
% 
% xdot = [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(1);
%     ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2)];

% xpplane = [-1*((v * x * y)/(k+y))+ a * y;
%     ((v * x * y)/(k + y)) - (r * y)/(y + l) - a*y];
%y = x2. k = k1, l = k2
%%
v = 0.5; % V1, infection rate (between 0 and 1)
k = [1000000,1000000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 20; % recovery rate

%x = [1e6,1]; % Initial: susceptible, infected individuals
a = 0.8; % Rate of reinfection/loss of immunity (hundreds place)
u = [0,3e3]; % Control inputs
IC1 = [1e6,10]; %Initial susceptible, Initial susceptible 
t = 0:10:1000;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))-a*x(2) + u(1);
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a *x(2) + u(2)];

options = odeset('Events', @events_function, 'NonNegative', [1, 2]);

[t, x1] = ode45(system, t, IC1, options);

fh1 = figure(1);
plot(t, x1(:,1), 'linewidth', 1.5);
hold on;
plot(t, x1(:,2), 'linewidth', 1.5);
title({'Zero-Input Simulation #1 (Time Based)',sprintf('V_{1} =%.1f, K_{1} =%.1f, K_{2} =%.1f, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)});
xlabel('Time (weeks)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on

% fh2 = figure(2);
% plot(x1(:,1), x1(:,2), 'linewidth', 1.5);
% title({'Zero-Input Simulation #1 (Trace)',sprintf('V_{1} =%.1f, K_{1} =%.1f, K_{2} =%.1f, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)});
% xlabel('X_{1}');
% ylabel('X_{2}');
% grid on


%% Save images
% filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. Control Systems\Case Study 1\Figure export";
% exportgraphics(fh1, fullfile(filepath, 'V Distribution.jpg'), 'resolution', 300);
% exportgraphics(fh2, fullfile(filepath, 'E Magnitude.jpg'), 'resolution', 300);
% exportgraphics(fh3, fullfile(filepath, 'E Field.jpg'), 'resolution', 300);

%% Functions 

% Event function to stop when x1 (susceptible population) reaches 0
function [value, isterminal, direction] = events_function(t, x)
    x1 = x(1); % Susceptible individuals

    % Stop the solver when x1 = 0
    value = x1;         % Detect when x1 = 0
    isterminal = 0;     % Stop the integration
    direction = -1;     % Only trigger when x1 is decreasing
end