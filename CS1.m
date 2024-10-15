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

%% Part 1: Model Sim using ODE45
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)
u = [0,0]; % Control inputs


IC1 = [1e6 - 10,10]; %Initial susceptible, Initial infected 
t = 0:1:150;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(1); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2)];

%options = odeset('Events', @events_function, 'NonNegative', [1, 2]);

[t, x1] = ode45(system, t, IC1);
%immune = ((IC1(1) + IC1(2)) * ones(length(x1),1)) - (x1(:,1) + x1(:,2));

%fh1 = figure(1)
figure;
plot(t, x1(:,1), 'linewidth', 1.5);
hold on;
plot(t, x1(:,2), 'linewidth', 1.5);
%plot(t, immune, 'linewidth', 1.5);
title({'Zero-Input Simulation #1 (Time Based)', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1f, K_{2} =%.1f, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    sprintf('u = [%.1f, %.1f]', u(1), u(2))})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on
% fh2 = figure(2);
% plot(x1(:,1), x1(:,2), 'linewidth', 1.5);
% title({'Zero-Input Simulation #1 (Trace)',sprintf('V_{1} =%.1f, K_{1} =%.1f, K_{2} =%.1f, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)});
% xlabel('X_{1}');
% ylabel('X_{2}');
% grid on

%% Part 1: Simulate System with Linearized Model
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)
t = 0:1:150;

IC1 = [1e6 - 10,10]; %Initial susceptible, Initial infected 

xeq =  [(r/k(2) + a)*(k(1)/v) + 0.1,0]; 

J = [0, (-1 * v * xeq(1))/k(1); 
     0 , (v * xeq(1))/k(1) - r/k(2) - a];

system = @(t, x) J * x;

[t, x] = ode45(system, t, IC1);

%fh2 = figure(1)
figure;
plot(t, x(:,1), 'linewidth', 1.5);
hold on;
plot(t, x(:,2), 'linewidth', 1.5);
%plot(t, immune, 'linewidth', 1.5);
title({'Linearized Simulation #1 (Time Based)', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    sprintf('u = [%.1f, %.1f]', u(1), u(2))})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on

%% Part 2: Implement Control Input
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)

IC1 = [1e6 - 10,10]; %Initial susceptible, Initial infected 
t = 0:1:150;

ivMax = 10000;
ivPeriod = 70; %days the intervention will last
ivStartTime = 20; % day you want to start control input

% Initialize control input (wave)
wave = zeros(length(t),1);

% Generate sine wave with given amplitude and period
sine = ivMax * sin(pi * (1/ivPeriod) * (0:(ivPeriod))); %pi instead of 2pi to get 1 half period

% Insert sine wave at specified start time
wave(ivStartTime:(ivStartTime + length(sine)-1),1) = sine;

% Control inputs
u = cat(2,wave(:), zeros(length(t),1));

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(round(t+1),1); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(round(t+1),2)];

%options = odeset('Events', @events_function, 'NonNegative', [1, 2]);

[t, x1] = ode45(system, t, IC1);
%immune = ((IC1(1) + IC1(2)) * ones(length(x1),1)) - (x1(:,1) + x1(:,2));

% fh1 = figure(1);
figure
plot(t, x1(:,1), 'linewidth', 1.5);
hold on;
plot(t, x1(:,2), 'linewidth', 1.5);
%plot(t, immune, 'linewidth', 1.5);
title({'Control Input Simulation #1 (Time Based)', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    sprintf('u = [%.1f, %.1f]', u(1), u(2))})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on

%% Part 3: Simulate Networked Model
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)

N = 4; % Number of regions
tspan = 0:1:150;

C = ones(N,length(tspan)); % Controls summation for x1 , rows regions, cols timestep

D = ones(N,length(tspan)); % Controls summation for x2 , rows regions, cols timestep

u = [0 0;
     0 0;
     0 0;
     0 0;]; % Control inputs, each row is a region

IC1 = [1e6 10;
       1e6 0;
       1e6 10;
       1e6 0]; % ICs for each region, First col susceptible, second infected

IC1 = IC1(:); % Convert to a vector for compatability with ODE45


[t, x] = ode45(@(t, x) networkedModel(t, x, N, v, k, r, a, u, C, D), tspan, IC1);

figure;
for i = 1:N
    subplot(N/2, 2, i);
    plot(t, x(:, 2*i-1), 'linewidth', 1.5);
    hold on
    plot(t, x(:, 2*i), 'linewidth', 1.5)
    title({sprintf('Region #%.1d', i), ...
           sprintf('u =[%.f,%.f]', u(i,1),u(i,2))});
    legend('Susceptible','Infected');
    grid on
end
sgtitle({'Networked Model Simulation #1', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)},...
    'FontSize', 12, 'FontWeight', 'bold')
% -------------------------------------------------------------------
%% Auxilliary Code : Meant Saving Information and Storing Old Code
% -------------------------------------------------------------------
%% Save images
% filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. Control Systems\Case Study 1\Figure export";
% exportgraphics(fh1, fullfile(filepath, 'V Distribution.jpg'), 'resolution', 300);
% exportgraphics(fh2, fullfile(filepath, 'E Magnitude.jpg'), 'resolution', 300);
% exportgraphics(fh3, fullfile(filepath, 'E Field.jpg'), 'resolution', 300);

%% Archive: Old Code That Might be Useful
% -------------------------------------------------------------------
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
%     ((v * x * y)/(k + y) - (r * y)/(y + l) - a*y];
%y = x2. k = k1, l = k2