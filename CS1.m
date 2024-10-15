%% ESE441 Case Study 1: 
%% Introduction
% * Authors:                  Will Burgess, Mack LaRosa
% * Class:                    ESE 441
% * Date:                     Created 10/10/2024, Last Edited 10/17/2024
%% Housekeeping
close all
clear
clc
code = "finished";
%% Part 1: Model Sim using ODE45
figure;
for i =1:4
    if i == 1
    v = 0.1; % V1, infection rate (between 0 and 1)
    k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
    r = 0.9; % recovery rate
    a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)
    elseif i == 2
    v = 0.2;
    k = [10000, 2000];
    r = 0.01;
    a = 0.12;
    elseif i == 3
    v = 0.2;
    k = [10000, 2000];
    r = 0.9;
    a = 0.03;
    elseif i == 4
    v = 0.2;
    k = [20000, 100];
    r = 0.9;
    a = 0.1;
    end

u = [0,0]; % Control inputs
IC1 = [1e6 - 10,10]; %Initial susceptible, Initial infected 
t = 0:1:150;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(1); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2)];

[t, x1] = ode45(system, t, IC1);

%fh1 = figure(1)

subplot(2,2,i)
plot(t, x1(:,1), 'linewidth', 1.5);
hold on;
plot(t, x1(:,2), 'linewidth', 1.5);
%plot(t, immune, 'linewidth', 1.5);
title({sprintf('Zero-Input Simulation #%.1d',i), ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    sprintf('u = [%.1f, %.1f]', u(1), u(2))})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on
end

%% Part 1: Simulate System with Linearized Model
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)
t = 0:1:150;

IC1 = [1e6 - 10,10]; %Initial susceptible, Initial infected 

xeq =  [(r/k(2) + a)*(k(1)/v) + 0.15,0]; 

J = [0, (-1 * v * xeq(1))/k(1) + a; 
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

ivMax = -10000;
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

% Plot control input
figure
plot(t,u(:,1), 'linewidth', 2.5)
hold on
plot(t,u(:,2), 'linewidth', 1.5)
ylim([ivMax - 1000, 1000])
title('Control Inputs (Delayed Sine Wave)')
ylabel('Amplitude (# of Individuals)')
xlabel('Time (days)')
legend('u_{1}(t)','u_{2}(t)', Location='southeast');
grid on

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(round(t+1),1); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(round(t+1),2)];

[t, x1] = ode45(system, t, IC1);

% fh1 = figure(1);
figure
plot(t, x1(:,1), 'linewidth', 1.5);
hold on;
plot(t, x1(:,2), 'linewidth', 1.5);
%plot(t, immune, 'linewidth', 1.5);
title({'Control Input Simulation (Delayed Sine Wave)', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
   'u = [wave, 0]'})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on
%% Part 2: Design Original Control Inputs (Operation Cannibal Zombies)
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)
beta = 0.15;

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

%These control inputs can control at what population an endemic state is
%reached by changing beta
system_endemic = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a*x(2) + beta*x(2); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) - beta*x(2)];

%Modeling cannibal zombies will drive disease to eradication
system_zombies = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a*x(2); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) - beta*x(2)];

[t_zombies, x1_zombies] = ode45(system_zombies, t, IC1);
[t_endemic, x1_endemic] = ode45(system_endemic, t, IC1);

figure
subplot(2, 1, 1)
hold on;
plot(t_zombies, x1_zombies(:,1), 'linewidth', 1.5);
plot(t_zombies, x1_zombies(:,2), 'linewidth', 1.5);
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
title(sprintf('Cannibal Zombies: u = [%.1f * x_{2}(t), - %.1f * x_{2}(t)]', beta, beta))
hold off;
grid on

subplot(2, 1, 2)
hold on;
plot(t_endemic, x1_endemic(:,1), 'linewidth', 1.5);
plot(t_endemic, x1_endemic(:,2), 'linewidth', 1.5);
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
title(sprintf('Endemic Control: u = [%.1d, -%.1f * x_{2}(t)]', 0, beta))

sgtitle({'Control Input Simulation (Cannibal Zombies)', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    %sprintf('u = [%.1f, %.1f]', u(1), u(2)), ...
    },'FontSize', 12, 'FontWeight', 'bold')
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

% fh2 = figure(2);
% plot(x1(:,1), x1(:,2), 'linewidth', 1.5);
% title({'Zero-Input Simulation #1 (Trace)',sprintf('V_{1} =%.1f, K_{1} =%.1f, K_{2} =%.1f, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)});
% xlabel('X_{1}');
% ylabel('X_{2}');
% grid on