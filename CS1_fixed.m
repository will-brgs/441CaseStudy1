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
%% General Notes
% Various pieces of code are repreated for different sections. This is done
% to allow flexibility for sections to be run independently. Overall code
% length will be shortened by removing these, but it will harm the ability
% to perform the different analysis techniques and strategiges used in this
% case study.
%% Part 1: Model Sim using ODE45
for i =1:4
    if i == 1
    v = 0.1; % V1, infection rate
    k = [100000,20000].';% Sat constant for: infection, recovery
    r = 0.9; % recovery rate
    a = 0.02; % Rate of reinfection/loss of immunity
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
IC = [1e6 - 10,10]; %Initial susceptible, Initial infected 
t = 0:1:150;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(1); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2)];

[t, x] = ode45(system, t, IC);

% Plot general system response
fh1 = figure(1);
subplot(2,2,i)
plot(t, x(:,1), 'linewidth', 1.5);
hold on;
plot(t, x(:,2), 'linewidth', 1.5);
title({sprintf('Zero-Input Simulation #%.1d',i), ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    sprintf('u = [%.1d, %.1d]', u(1), u(2))})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on
end

%Zoom in on oscilations at endemic state
x = x(75:151,:);
t = t(75:151);
fh2 = figure(2);
subplot(2,1,1)
plot(t, x(:,1), 'linewidth', 1.5);
xlabel('Time (days)');
ylabel('# of Individuals');
title('Susceptible Individuals')
grid on

colors = get(gca, 'ColorOrder');
orange = colors(2, :);

subplot(2,1,2)
plot(t, x(:,2), 'linewidth', 1.5, 'Color', orange);
xlabel('Time (days)');
ylabel('# of Individuals');
title('Infected Individuals')
sgtitle({'Zero-input Simulation #4 at Endemic State'},'FontSize', 12, 'FontWeight', 'bold')
grid on
%% Part 1: Simulate System with Linearized Model at Xeq
for i =1:2
   v = 0.1; % V1, infection rate
   k = [100,20000].';% Sat constant for: infection, recovery
   r = 0.9; % recovery rate
   a = 0.02; % Rate of reinfection/loss of immunity
t = 0:1:150000;
u = [0,0]; % control input
IC = [1e6 - 10,10]; %Initial susceptible, Initial infected
if i==1
   xeq =  [(r/k(2) + a)*(k(1)/v)+1,0]; % +1 will case the system to be unstable
else
   xeq =  [(r/k(2) + a)*(k(1)/v)-1,0]; % -1 will case the system to be stable
end

J = [0, ((-1 * v * xeq(1))/k(1)) + a;
    0 , ((v * xeq(1))/k(1)) - r/k(2) - a]; %Derrived jacobian matrix

system = @(t, x) J * x;

[t, x] = ode45(system, t, IC);

fh3 = figure(3);
subplot(2,1,i)
plot(t, x(:,1), 'linewidth', 1.5);
hold on;
plot(t, x(:,2), 'linewidth', 1.5);
title({'Linearized Simulation #1', ...
   sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
   sprintf('u = [%.1d, %.1d]', u(1), u(2))})
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
grid on
end
%% Part 2: Implement Control Input
v = 0.1; % V1, infection rate
k = [100000,20000].';% Sat constant for: infection, recovery
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity

IC = [1e6 - 10,10]; %Initial susceptible, Initial infected 
t = 0:1:150;

ivMax = -10000; %maximum peoeple treated per day in sine wave
ivPeriod = 70; %days the intervention will last
ivStartTime = 20; % day you want to start control input

wave = zeros(length(t),1); % zero padding
sine = ivMax * sin(pi * (1/ivPeriod) * (0:(ivPeriod))); %pi instead of 2pi to get 1 half period
wave(ivStartTime:(ivStartTime + length(sine)-1),1) = sine; % Insert sine wave at specified start time

u = cat(2,wave(:), zeros(length(t),1)); % Finalized control inputs

% Plot control input
fh5 = figure(5);
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

[t, x1] = ode45(system, t, IC);

fh6 = figure(6);
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

IC = [1e6 - 10,10]; %Initial susceptible, Initial infected 
t = 0:1:150;

%Modeling cannibal zombies will drive disease to eradication
sysZombies = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a*x(2); 
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) - beta*x(2)];

[t, xZombies] = ode45(sysZombies, t, IC);

fh7 = figure(7);
plot(t, xZombies(:,1), 'linewidth', 1.5);
hold on
plot(t, xZombies(:,2), 'linewidth', 1.5);
xlabel('Time (days)');
ylabel('# of Individuals');
legend('Susceptible','Infected');
title(sprintf('Cannibal Zombies: u = [%.1f * x_{2}(t), - %.1f * x_{2}(t)]', beta, beta))
grid on
title({'Control Input Simulation (Cannibal Zombies)', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a), ...
    sprintf('Cannibal Zombies: u = [%.1f * x_{2}(t), - %.1f * x_{2}(t)]', beta, beta)
    },'FontSize', 12, 'FontWeight', 'bold')
grid on

%% Part 3: Networked Model Simulation: No Control Input
%parameters
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)

N = 3; % number of regions

C = [0 0.2 0.1;
     0.2 0 0.1;
     0.1 0.1 0;]; %Coupling matrix for susceptible

D = [0    0.1  0.05;
     0.1  0    0.05;
     0.05 0.05 0;]; %Coupling matrix for Infected. Half as likely to travel

u = [0, 0, ...
     0, 0, ...
     0, 0]; % Control inputs, formatting helps vizualize regions(row) and sus/inf(col)

IC = [3e5 - 10 , 10,...
      2e5      , 0, ...
      5e5      , 0];
% Initial conditions: formatting helps vizualize regions(row) and sus/inf(col)
t = 0:1:150;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(1) - C(1,2)*(x(1)-x(3)) - C(1,3)*(x(1)-x(5)); %region 1 susceptible
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(2) - D(1,2)*(x(2)-x(4)) - D(1,3)*(x(2)-x(6)); %region 1 infected
    -1*((v * x(3) * x(4))/(k(1)+x(4)))+ a * x(4) + u(1) - C(2,1)*(x(3)-x(1)) - C(2,3)*(x(3)-x(5)); %region 2 susceptible
    ((v * x(3) * x(4))/(k(1) + x(4))) - (r * x(4))/(x(4) + k(2)) - a*x(4) + u(2) - D(2,1)*(x(4)-x(2)) - D(2,3)*(x(4)-x(6)); %region 2 infected
    -1*((v * x(5) * x(6))/(k(1)+x(6)))+ a * x(6) + u(1) - C(3,1)*(x(5)-x(1)) - C(3,2)*(x(5)-x(3)); %region 3 susceptible
    ((v * x(5) * x(6))/(k(1) + x(6))) - (r * x(6))/(x(6) + k(2)) - a*x(6) + u(2) - D(3,1)*(x(6)-x(2)) - D(3,2)*(x(6)-x(4))];

[t, x] = ode45(system, t, IC);

fh8 = figure(8);
for i = 1:N
    subplot(N, 1, i);
    plot(t, x(:, 2*i-1), 'linewidth', 1.5);
    hold on
    plot(t, x(:, 2*i), 'linewidth', 1.5)
    title({sprintf('Region #%.1d', i), ...
           'No Control Inputs'});
    legend('Susceptible','Infected');
    grid on
end
sgtitle({'Networked Model Simulation #1', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)},...
    'FontSize', 12, 'FontWeight', 'bold')
%% Part 3: Networked Model Simulation: Delayed Sine Input for Region 3
% Redefine delayed sine wave
ivMax = -10000;
ivPeriod = 70; %days the intervention will last
ivStartTime = 5; % day you want to start control input

wave = zeros(length(t),1); % zero padding
sine = ivMax * sin(pi * (1/ivPeriod) * (0:(ivPeriod))); %pi instead of 2pi to get 1 half period
wave(ivStartTime:(ivStartTime + length(sine)-1),1) = sine; % Insert sine wave at specified start time

% Control inputs
wave = cat(2,wave(:), zeros(length(t),1));
%parameters
v = 0.1; % V1, infection rate (between 0 and 1)
k = [100000,20000].';% Sat constant for: infection, recovery IMPORTANT CONSTRAINT
r = 0.9; % recovery rate
a = 0.02; % Rate of reinfection/loss of immunity (hundreds place)

N = 3; % number of regions

C = [0 0.2 0.1;
     0.2 0 0.1;
     0.1 0.1 0;]; %Coupling matrix for susceptible

D = [0    0.1  0.05;
     0.1  0    0.05;
     0.05 0.05 0;]; %Coupling matrix for Infected. Half as likely to travel

noInput = zeros(length(t),1); %Zero pad for time variant input
u = [noInput,noInput, ...
     wave   ,noInput, ...
     noInput,noInput]; % Control inputs, formatting helps vizualize regions(row) and sus/inf(col)

IC = [3e5 - 10 , 10,...
      2e5      , 0, ...
      5e5      , 0];
% Initial conditions: formatting helps vizualize regions(row) and sus/inf(col)
t = 0:1:150;

system = @(t, x) [-1*((v * x(1) * x(2))/(k(1)+x(2)))+ a * x(2) + u(round(t+1),1) - C(1,2)*(x(1)-x(3)) - C(1,3)*(x(1)-x(5)); %region 1 susceptible
    ((v * x(1) * x(2))/(k(1) + x(2))) - (r * x(2))/(x(2) + k(2)) - a*x(2) + u(round(t+1),2) - D(1,2)*(x(2)-x(4)) - D(1,3)*(x(2)-x(6)); %region 1 infected
    -1*((v * x(3) * x(4))/(k(1)+x(4)))+ a * x(4) + u(round(t+1),3) - C(2,1)*(x(3)-x(1)) - C(2,3)*(x(3)-x(5)); %region 2 susceptible
    ((v * x(3) * x(4))/(k(1) + x(4))) - (r * x(4))/(x(4) + k(2)) - a*x(4) + u(round(t+1),4) - D(2,1)*(x(4)-x(2)) - D(2,3)*(x(4)-x(6)); %region 2 infected
    -1*((v * x(5) * x(6))/(k(1)+x(6)))+ a * x(6) + u(round(t+1),5) - C(3,1)*(x(5)-x(1)) - C(3,2)*(x(5)-x(3)); %region 3 susceptible
    ((v * x(5) * x(6))/(k(1) + x(6))) - (r * x(6))/(x(6) + k(2)) - a*x(6) + u(round(t+1),6) - D(3,1)*(x(6)-x(2)) - D(3,2)*(x(6)-x(4))]; %region 3 infected

[t, x] = ode45(system, t, IC);

fh9 = figure(9);
for i = 1:N
    subplot(N, 1, i);
    plot(t, x(:, 2*i-1), 'linewidth', 1.5);
    hold on
    plot(t, x(:, 2*i), 'linewidth', 1.5)
    if i == 2
    title({sprintf('Region #%.1d', i), ...
           'Vaccine Program Intervention'});
    else
    title({sprintf('Region #%.1d', i), ...
           'No Control Inputs'});
    end
    legend('Susceptible','Infected');
    grid on
end
sgtitle({'Networked Model Simulation #2', ...
    sprintf('V_{1} =%.1f, K_{1} =%.1d, K_{2} =%.1d, r =%.1f, \\alpha =%.4f', v, k(1), k(2), r, a)},...
    'FontSize', 12, 'FontWeight', 'bold')
%% Save images
% filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. Control Systems\Case Study 1\Figure export";
% exportgraphics(fh1, fullfile(filepath, 'part1 different vars.jpg'), 'resolution', 300);
% exportgraphics(fh2, fullfile(filepath, 'part1 equilibrium zoom in.jpg'), 'resolution', 300);
% exportgraphics(fh3, fullfile(filepath, 'part1 linearized sim 1.jpg'), 'resolution', 300);
% fh4 no longer exists
% exportgraphics(fh5, fullfile(filepath, 'part2 delayed sine input u.jpg'), 'resolution', 300);
% exportgraphics(fh6, fullfile(filepath, 'part2 delayed sine sim.jpg'), 'resolution', 300);
% exportgraphics(fh7, fullfile(filepath, 'part2 zombies sim.jpg'), 'resolution', 300);
% exportgraphics(fh8, fullfile(filepath, 'part3 networked sim no control.jpg'), 'resolution', 300);
% exportgraphics(fh9, fullfile(filepath, 'part3 networked sim control.jpg'), 'resolution', 300);
