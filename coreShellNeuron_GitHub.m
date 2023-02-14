%% Induced current
filename = 'CouplingFactor_AlN_35-45um_extended.xlsx';
T = readtable(filename);
time_aux = T{:,4}*1e3;        %[ms]
amp_aux = T{:,8}'/1000;       %[\muA/cm^2; pA]; scaled amp

time = 0:1e-2:time_aux(end);   
i_induced = interp1(time_aux, amp_aux, time);
figure, plot(time,i_induced);
%% Parameter Set: ---
% NEURONAL CELLS
% Nominal values for electrical activity
% Voltages given in [mV]
% Conductances given in [nS or, equivalently, mS/cm^3]
% Time constants given in [ms]
% Capacitance given in [\muF/cm^2]
gK = 36;
gNa = 120;
gL = 0.3;
VNa = 50;
VK = -70;
VL = -54.4;
cm = 1;
%% Derivatives
% x(1) = Vm;
% x(2) = mNa;   -->m
% x(3) = mK;    -->n
% x(4) = hNa;   -->h
%-----------------2---------------------------------------------------------
f1 = @(t,x) [
    -(gK*x(3)^4*(x(1)-VK) + gNa*x(2)^3*x(4)*(x(1)-VNa) + gL*(x(1)-VL)-interp1(time,i_induced,t))/cm;
    0.1*(x(1)+40)/(1-exp(-(x(1)+40)/10))*(1-x(2))-4*exp(-(x(1)+65)/18)*x(2);
    0.01*(x(1)+55)/(1-exp(-(x(1)+55)/10))*(1-x(3))-0.125*exp(-(x(1)+65)/80)*x(3);
    0.07*exp(-(x(1)+65)/20)*(1-x(4))-1/(1+exp(-(x(1)+35)/10))*x(4);
    ];
% ODE with the initial condition
xinit = [-65 zeros(1,3)];
t_start = tic;
tspan_ode = [0 time(end)];
[t1, sol] = ode45(f1, time, xinit);
t_end = toc(t_start);
fprintf('\n     [ OK ]');
fprintf('\n     Elapsed %.4f s',t_end);
%% Induced Current & Plasma Membrane Potential
i_induced1 = interp1(time, i_induced, t1);
Vm = sol(:, 1);
figure,
subplot(2,2,[1,2]), plot(t1,i_induced1, 'LineWidth', 1); hold on; grid minor;
title('Induced current', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$i_{induced}$ [$\mu$A/cm$^2$]', 'Interpreter', 'latex', ...
    'FontSize', 14);
subplot(2,2,[3,4]), plot(t1, Vm, 'LineWidth', 1); hold on; grid minor;
title('Plasma membrane potential', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$v_m$ [mV]', 'Interpreter', 'latex', ...
    'FontSize', 14);