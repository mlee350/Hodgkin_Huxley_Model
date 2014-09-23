%%%Define Time

stepSize = 0.001;
time = linspace(1, (1/stepSize), 100); %Total Simulation time: 100ms

%%%Constants

g_K = 36; %mS/cm^2
g_Na = 120; %mS/cm^2
g_L = 0.3; %mS/cm^2

E_K = -12; %mV
E_Na = 115; %mV
E_L = 10.6; %mV

V_rest = -70; %mV

