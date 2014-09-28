%%%Define Time

steps = 0.001 ;
time = 0:steps:100; %Total Simulation time: 100ms

%%%Constants

g_K_max = 36; %mS/cm^2
g_Na_max = 120; %mS/cm^2
g_L_max = 0.3; %mS/cm^2

E_K = -12; %mV
E_Na = 115; %mV
E_L = 10.6; %mV

V_rest = 0; %mV
C_m = 1; %uf/cm^2

%%%Initial Condition

%Calculate initial conditions based on V_rest
a_m = 0.1*((25-V_rest)/(exp((25-V_rest)/10) - 1));
B_m = 4*exp(-V_rest/18);
a_n = .01 * ((10-V_rest)/(exp((10-V_rest)/10)-1));
B_n = .125*exp(-V_rest/80);
a_h = .07*exp(-V_rest/20);
B_h = 1/(exp((30-V_rest)/10) + 1);

m_0 = a_m/(a_m+B_m);
n_0 = a_n/(a_n+B_n);
h_0 = a_h/(a_h+B_h);

%Create membrane potential vector and initialize first element to 0
V_m = zeros(size(time));
V_m(1) = 0;

%Calculate required values for the first iteration
m = m_0;
n = n_0;
h = h_0;

I_Na = m^3*g_Na_max*h*(V_m(1)-E_Na); 
I_K = n^4*g_K_max*(V_m(1)-E_K);
I_L = g_L_max*(V_m(1)-E_L);
I_ion = -I_Na - I_K - I_L;
    
g_Na_vec = zeros(size(time));
g_K_vec = zeros(size(time));
g_L_vec = zeros(size(time));

g_Na = m^3*h*g_Na_max;
g_K = n^4*g_K_max;
g_L = g_L_max;

g_Na_vec(1) = g_Na;
g_K_vec(1) = g_K;
g_L_vec(1) = g_L;

%%%Iterate + Eulers

%Start at second iteration since first iteration is already calculated
%based on initial conditions
for i = 2:length(time)

    %Eulers on membrane potential
    V_m(i) = V_m(i-1) + (steps*(I_ion / C_m));
    
    %Calculate Current
    I_Na = m^3*g_Na_max*h*(V_m(i)-E_Na); 
    I_K = n^4*g_K_max*(V_m(i)-E_K);
    I_L = g_L_max*(V_m(i)-E_L);
    I_ion = 5-I_Na - I_K - I_L;
  
    a_m = 0.1*((25-V_m(i))/(exp((25-V_m(i))/10) - 1));
    B_m = 4*exp(-V_m(i)/18);
    a_n = .01 * ((10-V_m(i))/(exp((10-V_m(i))/10)-1));
    B_n = .125*exp(-V_m(i)/80);
    a_h = .07*exp(-V_m(i)/20);
    B_h = 1/(exp((30-V_m(i))/10) + 1);
    
    %Eulers on m,n,h variables
    m = m + (steps*(a_m*(1-m)-B_m*m));
    n = n + (steps*(a_n*(1-n)-B_n*n));
    h = h + (steps*(a_h*(1-h)-B_h*h));
    
    %Update conductances
    g_Na = m^3*h*g_Na_max;
    g_K = n^4*g_K_max;
    g_L = g_L_max;
    
    g_Na_vec(i) = g_Na;
    g_K_vec(i) = g_K;
    g_L_vec(i) = g_L;
end

%%%Plot

%Compensate for -70 mV since V_rest initialized as 0
V_m = V_m - 70;

figure
plot(time,V_m)
ylim([-100, 70]);

figure
plot(time, g_Na_vec, time, g_K_vec, time, g_L_vec) 
axis([0, 100, 0, 50]);

