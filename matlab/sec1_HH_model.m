% Section 1: Hudgkin Huxley model
clear;
warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Set up directory (check if it exists)
foldername = 'sec1_HH_model' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

%===simulation time===
simulationTime = 100; %in milliseconds
deltaT=.01;
t=0:deltaT:simulationTime;


%===specify the external current I===
changeTimes = [0]; %in milliseconds
currentLevels = [50]; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)

%Set externally applied current across time
%Here, first 500 timesteps are at current of 50, next 1500 timesteps at
%current of zero (resets resting potential of neuron), and the rest of
%timesteps are at constant current
% I(1:500) = currentLevels; I(501:2000) = 0; I(2001:numel(t)) = currentLevels;
I(1:500) = 0; I(501:1000) = currentLevels; I(2001:4000) = 0; I(4001:numel(t)) = currentLevels;
%Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
%I(1:numel(t)) = currentLevels;


%===constant parameters===%
%All of these can be found in Table 3
gbar_K=36; gbar_Na=120; g_L=.3;
E_K = -12; E_Na=115; E_L=10.6;
C=1;


%===set the initial states===%
V=0; %Baseline voltage
alpha_n = .01 * ( (10-V) / (exp((10-V)/10)-1) ); %Equation 12
beta_n = .125*exp(-V/80); %Equation 13
alpha_m = .1*( (25-V) / (exp((25-V)/10)-1) ); %Equation 20
beta_m = 4*exp(-V/18); %Equation 21
alpha_h = .07*exp(-V/20); %Equation 23
beta_h = 1/(exp((30-V)/10)+1); %Equation 24

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18


for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
   
    %---calculate the coefficients---%
    %Equations here are same as above, just calculating at each time step
    alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
   
   
    %---calculate the currents---%
    I_Na = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
    I_L = g_L *(V(i)-E_L); %Equation 5
    I_ion = I(i) - I_K - I_Na - I_L;
   
   
    %---calculate the derivatives using Euler first order approximation---%
    V(i+1) = V(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i)); %Equation 16

end


V = V-70; %Set resting potential to -70mv

%===plot Voltage===%
plot(t,V,'LineWidth',3)

% Save
filename = 'trial3.mat';
dir_file = fullfile(dir_folder, filename);
save(dir_file, 't', 'V', 'n', 'm', 'h');
