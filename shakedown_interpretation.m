clear
close all
clc

%read the data file
data_with_header = xlsread('shakedown_prep');
%sample data
data_with_header = csvread('sample_data.csv');
data = data_with_header(11:end,:);
t = data(:,1);
v = data(:,2);
num_channels = data_with_header(1,2)
num_samples = data_with_header(2,2)
acquisition_date = data_with_header(3,2)
acquisition_time = data_with_header(4,2)

%% measured value definition
m_measured=0;
r_measured=0;
d_measured=0;
s_measured=0;
rpm_measured = 0;
omega_measured = rpm_measured *(1/60)*2*pi;  %rad/s
U_measured = m_measured*r_measured*d_measured;
F_experimental_expected = U_measured*(omega_measured^2)/s_measured;
%% Expected value definition
U_expected = 0;
m_expected = 0;
r_expected=0;
d_expected=0;
s_expected=0;
rpm_expected= 0;
omega_expected = rpm_expected *(1/60)*2*pi;  %rad/s
F_theoretical_expected = U_expected*(omega_expected^2)/s_expected;
%% sensor parameters
% excitation_voltage = 10;
%LCL-005
% rated_output = 2e-3;   %mV/V
% max_voltage = rated_output*excitation_voltage;
% min_voltage = 0;    %assume 0 volts until zero balance
% max_load = 0.113*9.80665;    %0.113kgf*8.80665 N/kgf = N
% min_load = 0;   %can we go to negative load?
% zero_balance = 0.3e-3; %mV/V
% min_voltage = zero_balance*excitation_voltage;
% sensitivity = (max_voltage-min_voltage)/(max_load-min_load);
[excitation_voltage,rated_output,max_voltage,min_voltage,max_load,min_load,zero_balance,sensitivity] = lcl_005(10,2,0.113,0.3);
%% get output
output_newtons = v*(1/sensitivity);
sample_rate = 1/(t(2)-t(1));    %Hz
%% visualize data
close all
data_visualize = figure;
plot(t,output_newtons,'b')
xlabel('time (s)')
ylabel('force (N)')
hold on
yline(0);
yline(max_load,'--r');
yline(min_load,'--k');
%% get a dft of the force data
close all
y = output_newtons;
%sample rate (Hz) is known
fs = sample_rate;

N = length(y);
dt = 1/fs; % time increment
% t = [0:dt:(N-1)*dt]; % vector of times

% create a figure with two subplots
dft_plots = figure;

subplot(3,1,1) % the first subplot is the data y(t)
plot(t,y)
xlabel('Time(s)'), ylabel('Force N'), grid on

% compute the array of complex dft coefficients
ck_complex = (1/N)*fft(y);

% obtain the magnitudes of the complex numbers
ck = abs(ck_complex);

% create the frequency array fk corresponding to the ck array
Tr = N*dt; % record length
df = 1/Tr; % frequency increment
fk = [0:df:(N-1)*df]; % frequency array of length N

subplot(3,1,2); % second subplot is full spectrum
stem(fk,ck, 'filled', 'MarkerSize', 4)
xlabel('Frequency (Hz)'), ylabel('Amplitude (c_k) [N]'), grid on

% f is the array of physically meaningful frequencies
f = fk(1:N/2); % frequency array of length N/2

% A is the array of physically meaningful amplitudes of length N/2
A(1) = ck(1); % first element of the amplitude array
A(2:N/2) = 2*ck(2:N/2); % remaining elements of array are doubled
subplot(3,1,3); % graph the corrected spectrum
stem(f,A, 'filled', 'MarkerSize', 4)
xlabel('Frequency (Hz)'), ylabel('Force [N]'), grid on
axis([0 310 0 max(A)*1.01])
%% Get force at 300 RPM
target_rpm = 300;
omega_300=target_rpm*(1/60)*2*pi;  %rad/s
omega_300_hz = 1/(2*pi)*omega_300;  %Hz
[closest_freq,target_rpm_index]=min(abs(f-omega_300_hz));
% target_rpm_index = find(f==omega_300_hz);
f_measured = A(target_rpm_index);
%% Get Force vs/ Unbalance calibration curve
measured_unbalances = [1500,1800,1950,2250,2275,2333,2799,3033,3499];   %kg-mm^2
measured_forces = [4.486183820, 5.383420583, 5.832038965, 6.729275729, 6.804045460, 6.977511234, 8.371219007, 9.071063683, 10.46477146];    %N
coefficients = polyfit(measured_unbalances, measured_forces, 1);
unbalanceFit = linspace(min(measured_unbalances), max(measured_unbalances), 1000);
forceFit = polyval(coefficients , unbalanceFit);

%visualize calibration curve and data
figure
plot(measured_unbalances, measured_forces, 'b.', 'MarkerSize', 15);
hold on;
plot(unbalanceFit, forceFit, 'r-', 'LineWidth', 2);
grid on;
title("Unbalance Calibration curve")
xlabel("Unbalance (kg-mm^2)")
ylabel("Force (N)")

%% Calculate uncertainty in measured parameters

%uncertainty in Force
w_f = 1.067; %N
wf_rel = (100*w_f)/f_measured

%uncertainty in Force due to DFT
w_dft = sqrt((1/(2*N))*w_f^2)

%uncertainty in mass
w_m = sqrt(2)*0.005;  %g
w_m_rel = 100*w_m/m_measured

%uncertainty in radius
w_r = 0.014;  %mm
w_r_rel = 100*w_r/r_measured

%uncertainty in couple distance
w_d = 0.014;  %mm
w_d_rel = 100*w_d/d_measured

%uncertainty in bearing distance
w_s = 0.014;  %mm
w_s_rel = 100*w_s/s_measured

%uncertainty in unbalance
dudr = m_measured*d_measured;
dudm = r_measured*d_measured;
dudd = m_measured*r_measured;
w_u = (dudr*w_r)^2+(dudm*w_m)^2+(dudd*w_d)^2

%uncertainty in slope
mean_force = mean(measured_forces);
inner_sum = 0;
for i=1:length(measured_forces)
   inner_sum = inner_sum + (measured_forces(i)-mean_force)^2;
end
w_slope = (1/inner_sum)*w_f
%% Percent difference between expected and experimental values
%Percent diff in expected/measured force
f_error = abs(f_measured - F_experimental_expected)/((f_measured+F_experimental_expected)/2)*100
%Percent diff in expected/measured unbalance
u_error = abs(U_measured - U_expected)/((U_measured+U_expected)/2)*100