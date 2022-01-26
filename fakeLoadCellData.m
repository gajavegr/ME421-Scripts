%create sample data for load cell interpretation
%For an unbalance of 2275kg-mm^2 we expect:
% Fc = 6.804045460N (~6.8N)

%record Length: 10s
clear
close all
clc
Tr = 10;    %s
fs = 1000;  %Hz
rpm=300;    %rpm
omega=rpm*(1/60)*2*pi;  %rad/s
expected_force = 6.8;   %N
%% define sensor parameters
[excitation_voltage,rated_output,max_voltage,min_voltage,max_load,min_load,zero_balance,sensitivity] = lcl_005(10,2,0.113,0.3);
%% fake/noisy data time
t = 0:1/fs:Tr;
omega_hz = 1/(2*pi)*omega;  %Hz
b = 2*pi*omega_hz;
measurand = expected_force*sin(b*t)+0.5*rand(size(t));   %N
v = measurand*sensitivity;
figure
plot(t,measurand)
xlabel("time (s)")
ylabel("Force (N)")
figure
plot(t,v*0.001)
xlabel("time (s)")
ylabel("Voltage (V)")
%% create data file
% datafile =   [
%     ["Channels","1"];
%     ["Samples","37800"];
%     ["Date","1/17/2022"];
%     ["Time","15:26:56"];
%     ["Y_Unit_Label","Voltage (V)"];
%     ["X_Dimension","Time (s)"];
%     ["X0","26:00.0"];
%     ["Delta_X","0.001"];
%     ["***End_of_Header***",""];
%     ["X_Value","01/17/2022 03:25:59 PM - Voltage - cDAQ4Mod1_ai0"];
%     ];
datafile = zeros(10,2);
datafile = vertcat(datafile,horzcat(t',v'));
csvwrite('sample_data.csv',datafile);
%% get output

% output_newtons = v*(1/sensitivity);
% sample_rate = 1/(t(2)-t(1));    %Hz
