function [excitation_voltage,rated_output,max_voltage,min_voltage,max_load,min_load,zero_balance,sensitivity] = lcl_005(v_excitation,output,max_ds_load,zb)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%   lcl_005(10,2,0.113,0.3)
excitation_voltage = v_excitation;
%LCL-005
rated_output = output;   %mV/V
max_voltage = rated_output*excitation_voltage;
min_voltage = 0;    %assume 0 volts until zero balance
max_load = max_ds_load*9.807;    %0.113kgf*9.80665 N/kgf = N
min_load = 0;   %can we go to negative load?
zero_balance = zb; %mV/V
min_voltage = zero_balance*excitation_voltage;
sensitivity = (max_voltage-min_voltage)/(max_load-min_load);
params = [excitation_voltage,rated_output,max_voltage,min_voltage,max_load,min_load,zero_balance,sensitivity];
end

