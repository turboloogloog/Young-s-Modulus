clc;clear all;close all;
%% set parameters
alfa = 3;
c = 3.3;
Lamda_ref = 785e-9;
Lamda = 785e-9;
A0 = Lamda_ref/2*2;     % Vibration amplitude (m)
ft = 40;                % Vibration frequency ( Hz)
fs = 20000;             % Sampling frequency (Hz)
n = fs/ft*8;            % Sampling length
t = (0:n-1)/fs;         % time index

%% simulating

[emit_laser, result_y, result_g] = SMI(alfa,c,Lamda,ft,fs,n);
    
%% figrue sketch

sketch_all(emit_laser, result_y, result_g,n,fs,t)

sketch_emit_laser(emit_laser,n,t)   
spectrum_emit_laser(emit_laser,fs,n)

sketch_SMI(result_g,t)
spectrum_SMI(result_g,n,fs)