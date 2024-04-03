clear;
spec_final;

% baseline
clear; close all; bias_adjust=0; break_debt=0;
main_baseline;
clear; close all; bias_adjust=1; break_debt=0;
main_baseline;

% Simulation
main_simulated_bias;
