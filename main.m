%%%%%%%%%%%%%%%%%%%%%
%%% Main Function %%%
%%%%%%%%%%%%%%%%%%%%%

%% Usage %%
% Simply run the program
% Most editable values to change the conditions can be found in 'setup.m'

%%%%% Problem Statement %%%%%%
% Solve 2-D compressible Navier-Stokes using MacCormack's method.
% 1.  Mach 4 flow over a flat plate with a sharp leading edge
% 2a. Generate a numerical schlieren
% 2b. Show that the simulation faithfully reproduces the theoretical Mach
%     angle trend
% 2c. Implement the adiabatic wall condition

%% Setup %%

clear; clc;     % Clear workspace
addpath('lib\') % Add helper function directory if not existing

% Setup function
setup

%% Solver %%

% Solve 2-D Compressible Navier-Stokes using MacCormack's method
maccormack

%% Graphics %%

% Animate solution
animate

%% Notes %%
% Reshaping is slightly faster than squeeze if size of array is known
% because squeeze recalculates the size of the input array