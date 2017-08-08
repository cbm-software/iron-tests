clc
clear all
close all

%% FILE DEFINITIONS
% 
LOADING = 'SHEAR'; % SHEAR; UNIAX
REFINEMENT = 'MEDIUM'; % COARSE; MEDIUM; FINE
DIMENSION = '2D'; % 2D; 3D
INTERPOLATION = 1; %1; 2
CONTROL = 'DISPLACEMENT'; % DISPLACEMENT; FORCE

addpath('../../scripts/compare_abaqus_opencmiss/')
compare_solutions_main