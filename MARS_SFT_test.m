%% Test of Robust MARS-SFT
%
% Copyright (c) 2018 Shaogang Wang
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation. This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY.

clear;
clc;
close all;

N0=256;
N1=256;
N = N0*N1;
K = 10;  % signal sparsity
L = lcm(N0,N1); % line length
epsilon = 1e-10; % Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
gamma = 1e-10; % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
n_s = 1;
n_d = 1;
T = 10;

%% Generate signal
SNR_dB = 20; %signal SNR in dB
SNR = 20^(SNR_dB/20);
sigma_n = 0; % set this as 0 for noiseless signals, otherwise set as 1
a_min = 20^(SNR_dB/20); % the minimum signal power 

% Generate on-grid frequencies
%[Sig, A, u, v] = genSig2(N0,N1,K,a_min,sigma_n); 

% Generate off-grid frequencies
%[Sig, A, u, v] = genSigOffGrid(N0,N1,K,a_min,sigma_n); 

% Generate on-grid frequency clusters
[Sig, A, u, v] = genSigClusters(N0,N1,K,2,a_min,sigma_n); K = K*25;

Omega_gt = [u,v];  % ground truth frequency locations

%% Generate window; genearate rect window for on-grid cases
%% Chebyshev window
% att = 70; % PSR in dB
% rho_w = 10^(att/10);
% win0 = chebwin(N0,att);
% win1 = chebwin(N1,att);

%% Rectangular window
win0 = ones(N0,1);
win1 = ones(N1,1);
Win = win0*win1.';

[Omega, hA, P] = MARS_SFT(Sig, Win, N0, N1, T, epsilon, gamma, n_s, n_d);

% delta = 0.001;
% [I,hA,P] = DL_SFT_Phase_Reduction(Sig, N0,N1 ,epsilon, gamma, delta);

[Omega, ind] = sortrows(Omega);
hA = hA(ind);

Omega_gt = round(Omega_gt);
visual_localization(N0,N1,Omega_gt,Omega);

