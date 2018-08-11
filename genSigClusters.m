function [Sig, A, u, v] = genSigClusters(N0, N1, K, c_size, a_min, sigma_n)
% Generate 2-D signals with on-grid clustered frequencies
% Copyright (c) 2018 Shaogang Wang
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation. This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY.

% Input: 
%   N0, N1: signal length of the two dimensions
%   K: [CAUSION] this is the number of clusters
%   c_size: cluster size; the extent of each direction from the center. If
%   c_size = 1; each cluster contains 9(3-by-3) frequencies
%   a_min: the minimum sigal amplitude
%   sigma_n: noise variance 

% Output: 
%   Sig: signal matrix (N0-by-N1)
%   A: amplitude of each frequency
%   u,v: frequency locations

    u1 = randi(N0,[K,1])-1;
    v1 = randi(N1,[K,1])-1; 
    u = mod(u1-c_size:u1+c_size, N0);
    v = mod(v1-c_size:v1+c_size, N1);    
    ind = [];
    for i = 1:length(u1)
        for j = 1:2*c_size+1
            for k = 1:2*c_size+1
                in = [mod(u1(i)-j+c_size,N0), mod(v1(i)-k+c_size,N1)];
                ind = [ind; in];
            end
        end
    end
    
    uInd = unique(ind,'rows');
    K = size(uInd,1);
    rdI = randperm(size(uInd,1));
    uInd = uInd(rdI(1:K),:);
    oriInd = sortrows(uInd);
    u = oriInd(:,1);
    v = oriInd(:,2);
    fu = u/N0;
    fv = v/N1;
    Sig=zeros(N0,N1);
    A = a_min*exp(1i*2*pi*rand(K,1));
    for ii=1:K
        sig=A(ii)*exp(1i*2*pi*fu(ii)*(0:N0-1));
        h=exp(1i*2*pi*fv(ii)*(0:N1-1));
        Sig = Sig+sig.'*h;
    end
    noise = sigma_n*sqrt(2)/2*(randn(N0,N1)+1i*randn(N0,N1));
    Sig = Sig+noise;
end

    
