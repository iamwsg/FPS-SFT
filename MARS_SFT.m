function [Omega,A,P] = MARS_SFT(Sig,Win,N0,N1,T ,epsilon, gamma, n_d, n_s)
% Robutst MARS_SFT algorithm; for FPS_SFT, set Win as rect window, epsilon =
% gamma = 1e-10; n_d = n_s = 1
%
% Copyright (c) 2018 Shaogang Wang
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation. This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY.

% Input:
%   Sig: singal matrix 
%   Win: window matrix
%   N0, N1: signal length of the two dimensions
%   epsilon: the threshold of detecting significant frequencies in a slice
%   gamma: the threshold for 1-sparsity detection
%   n_d, n_s: parameters of $n_d$-out-of-$n_s$ detection

% Output: 
%   Omega: frequency locations
%   A: complex amplitude of each frequency 
%   P: debug info

Omega = [];
A = [];
i = 0;
while i<T
    i = i+1;
    [ind,ha,P1] = SFT_INNER(Sig,Win,Omega,A,N0,N1 ,epsilon, gamma, n_d, n_s);
    
    if ~isempty(ha)
        if isempty(A)
            Omega = [Omega;ind];
            A = [A; ha];
        else
            [C,ia,ib] = intersect(Omega,ind,'rows');
            if ~isempty(C)
                A(ia) = A(ia)+ha(ib);

                df = setdiff(1:size(ind,1), ib);
                Omega = [Omega; ind(df,:)];
                A = [A; ha(df)];
            else
                Omega = [Omega;ind];
                A = [A; ha];
            end            
            
        end
        
    end

    P{i,1} = ind;
    P{i,2} = ha;
    P{i,3} = mean(cell2mat(P1(:,5)));    
end
end