function [Sig, A, u, v] = genSig2(N0,N1,K,a_min,sigma_n)
% Generate 2-D signals with on-grid randomly distributed frequencies
% Copyright (c) 2018 Shaogang Wang
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation. This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY.

% Input: 
%   N0, N1: signal length of the two dimensions
%   K: sparsity level
%   a_min: the minimum sigal amplitude
%   sigma_n: noise variance 

% Output: 
%   Sig: signal matrix (N0-by-N1)
%   A: amplitude of each frequency (K-length vector)
%   u,v: frequency locations

    u = randi(N0,[K,1])-1;
    v = randi(N1,[K,1])-1;
    ind = [u,v];
    uInd = unique(ind,'rows');
    while size(uInd,1)<K
        u = randi(N0,[K,1])-1;
        v = randi(N1,[K,1])-1;
        uInd = [uInd; [u,v]];
        uInd = unique(uInd,'rows');
    end
    rdI = randperm(length(uInd));
    uInd = uInd(rdI(1:K),:);
    oriInd = sortrows(uInd);
    u = oriInd(:,1);
    v = oriInd(:,2);
    A = a_min*exp(1i*2*pi*rand(K,1));    
    Sig=zeros(N0,N1);
    Sig(sub2ind([N0,N1],u+1,v+1)) = A;
    Sig = N0*N1*ifft2(Sig);
    noise = sqrt(sigma_n)*sqrt(2)/2*(randn(N0,N1)+1i*randn(N0,N1));
    Sig = Sig+noise;

%% Visualization       
%     figure,mesh(0:N1-1,0:N0-1,abs(fft2(Sig)));
%     xlabel('$v$','Interpreter','LaTex'), ylabel('$u$','Interpreter','LaTex'), zlabel('Amplitude');
%     axis tight;
%     set(gca,'FontSize',20);
end