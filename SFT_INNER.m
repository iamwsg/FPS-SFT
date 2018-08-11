function [Omega, A, P] = SFT_INNER(Sig, Win, Omega1, A1, N0, N1, epsilon, gamma, n_d, n_s)
% SFT_INNER
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
%   Omega1, A1: location and amplitude of previously recovered frequenceis
%   N0, N1: signal length of the two dimensions
%   epsilon: the threshold of detecting significant frequencies in a slice
%   gamma: the threshold for 1-sparsity detection
%   n_d, n_s: parameters of $n_d$-out-of-$n_s$ detection

% Output: 
%   Omega: frequency locations
%   A: complex amplitude of each frequency 
%   P: debug info

S = Sig; 
P = [];
Omega = [];
A = [];
oc = [];
it = 0;
L = round(lcm(N0,N1));
c0 = L/N0; c1 = L/N1;
while it < n_d
    it = it+1;    
    alpha0 = randi(N0,1)-1; 
    alpha1 = randi(N1,1)-1;     
    while gcd(alpha0,alpha1) ~= 1 || gcd(alpha0,c1) ~= 1 || gcd(alpha1,c0) ~= 1
        alpha0 = randi(N0,1)-1; 
        alpha1 = randi(N1,1)-1; 
    end
    
    tau0 = randi(N0,1)-1;
    tau1 = randi(N1,1)-1;
    
    [ind,ha, res] = triSlicing2(S, Win, Omega1, A1 ,alpha0, alpha1, tau0,tau1, epsilon, gamma,L);
        
    P{it,1} = [alpha0,alpha1];
    P{it,2} = [tau0, tau1];
    P{it,3} = ind;
    P{it,4} = ha;
    P{it,5} = res;
    
    if ~isempty(ha)
        if isempty(A)
            Omega = [Omega;ind];
            A = [A; ha.'];
            oc = [oc; ones(length(ha),1)];
        else
            [C,ia,ib] = intersect(Omega,ind,'rows');
            if ~isempty(C)
                A(ia) = A(ia)+ha(ib).';
                oc(ia) = oc(ia)+1;

                df = setdiff(1:size(ind,1), ib);
                Omega = [Omega; ind(df,:)];
                A = [A; ha(df).'];
                oc = [oc; ones(length(df),1)];
            else
                Omega = [Omega;ind];
                A = [A; ha.'];
                oc = [oc; ones(length(ha),1)];
            end            
            
        end
        
    end
    
end

%visual_votes(V,3);
i = find(oc>=n_s);
Omega = Omega(i,:);
A = A(i,:)./oc(i);

end

function [ind, ha, res] = triSlicing2(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1, epsilon, gamma,L)
    [N0,N1] = size(Sig);
    [in0, hs0, res0] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1, epsilon,L);
    [in1, hs1, res1] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0+1, tau1, epsilon,L);
    [in2, hs2, res2] = SLICING(Sig, Win, I1, hA1, alpha0, alpha1, tau0, tau1+1, epsilon,L);
    
    res = (res0+res1+res2)/3;
    in = intersect(in0, in1);
    in = intersect(in, in2);
    hs = zeros(3, length(in));
    
    if length(in) ~= length(in0) || length(in) ~= length(in1) || length(in) ~= length(in2)
        
        [~,ia0,~] = intersect(in0,in);
        [~,ia1,~] = intersect(in1,in);
        [~,ia2,~] = intersect(in2,in);
        
        hs(1,:) = hs0(ia0);
        hs(2,:) = hs1(ia1);
        hs(3,:) = hs2(ia2);
        
    else
        
        hs(1,:) = hs0;
        hs(2,:) = hs1;
        hs(3,:) = hs2;
    end
        
    u = [];
    v = [];
    ha = [];
    ahs = abs(hs);
    for i = 1:length(in)
        % 1-sparse detection
        if var(ahs(:,i))< gamma
            u0 = wrapTo2Pi(angle(hs(2,i)/hs(1,i)))*N0/(2*pi);
            u(end+1) = mod(round(u0), N0);
            v0 = wrapTo2Pi(angle(hs(3,i)/hs(1,i)))*N1/(2*pi);
            v(end+1) = mod(round(v0), N1); 
            %ha(end+1) = hs(1,i)*exp(-1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1));
            ha(end+1) = ((hs(1,i)+ hs(2,i)*exp(-1i*2*pi*u(end)/N0) + ...
                hs(3,i)*exp(-1i*2*pi*v(end)/N1)))*exp(-1i*2*pi*(u(end)*tau0/N0+v(end)*tau1/N1))/3;
        end
    end

    ind = [u; v]'; 
                
end

function [in, hs, res] = SLICING(Sig, Win, I, hA,alpha0, alpha1, tao0,tao1, epsilon,L)
    [N0, N1]=size(Sig);
    % sampling on a line
    l = 0:L-1;
    x = mod(alpha0*l+tao0,N0)+1;
    y = mod(alpha1*l+tao1,N1)+1;
    ind = sub2ind([N0,N1], x, y);
    sl = Sig(ind);
    
    % add window
    win = Win(ind);
    sl = sl.*win;
    
    % construct the line from found frequencies
    cl = zeros(1, L);  
    if ~isempty(hA)
        cl=CONSTRUCTION(I,hA,alpha0,alpha1,N0,N1,tao0,tao1,L);
    end
    
    
    dl = sl-cl;
    %dl = sl;
    res = norm(dl)/L;
    
    % detection of significant frequencies
    fsl = 1/L*fft(dl);
    ft= abs(fsl);
    
%     figure,stem(0:N0-1,ft,'LineWidth',2),grid;
%     xlabel('$\omega$','Interpreter','LaTex');
%     ylabel('Amplitude');
%     %legend('DFT grids','Ground truth','Estimated');
%     set(gca,'FontSize',20);
%     axis tight;
    
    in = find(ft>=epsilon); % the first detection
    hs = fsl(in);
    in  = in-1;
end


function f = CONSTRUCTION(I,f_hat,alpha0,alpha1,N0,N1,tao0,tao1,L)
    f_hat = f_hat.* exp(1i*2*pi*(I(:,1)/N0*tao0 + I(:,2)/N1*tao1));
    k= mod(round(I*[alpha0*L/N0 alpha1*L/N1]'),L)+1;
    fhat1=accumarray(k,f_hat,[L,1]);
    f=L*ifft(fhat1);
    f= f.';
end





