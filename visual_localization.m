function [] = visual_localization(N0,N1,Omega_gt,Omega)    
% visual_localization
%
% Copyright (c) 2018 Shaogang Wang
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation. This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY.

% Input:
%   N0, N1: signal length of the two dimensions
%   Omega_gt: groud truth locations of frequencies
%   Omega: estimated locations of frequencies

% Output: 

    V = zeros(N0,N1);
    [x0,y0] = ind2sub(size(V),find(V>=0));
    figure,scatter(x0-1,y0-1,5,'fill');
    hold on;scatter(Omega_gt(:,1),Omega_gt(:,2),100,'Or');
    if ~isempty(Omega)
        hold on;scatter(Omega(:,1),Omega(:,2),50,'black','fill');
    end
    hold off;
    xlabel('$m_0$','Interpreter','LaTex');
    ylabel('$m_1$','Interpreter','LaTex');
    legend('DFT grid','Ground truth','Estimated');
    set(gca,'FontSize',20);
    axis tight;
    
    %% Add a small detailed graph
%     axes('Position',[.18 .18 .3 .3]);
%     x = 70:90;
%     y = 175:195;
%     [X,Y] = meshgrid(x,y);
%     
%     box on
%     scatter(reshape(X,[1,numel(X)]),reshape(Y,[1,numel(Y)]),1,'fill'),grid;
%     hold on;scatter(Ori(:,1),Ori(:,2),100,'Or');
%     hold on;scatter(I(:,1),I(:,2),50,'black','fill');
%     axis([x(1), x(end), y(1), y(end)]);
end