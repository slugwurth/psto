function [elPot] = elementPotential(nndx,nndy)
%% Element Potential Matrix Generation
%
%	Matt Ireland 2019
%
% Accepts
%   a grid size in nodal dimensions
% Returns
%   the element potential matrix
%       where:      (:,1) is the design variable number
%                   (:,2:5) are the clock-wise ordered nodes corresponding
%                   to each design variable
%                   

% Find element count
nel = (nndx-1) * (nndy-1);

% Initialize element potential matrix
elPot = zeros(nel,5);
elPot(:,1) = 1:nel;

% Populate element potential matrix
kk = 1;
for jj = 0:nndy-2
    for ii = 0:nndx-2
        % The four available nodes in an element potential space
        elPot(kk,2) = kk + jj;
        elPot(kk,3) = kk + jj + 1;
        elPot(kk,4) = kk + jj + nndx;
        elPot(kk,5) = kk + jj + nndx + 1;
        kk = kk + 1;
    end
end

end

