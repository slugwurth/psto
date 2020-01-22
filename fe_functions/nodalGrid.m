function [grid] = nodalGrid(nndx,nndy,scale)
%% Generate Nodal Coordinates
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Builds a unit-spaced nodal grid according to prescribed dimensions
%
% Accepts 
%   a grid size in nodal dimensions
%   a scaling factor
%
% Returns 
%   the properly formatted nodal coordinate matrix

% Initialize coordinate array
grid = zeros(nndx*nndy,3);

% Increment along x and y modulo grid size
kk = 1;
for ii = 1:nndy
    for jj = 1:nndx
        grid(kk,1) = kk;
        grid(kk,2) = mod(jj - 1,nndx);
        grid(kk,3) = mod(ii - 1,nndy);
        kk = kk + 1;
    end
end

% Scale the grid dimnensions
grid(:,2:3) = grid(:,2:3).*scale;

end
