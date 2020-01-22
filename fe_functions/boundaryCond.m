function [nodalFreeMat,numFreeDOF,numNodes,imposedLoad] = boundaryCond(nndx,nndy,scale,nodof,nodalCoords,type,pload)
%% Boundary Conditions
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Accepts 
%   a grid size in nodal dimensions
%   nodal degrees of freedom (count)
%   a concise nodal position matrix
%   the type of boundary conditions
% Returns 
%   a nodal freedom matrix
%   the number of free dof
%   the number of nodes
%   a loading matrix

% DOF structure: [x y r] for nodof = 3

%% BC
% 1 is unconstrained, 0 is constrained

% Build a nodal location grid
grid = nodalGrid(nndx,nndy,scale);

% Find the number of nodes
numNodes = length(nodalCoords);
   
% Initialize nodal freedom matrix
nodalFreeMat = ones(numNodes,nodof);

% Apply load case-specific BC
if type == 1% MBB load case
    % Midplane rollers will always be at position x = 0
    
    % Index by node number into freedom array
    nodalFreeMat(nodalCoords(:,2) == 0,1) = 0;            % Midplane     rollers in Y
    nodalFreeMat(nodalCoords(:,1) == grid(nndx,1),2) = 0; % Lower-right  roller in X
end

% Count the free degrees of freedom
numFreeDOF = 0;
for i = 1:numNodes
    for j = 1:nodof
        if nodalFreeMat(i,j) ~= 0
            numFreeDOF = numFreeDOF+1;
            nodalFreeMat(i,j) = numFreeDOF;
        end
    end
end

%% Loading
% Initialize load matrix
imposedLoad = zeros(numNodes,3); % [x y r]

if type == 1                        % MBB load case
    % Find node number
    L1 = grid(nndx*(nndy-1)+1,1);   % Top-left
    
    % Find index of load node
    idxL1 = nodalCoords(:,1) == L1;
    
    % Index by node number into load array
    imposedLoad(idxL1,:) = [0 -pload 0];  % Top-Left, down in Y
end

end
