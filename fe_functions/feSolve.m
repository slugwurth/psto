function [delta,node_disp] = feSolve(numNodes,nodof,nodalFreeMat,KK,F)
%% Finite Element Solver
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Accepts 
%   the number of nodes, numNodes
%   the number of nodal DOF, nodof
%   the global nodal freedom matrix, nodalFreeMat
%   the global stiffness matrix, KK
%   the global forcing vector, F
%   
% Returns 
%   the unknown nodal displacements, delta
%   the global nodal displacements, node_disp

% Initialize
node_disp = zeros(numNodes,nodof);
            
% Solve for unknown displacements
delta = KK\F ;
% Extract nodal displacements
for i = 1:numNodes
    for j = 1:nodof
        node_disp(i,j) = 0;
        if nodalFreeMat(i,j) ~= 0
            node_disp(i,j) = delta(nodalFreeMat(i,j));
        end
    end
end

end

