function[F] = formGlobalF(numFreeDOF,nnd,nodof,nf,load)
%% Form Global Force Vector
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Builds the global force vector from boundary conditions and prescribed
% loads
%
% Accepts 
%   number of unconstrained DOF, numFreeDOF
%   the total number of nodes, nnd
%   the number of nodal degrees of freedom, nodof
%   the nodal freedom matrix, nf
%   the loading array, load
%
% Returns 
%   the global forcing vector, F

% Initialize             
F = zeros(numFreeDOF,1);
% Build
for i = 1:nnd
    for j = 1:nodof
        if nf(i,j) ~= 0
            F(nf(i,j)) = load(i,j);
        end
    end
end
end