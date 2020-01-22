function[g] = frameSteering(i,elConnec,nf,nodalCoords)
%% Element Steering Vector
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Forms the steering vector for element i, which ties the nodes in an
% element to their degree of freedom number
%
% Accepts 
%   the element number
%   the global connectivity matrix
%   the global nodal freedom matrix
%   the nodal coordinate matrix
%
% Returns 
%   the global stiffness matrix
%

% retrieve the nodes of element i
node1 = elConnec(i,1);
node2 = elConnec(i,2);

% form the steering vector from element's degrees of freedom
g=[nf(nodalCoords(:,1)== node1,1); ...
   nf(nodalCoords(:,1)== node1,2); ...
   nf(nodalCoords(:,1)== node1,3); ...
   nf(nodalCoords(:,1)== node2,1); ...
   nf(nodalCoords(:,1)== node2,2); ...
   nf(nodalCoords(:,1)== node2,3)];
end
