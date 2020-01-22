function[kl] = frameLocalStiffness(i,nodalCoords,elConnec,prop)
%% Form Local Element Stiffness Matrix
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Forms the beam-column element stiffness in local coordinates
%
% Accepts 
%   the element number
%   the global nodal coordinate matrix
%   the global element connectivity matrix
%   the global element property matrix
%
% Returns 
%   the local element stiffness matrix

% retrieve the nodes of element i
node1 = elConnec(i,1);
node2 = elConnec(i,2);

% Retrieve the x and y coordinates of nodes 1 and 2
x1 = nodalCoords(nodalCoords(:,1) == node1,2); y1 = nodalCoords(nodalCoords(:,1) == node1,3);
x2 = nodalCoords(nodalCoords(:,1) == node2,2); y2 = nodalCoords(nodalCoords(:,1) == node2,3);

% Evaluate length of element i
L = sqrt((x2-x1)^2 + (y2-y1)^2);

% Retrieve section properties of element i
E = prop(i,1); 
A = prop(i,2); 
I = prop(i,3);

% Rename for clarity below
EA = E*A; EI = E*I;

% Calculate element stiffness matrix in its local coordinates   
kl =   [EA/L       0           0           -EA/L        0            0      ; ...
         0      12*EI/L^3     6*EI/L^2       0      -12*EI/L^3     6*EI/L^2 ; ...
         0      6*EI/L^2      4*EI/L         0       -6*EI/L^2     2*EI/L   ; ...
       -EA/L       0             0          EA/L        0            0      ; ...
         0     -12*EI/L^3    -6*EI/L^2       0       12*EI/L^3     -6*EI/L^2 ; ...
         0      6*EI/L^2      2*EI/L         0       -6*EI/L^2      4*EI/L   ];

end

