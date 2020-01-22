function[C] = frameTransform(i,nodalCoords,elConnec)
%% Element Transformation Matrix
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% This function forms the transformation between local and global coordinates
%
% Accepts 
%   the element number
%   the nodal coordinate matrix
%   the element connectivity matrix
%
% Returns 
%   the element transformation matrix

% retrieve the nodes of element i
node1 = elConnec(i,1);
node2 = elConnec(i,2);

% Retrieve the x and y coordinates of nodes 1 and 2
x1 = nodalCoords(nodalCoords(:,1) == node1,2); y1 = nodalCoords(nodalCoords(:,1) == node1,3);
x2 = nodalCoords(nodalCoords(:,1) == node2,2); y2 = nodalCoords(nodalCoords(:,1) == node2,3);

% Evaluate the angle that the member makes with the global axis X
if (x2 - x1) == 0
    if (y2 > y1)
        theta = 2*atan(1);
    else
        theta = -2*atan(1);
    end
else
    theta = atan((y2-y1)/(x2-x1));
end

% Construct the transformation matrix
C = [cos(theta)   -sin(theta)    0         0            0          0 ; ...
     sin(theta)    cos(theta)    0         0            0          0 ; ...
         0             0         1         0            0          0 ; ...
         0             0         0     cos(theta)   -sin(theta)    0 ; ...
         0             0         0     sin(theta)    cos(theta)    0 ; ...
         0             0         0          0           0          1 ];
end
