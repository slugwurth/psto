function [force] = findElementForces(nel,nodalCoords,elConnec,prop,nodalFreeMat,eldof,delta)
%% Find Element Forces
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Calculates the forces in each element resulting from prescribed loading
%
% Accepts 
%   the number of elements
%   the global nodal coordinate matrix
%   the global element connectivity matrix
%   the global element property matrix
%   the global nodal freedom matrix
%   the number of elemental degrees of freedom
%   the nodal displacement matrix
%
% Returns 
%   the global elemental force vector

%% Initialize
% element displacement array
edg = zeros(eldof);
% force output vector
force = zeros(nel);  

%% Find Forces
% step through the elements
 for i = 1:nel
    % Form element matrix in local xy
    kl = frameLocalStiffness(i,nodalCoords,elConnec,prop);
    % Form transformation matrix
    C = frameTransform(i,nodalCoords,elConnec);
    % Transform the element matrix from local to global coordinates
    kg = C*kl*C'; 
    % Retrieve the element steering vector
    g = frameSteering(i,elConnec,nodalFreeMat,nodalCoords) ; 
    
    for j = 1:eldof
        if g(j) == 0
            % displacement = 0 for restrained freedom 
            edg(j) = 0;  
        else
            edg(j) = delta(g(j));
        end
    end
    
    % Element force vector in global XY 
    fg = kg*edg';
    % Element force vector in local  xy 
    fl = C'*fg ;
    % Write to output vector
    force(i) = fl(3);
 end
end

