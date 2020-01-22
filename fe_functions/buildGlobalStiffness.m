function [KK] = buildGlobalStiffness(numFreeDOF,nel,nodalCoords,elDist,prop,nodalFreeMat,eldof)
%% Build the Global Stiffness Matrix
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Adds to the global stiffness matrix by element properties passed to it
%
% Accepts 
%   the number of free degrees of freedom
%   the number of elements
%   the global nodal coordinate matrix
%   the global element distribution
%   the global element property matrix
%   the global nodal freedom matrix
%   the number of elemental degrees of freedom
%
% Returns 
%   the global stiffness matrix

%% Initialize 
KK = zeros(numFreeDOF); 
elConnec = elDist(:,2:3);

%% Populate global stiffness matrix
for i = 1:nel
    % Form element matrix in local xy 
    kl = frameLocalStiffness(i,nodalCoords,elConnec,prop);       
    
    % Form transformation matrix
    C = frameTransform(i,nodalCoords,elConnec);       
    
    % Transform the element coordinates from local to global
    kg = C*kl*C';
    
    % Retreive the element steering vector                 
    g = frameSteering(i,elConnec,nodalFreeMat,nodalCoords);
    
    % Assemble the global stiffness matrix
    KK = placeKL(KK,kg,g,eldof);
end
end

