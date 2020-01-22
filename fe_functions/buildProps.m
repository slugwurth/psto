function [prop] = buildProps(nel,E,A,I)
%% Build Element Property Matrix
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Accepts 
%   the number of elements
%   the nominal element modulus
%   the element sectional area
%   the element moment of inertia
%
% Returns 
%   the global element property matrix

% Initialize prop matrix
prop = ones(nel,3);
% Apply mechanical props to columns
prop(:,1) = prop(:,1) * E;
prop(:,2) = prop(:,2) * A;
prop(:,3) = prop(:,3) * I;

end

