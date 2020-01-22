function[KK] = placeKL(KK,kg,g,eldof)
%% Add to Global Stiffness Matrix
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Adds to the global stiffness matrix by element properties passed to it
%
% Accepts 
%   the global stiffness matrix
%   the element matrix in global coordinates
%   the element steering vector
%   the number of nodal degrees of freedom
%
% Returns 
%   the global stiffness matrix

for i = 1:eldof
    if g(i) ~= 0
        for j = 1:eldof
            if g(j) ~= 0
                KK(g(i),g(j)) = KK(g(i),g(j)) + kg(i,j);
            end
        end
    end
end
end

