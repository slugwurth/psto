function [dVar] = dvarParse(elPot,dvarDim)
%% Design Variable Parser
%
%	Matt Ireland 2019
%
% Accepts
%   an array of nodes for the unit cell
%   a value of the design variable for the unit cell
% Returns
%   nodal connections within the unit cell corresponding to the value of
%   the design variable
%       where:      (:,1) is the element number
%                   (:,2:3) are the node pairs corresponding to each line
%                   element in the unit cell
%                   

%% Assign Connectivity
% Values between 1 and 56 correspond to different numbers of elements:
%   1:28    0 elements
%   29:56   6 elements%
%
%   dVar will have varying numbers of rows depending on the number of
%   elements, but will always have two columns: [node1 node2]

if dvarDim <= 28
    % No nodal connections
    dVar = [0 0];
else
    if dvarDim > 28
        % All nodes connected
        dVar = [elPot(2) elPot(3);
            elPot(3) elPot(4);
            elPot(4) elPot(5);
            elPot(5) elPot(2);
            elPot(2) elPot(4);
            elPot(3) elPot(5)];
        % Add a column of element numbers
        elNum = 1:1:size(dVar,1);
        dVar = [elNum.',dVar];
    end
end


 

end

