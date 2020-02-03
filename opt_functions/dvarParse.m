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
%   1       0 elements
%   2:7     1 element
%   8:22    2 elements
%   23:42   3 elements
%   43:49   4 elements
%   50:55   5 elements
%   56      6 elements%
%
%   dVar will have varying numbers of rows depending on the number of
%   elements, but will always have two columns: [node1 node2]

if dvarDim == 1
    % No nodal connections
    dVar = [0 0 0];
else
    if dvarDim == 56
        % All nodes connected
        dVar = [elPot(2) elPot(3);
            elPot(3) elPot(4);
            elPot(4) elPot(5);
            elPot(5) elPot(2);
            elPot(2) elPot(4);
            elPot(3) elPot(5)];
    else
        if dvarDim <= 7
            % One pairing of nodes
            group1 = [elPot(2) elPot(3);
                elPot(3) elPot(4);
                elPot(4) elPot(5);
                elPot(5) elPot(2);
                elPot(2) elPot(4);
                elPot(3) elPot(5)];
            dim = dvarDim - 1;
            dVar = group1(dim,:);
        else
            if dvarDim <= 22
                % Two pairings of nodes
                group2 = [elPot(3) elPot(4) elPot(4) elPot(5);
                    elPot(4) elPot(5) elPot(5) elPot(2);
                    elPot(5) elPot(2) elPot(2) elPot(3);
                    elPot(2) elPot(4) elPot(4) elPot(3);
                    elPot(2) elPot(4) elPot(4) elPot(5);
                    elPot(5) elPot(2) elPot(2) elPot(4);
                    elPot(3) elPot(2) elPot(2) elPot(4);
                    elPot(2) elPot(5) elPot(3) elPot(4);
                    elPot(2) elPot(3) elPot(3) elPot(4);
                    elPot(2) elPot(3) elPot(4) elPot(5);
                    elPot(4) elPot(3) elPot(3) elPot(5);
                    elPot(4) elPot(5) elPot(5) elPot(3);
                    elPot(2) elPot(5) elPot(5) elPot(3);
                    elPot(2) elPot(3) elPot(3) elPot(5);
                    elPot(2) elPot(4) elPot(3) elPot(5)];
                dim = dvarDim - 7;
                dVar = group2(dim,:);
                dVar = reshape(dVar,[2,2])';
            else
                if dvarDim <= 42
                    % Three pairings of nodes
                    group3 = [elPot(3) elPot(4) elPot(4) elPot(5) elPot(5) elPot(2);
                        elPot(4) elPot(5) elPot(5) elPot(2) elPot(2) elPot(3);
                        elPot(5) elPot(2) elPot(2) elPot(3) elPot(3) elPot(4);
                        elPot(2) elPot(3) elPot(3) elPot(4) elPot(4) elPot(5);
                        elPot(2) elPot(4) elPot(2) elPot(3) elPot(3) elPot(4);
                        elPot(2) elPot(4) elPot(3) elPot(4) elPot(5) elPot(4);
                        elPot(2) elPot(4) elPot(2) elPot(5) elPot(5) elPot(4);
                        elPot(2) elPot(4) elPot(2) elPot(5) elPot(2) elPot(3);
                        elPot(3) elPot(5) elPot(3) elPot(2) elPot(3) elPot(4);
                        elPot(3) elPot(5) elPot(3) elPot(4) elPot(4) elPot(5);
                        elPot(3) elPot(5) elPot(5) elPot(2) elPot(5) elPot(4);
                        elPot(3) elPot(5) elPot(5) elPot(2) elPot(3) elPot(2);
                        elPot(2) elPot(4) elPot(3) elPot(5) elPot(3) elPot(4);
                        elPot(2) elPot(4) elPot(3) elPot(5) elPot(4) elPot(5);
                        elPot(2) elPot(4) elPot(3) elPot(5) elPot(5) elPot(2);
                        elPot(2) elPot(4) elPot(3) elPot(5) elPot(2) elPot(3);
                        elPot(2) elPot(4) elPot(2) elPot(3) elPot(4) elPot(5);
                        elPot(2) elPot(4) elPot(4) elPot(3) elPot(2) elPot(5);
                        elPot(3) elPot(5) elPot(3) elPot(2) elPot(5) elPot(4);
                        elPot(3) elPot(5) elPot(5) elPot(2) elPot(3) elPot(4)];
                    dim = dvarDim - 22;
                    dVar = group3(dim,:);
                    dVar = reshape(dVar,[2,3])';
                else
                    if dvarDim <= 49
                        % Four pairings of nodes
                        group4 = [elPot(2) elPot(3) elPot(3) elPot(4) elPot(4) elPot(5) elPot(5) elPot(2);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(3) elPot(4) elPot(4) elPot(5);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(4) elPot(5) elPot(5) elPot(2);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(5) elPot(2) elPot(2) elPot(3);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(2) elPot(3) elPot(3) elPot(4);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(5) elPot(2) elPot(3) elPot(4);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(2) elPot(3) elPot(4) elPot(5)];
                        dim = dvarDim - 42;
                        dVar = group4(dim,:);
                        dVar = reshape(dVar,[2,4])';
                    else
                        % Five pairings of nodes
                        group5 = [elPot(2) elPot(3) elPot(3) elPot(4) elPot(4) elPot(5) elPot(5) elPot(2) elPot(2) elPot(4);
                            elPot(2) elPot(3) elPot(3) elPot(4) elPot(4) elPot(5) elPot(5) elPot(2) elPot(3) elPot(5);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(4) elPot(5) elPot(5) elPot(2) elPot(2) elPot(3);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(5) elPot(2) elPot(2) elPot(3) elPot(3) elPot(4);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(2) elPot(3) elPot(3) elPot(4) elPot(4) elPot(5);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(3) elPot(4) elPot(4) elPot(5) elPot(5) elPot(2)];
                        dim = dvarDim - 49;
                        dVar = group5(dim,:);
                        dVar = reshape(dVar,[2,5])';
                    end
                end
            end
        end
    end
    % Add a column of element numbers
    elNum = 1:1:size(dVar,1);
    dVar = [elNum.',dVar];
end
               
                
        

end

