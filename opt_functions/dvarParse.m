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
% Values between 1 and 140 correspond to different numbers of elements:
%   1:<20        0 elements
%   20:<40       1 element
%   40:<60       2 elements
%   60:<80       3 elements
%   80:<100      4 elements
%   100:<120     5 elements
%   120:140      6 elements
%
%   dVar will have varying numbers of rows depending on the number of
%   elements, but will always have two columns: [node1 node2]

if dvarDim < 20
    % No nodal connections
    dVar = [0 0 0];
else
    if dvarDim >= 120
        % All nodes connected
        dVar = [elPot(2) elPot(3);
            elPot(3) elPot(4);
            elPot(4) elPot(5);
            elPot(5) elPot(2);
            elPot(2) elPot(4);
            elPot(3) elPot(5)];
    else
        if dvarDim < 40
            % One pairing of nodes
            group1 = [elPot(2) elPot(3);
                elPot(3) elPot(4);
                elPot(4) elPot(5);
                elPot(5) elPot(2);
                elPot(2) elPot(4);
                elPot(3) elPot(5)];
            dim = dvarDim - 20;
            for ii = 1:6
               if dim < ii * (20/6)
                   dim = ii;
                   break
               end
            end
            dVar = group1(dim,:);
        else
            if dvarDim < 60
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
                dim = dvarDim - 40;
                for ii = 1:15
                    if dim < ii * (20/15)
                        dim = ii;
                        break
                    end
                end
                dVar = group2(dim,:);
                dVar = reshape(dVar,[2,2])';
            else
                if dvarDim < 80
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
                    dim = dvarDim - 60;
                    for ii = 1:20
                        if dim < ii * (20/20)
                            dim = ii;
                            break
                        end
                    end
                    dVar = group3(dim,:);
                    dVar = reshape(dVar,[2,3])';
                else
                    if dvarDim < 100
                        % Four pairings of nodes
                        group4 = [elPot(2) elPot(3) elPot(3) elPot(4) elPot(4) elPot(5) elPot(5) elPot(2);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(3) elPot(4) elPot(4) elPot(5);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(4) elPot(5) elPot(5) elPot(2);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(5) elPot(2) elPot(2) elPot(3);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(2) elPot(3) elPot(3) elPot(4);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(5) elPot(2) elPot(3) elPot(4);
                            elPot(2) elPot(4) elPot(3) elPot(5) elPot(2) elPot(3) elPot(4) elPot(5)];
                        dim = dvarDim - 80;
                        for ii = 1:7
                            if dim < ii * (20/7)
                                dim = ii;
                                break
                            end
                        end
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
                        dim = dvarDim - 100;
                        for ii = 1:6
                            if dim < ii * (20/6)
                                dim = ii;
                                break
                            end
                        end
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

