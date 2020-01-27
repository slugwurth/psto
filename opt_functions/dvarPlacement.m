function [elDist,randVar,g,nodalCoords] = dvarPlacement(nndx,nndy,scale,type)
%% Viable Random Element Distribution
%
% Matt Ireland 2019
%
% Accepts
%   a grid size in nodal dimensions
%   the unit cell scaling factor
%   the BC type
% Returns
%   -a global arrangement of node connections according to a random
%   selection of design variable values
%       that:       are fully connected (no islands)
%                   reach and connect BC nodes
%   -the corresponding random design variable values
%   -a graph/network representation of the elements
%   nodal numbers and positions for the connections

% Initialize the unit cell potential vector
uniPot = unitCellPotential();% [n 1 2 3 ... 57]

% Initialize the element potential matrix
elPot = elementPotential(nndx,nndy);% [node1 node2 node3 node4]

% Initialize the design variable value matrix
randVar = ones(size(elPot,1),2);% [n val]
randVar(:,1) = elPot(:,1);% first column corresponds to unit cell number

% Initialize the element distribution matrix
elDist = [];% unknown-by-three: [Nel node1 node2]

%% Find viable combinations
% Build a random distribution of elements according to a random selection
% of design value for each unit cell, then make sure the BC nodes
% are connected before proceeding


% Iterative search for viable combinations
viable = 0; iter = 1;
while viable == 0
    
    % for all of the unit cells
    for ii = 1:size(elPot,1)
        % Select randomly from unit cell potential
        randVar(ii,2) = datasample(uniPot,1,'Replace',false);
        % Parse design variable value into line elements for the unit cell
        unitElements = dvarParse(elPot(ii,:),randVar(ii,2));
        % For the nonzero case
        if unitElements(1,1) ~= 0
            % Change the local element numbering to the global values
            if ii > 1 && size(elDist,1) > 0
                unitElements(:,1) = unitElements(:,1) + elDist(end,1);
            end
            % Concatenate to the element distribution matrix
            elDist = [elDist; unitElements];
        end
    end
    
    % Loop again if selecting the null case
    if isempty(elDist)
        continue
    end
        
    % Represent connectivity as a network object
    g = graph(elDist(:,2),elDist(:,3));
    % Create node variables to check for connectivity
    if type == 1% (MBB loading)
        n1 = 1; n2 = nndx; n3 = nndx *(nndy-1) + 1;
    else
        if type == 2% (cantilever midplane tip loading)
            n1 = 1; n2 = nndx * ((nndy-1)/2 + 1); n3 = nndx *(nndy-1) + 1;
        end
    end
    % Put network components into bins
    [gbins,gbinsizes] = conncomp(g);
    % Build logical array of invalid components
    h = gbinsizes < max(gbinsizes) & gbinsizes > 1;
    
    % Make sure component is large enough
    if length(gbins) >= n3
        % Check if vital nodes are in the same component (connected)
        if gbins(n1) == gbins(n2) == gbins(n3)
            % End the search
            viable = 1;
        end
        % Search through nodes and check for disconnected components
        for ii = 1:length(gbins)
            % If node belongs to disconnected component
            if h(gbins(ii)) == 1
                % Mark element distribution as invalid
                viable = 0;
            end
        end
    end
    
    % Keep searching
    if viable == 0
        % Update the iteration count
        iter = iter + 1;
        % Clear element distribution
        elDist = [];
    end
end





%% Trim nodal grid
% Keep node numbering consistent with a full grid, but don't include unused
% nodes in the nodal coordinate matrix

% Build a nodal location grid
grid = nodalGrid(nndx,nndy,scale);

% Build vector of nodes connected by edges
nodeList = cat(1,elDist(:,2),elDist(:,3));

% Remove duplicate node entries
nodeList = unique(nodeList);

% Sort the nodes
nodeList = sort(nodeList);

% Initialize a location matrix
nodalCoords = zeros(size(nodeList,1),3);
nodalCoords(:,1) = nodeList;

% Search through the location grid and apply coordinates
for ii = 1:size(nodeList,1)
    nodalCoords(ii,2:3) = grid(nodeList(ii),2:3);
end

end

