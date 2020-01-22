function [randEl,g,nodeLocs] = elPlacement(nndx,nndy,elPot)
%% Viable Random Element Distribution
%
%	Matt Ireland 2019
%
% Accepts
%   a grid size in nodal dimensions
%   the element potential matrix
% Returns
%   a random arrangement of line elements
%       that:       are non-overlapping
%                   are fully connected
%                   reach and connect MBB BC nodes
%   a graph/network representation of the elements
%   nodal numbers and positions for the line elements

% Initialize random element matrix
randEl = zeros(size(elPot,1),3);
randEl(:,1) = elPot(:,1);

%% Find viable combinations

% Iterative search for viable combinations
viable = 0; iter = 1;
while viable == 0
    
    % Select randomly from element potential
    for ii = 1:size(elPot,1)
        randEl(ii,2:3) = datasample(elPot(ii,2:end),2,'Replace',false);
    end
    
    % Represent connectivity as a network object
    g = graph(randEl(:,2),randEl(:,3));
    
    % Search through for duplicates and reassign
    overlap = 1;
    while overlap == 1
        
        % Initial check and reassignment
        for jj = 1:size(elPot,1)
            % Temp variable for test edge
            edg = randEl(jj,2:3);
            % Temp Variable for new edge potential
            npot = elPot(jj,2:end);
            % Temp variable for new choices
            nedg = setdiff(npot,edg);
            
            % Search network for duplicate edge
            if edgecount(g,edg(1),edg(2)) > 1
                randEl(jj,2) = nedg(1);
                randEl(jj,3) = nedg(2);
            end
        end
        
        % Refresh network definition
        g = graph(randEl(:,2),randEl(:,3));
        
        % Check for duplicate edges after reassignment
        for kk = 1:size(elPot,1)
            if edgecount(g,randEl(kk,2),randEl(kk,3)) > 1
                break
            else
                if edgecount(g,randEl(kk,2),randEl(kk,3)) == 1 && kk == size(elPot,1)
                    overlap = 0;
                end
            end
        end
        
    end
    
    % Re-define new network object
    g = graph(randEl(:,2),randEl(:,3));
    
    % Create node variables to check for connectivity
    n1 = 1; n2 = nndx; n3 = nndx *(nndy-1) + 1;
    
    % Put network components into bins
    [gbins,gbinsizes] = conncomp(g);
    
    % Make sure component is large enough
    if length(gbins) >= n3
        % Check if vital nodes are in the same component (connected)
        if gbins(n1) == gbins(n2) == gbins(n3)
            % End the search
            disp('Viable Candidate Found');
            disp(['Iteration: ' num2str(iter)]);
            viable = 1;
        end
    end
    
    if viable == 0
        % Keep searching
        disp('Searching for viable candidate');
        disp(['Iteration: ' num2str(iter)]);
        % Update the iteration count
        iter = iter + 1;
    end
end

%% Remove invalid components

% Build logical array of invalid components
h = gbinsizes < max(gbinsizes) & gbinsizes > 1;

% Search through nodes and remove edges from invalid components
for ii = 1:length(gbins)
    % If node belongs to invalid component
    if h(gbins(ii)) == 1
        
        % Rename variable for clarity
        node = ii;
        % Convert edge table to array for search
        edgeArray = table2array(g.Edges);
        % Find the edge index for the corresponding node
        [idx,~] = find(edgeArray == node);
        % Remove the edge
        g = rmedge(g,idx);
    end
end

% Update the edge array to the correct length
edgeArray = table2array(g.Edges);

% Replace the element matrix with a single valid component
randEl = zeros(size(edgeArray,1),3);
randEl(:,1) = 1:1:size(edgeArray,1);
randEl(:,2:3) = edgeArray;

% Re-define the network object
g = graph(randEl(:,2),randEl(:,3));

%% Trim nodal grid

% Build a nodal location grid
grid = nodalGrid(nndx,nndy);

% Build vector of nodes connected by edges
nodeList = cat(1,randEl(:,2),randEl(:,3));

% Remove duplicate node entries
nodeList = unique(nodeList);

% Sort the nodes
nodeList = sort(nodeList);

% Initialize a location matrix
nodeLocs = zeros(size(nodeList,1),3);
nodeLocs(:,1) = nodeList;

% Search through the location grid and apply coordinates
for ii = 1:size(nodeList,1)
    nodeLocs(ii,2:3) = grid(nodeList(ii),2:3);
end

end
