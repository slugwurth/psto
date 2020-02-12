function [prop] = pInterp(prop,randEl,nodalCoords,scanDir,maxRed)
%% Process Interpreter
%
%	Matt Ireland 2019
%
% Accepts 
%   a matrix Nx3 of element material properties
%   a nodal connectivity matrix (N-1)x3
%   a nodal coordinate matrix Nx3
%   the scanning direction vector 1x2
%   the scalar maximum percent (decimal) reduction in stiffness
% Returns 
%   the modulated element material property matrix Nx3

% Note:
%   The row indices in prop and randEl correspond to elements
%   and values in randEl(:,2:3) correspond to row indices in nodalCoords

%% Validate
% Check the scan direction vector
% if scanDir(1) == 1 && scanDir(2) == 0
%     disp('Valid direction vector');
% else
%     if scanDir(1) == 0 && scanDir(2) == 1
%         disp('Valid direction vector');
%     else
%         error('Invalid direction vector');
%     end
% end

%% Initialize
% Find the midpoints of the elements

% Build an element coordinate data structure
elCoords = zeros(size(randEl,1),5);% Nx5 [n x1 y1 x2 y2]
% Populate element coordinates
elCoords(:,1) = randEl(:,1);% n
for ii = 1:2% node 1 or 2
    for kk = 1:size(randEl,1)% by element
        row = nodalCoords(:,1) == randEl(kk,ii+1);
        elCoords(kk,(2*ii):(2*ii)+1) = nodalCoords(row,2:3);
    end
end

% Find the midpoint coordinates for each element
elMids = zeros(size(randEl,1),3);% Nx3 [n x y]
elMids(:,1) = elCoords(:,1);% n
elMids(:,2) = (elCoords(:,2) + elCoords(:,4))/2;% x midpoint
elMids(:,3) = (elCoords(:,3) + elCoords(:,5))/2;% y midpoint

%% Sort Elements
% Use the scan direction vector [x y] to decide a sort order

% Sort by scan direction then by scan order
if scanDir(2) == 1 % Scanning in the Y direction in X order
    [elMids,index] = sortrows(elMids,[3 2]);
else% Scanning in the X direction in Y order
    [elMids,index] = sortrows(elMids,[2 3]);
end

% Write out sorted elements
% fprintf('[Element Xmid Ymid] \nsorted by scan direction:\n');
% disp(num2str(elMids));

% Sort the element coordinate matrix to match
elCoords = elCoords(index,:);

%% Assign Orientation Weights
% Use the scan direction to assign weights to parallel, perpindicular, or
% oblique elements

% Normalized oblique weight value
wMid = 0.5;

% Initialize an element weight matrix
weight = zeros(size(elMids,1),2);% Nx2 [n w]
% Populate element numbers
weight(:,1) = elMids(:,1);% n

% Assign weights according to parallelism
for ii = 1:size(elCoords,1)
    if elCoords(ii,2) == elCoords(ii,4)% vertical element
        % Check parallelism
        if scanDir(1) == 1% horizontally travelling scan
            weight(ii,2) = 1;% element parallel to scanline
        end
    else
        if elCoords(ii,3) == elCoords(ii,5)% horizontal element
            % Check parallelism
            if scanDir(2) == 1% vertically travelling scan
                weight(ii,2) = 1;% element parallel to scanline
            end
        else % off-axis element
            weight(ii,2) = wMid;% element oblique to scanline
        end
    end
end

%% Find Average Inter-Element Lag
% Use midpoint locations

% Initialize a lag matrix
lag = zeros(size(elMids,1),2);% Nx2 [n d]
% Populate element numbers
lag(:,1) = elMids(:,1);% n

% Horizontal, Scanning in the X direction in Y order
if scanDir(1) == 1
    for ii = 1:size(elMids,1)
        x = elMids(ii,2);% current X coordinate
        y = elMids(elMids(:,2) == x,3);% y coordinates of in-line elements
        % Find number of in-line elements
        nel = size(y,1);
        % Check if solo element
        if nel == 1
            % Solo elements have no inter-element distance
            continue
        end
        % Find inter-element distances
        z = circshift(y,1);
        id = y(2:end) - z(2:end);
        % Maximum inter-element distance
        dMax = max(id);
        % Check for overlapping elements
        if dMax == 0
            % Overlapping elements have no inter-element distance
            continue
        end
        % Find gap count fraction
        lag(ii,2) = (sum(id>0)/nel);
    end
end

% Vertical, Scanning in the Y direction in X order
if scanDir(2) == 1
    for ii = 1:size(elMids,1)
        y = elMids(ii,3);% current Y coordinate
        x = elMids(elMids(:,3) == y,2);% x coordinates of in-line elements
        % Find number of in-line elements
        nel = size(x,1);
        % Check if solo element
        if nel == 1
            % Solo elements have no inter-element distance
            continue
        end
        % Find inter-element distances
        z = circshift(x,1);
        id = x(2:end) - z(2:end);
        % Maximum inter-element distance
        dMax = max(id);
        % Check for overlapping elements
        if dMax == 0
            % Overlapping elements have no inter-element distance
            continue
        end
        % Find gap count fraction
        lag(ii,2) = (sum(id>0)/nel);
    end
end

%% Modify Element Stiffness

% Re-calculate weights with inter-element lag
for ii = 1:size(lag,1)
    if lag(ii,2) ~= 0% only for groups of elements
        if weight(ii,2) ~= 0% combine orientation and inter-element weights
            weight(ii,2) = (weight(ii,2) + lag(ii,2))/weight(ii,2);
        else% only use inter-element weights
            weight(ii,2) = lag(ii,2);
        end
    end
end

% Re-size the weights according to maximum reduction
weight(:,2) = weight(:,2) * maxRed;

% Re-assign stiffness values
for ii = 1:size(prop,1)
    prop(weight(ii,1),1) = prop(weight(ii,1),1)*(1 - weight(ii,2));
end

end
