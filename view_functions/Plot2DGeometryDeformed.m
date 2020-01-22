function Plot2DGeometryDeformed(elConnec,nodalCoords,NDOF,d,SFactor)
%% Plot Deformed 2D Geometry
%*********************************************************************
%            Finite Element Learning Toolbox (FELT)
%                       Senthil Vel
%					modified by Matt Ireland 2019
%
%                Function: Plot2DGeometryDeformed
%            Plot the deformed 2D geometry
%
%  Input:  ENODES - Array of element connectivity
%          NODALCOORDINATES - Array containing the node numbers and the
%                             corresponding x and y coordinates
%          d - Solution array for the nodal displacements 
%          SFactor - Scale (i.e., amplification) factor for plotting the
%                    the displacements
%  Output: None.  The deformed shape will be plotted in the current
%          figure window.
%*********************************************************************

% First plot the undeformed shape
Plot2DGeometryUndeformed(elConnec,nodalCoords,0,'Undeformed');
hold on;

% Determine the number of elements and nodes per element
[rows,cols] = size(elConnec);
nel = rows;
nne = cols-1;

% Draw the deformed elements by connecting the displaced nodes
for e = 1: nel
    x = [];
    y = [];
    for n = 1: nne 
       node = elConnec(e,1+n);
       x = [x, nodalCoords(nodalCoords(:,1) == node,2) + ...
           SFactor * d(nodalCoords(:,1) == node,1)];
       y = [y, nodalCoords(nodalCoords(:,1) == node,3) + ...
           SFactor * d(nodalCoords(:,1) == node,2)];    
    end
    x = [x,x(1)];
    y = [y,y(1)];
    H = line(x,y);
    set(H,'LineWidth',1.5);        % Set line thickness
    set(H,'Color','k');             % Draw in grey color
     
    % Insert a dot at each node
    H = plot(x,y,'k.');
    set(H,'MarkerSize',20);             % Increase the size of the dot   

end

% Enlarge axes limits to include some additional empty space on the sides
axis tight;
H = gca;
XLim = get(H,'XLim');
YLim = get(H,'YLim');
XRange = XLim(2) - XLim(1);
YRange = YLim(2) - YLim(1);
set(H,'XLim',[XLim(1)-0.1*XRange, XLim(2)+0.1*XRange]);
set(H,'YLim',[YLim(1)-0.1*YRange, YLim(2)+0.1*YRange]);

% Set aspect ratio
H = gca;
set(H,'DataAspectRatio',[1 1 1]);
set(H,'PlotBoxAspectRatioMode','manual')
set(H,'PlotBoxAspectRatio',[1 1 1])

% Insert figure caption
title('Deformed shape');
end
