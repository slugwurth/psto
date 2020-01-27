function Plot2DGeometryUndeformed(varargin)
%% Plot Undeformed 2D Geometry
% ********************************************************************
%            Finite Element Learning Toolbox (FELT)
%                       Senthil Vel
%					modified by Matt Ireland 2019
%
%                Function: Plot2DGeometryUndeformed
%       Plot the undeformed 2D geometry (plane stress or plane strain)
%
%  Input:  ENODES - Array of element connectivity
%          NODALCOORDINATES - Array containing the node numbers and the
%                             corresponding x and y coordinates
%          numToggle - 0 for no numbers, 1 to plot node and element numbers
%          props - a matrix of element stiffnesses
%  Output: None.  The undeformed shape will be plotted in the current
%          figure window, either with element numbering or normalized
%          element stiffnesses superimposed on the elements.
%*********************************************************************

if nargin == 5
    [elConnec,nodalCoords,numToggle,caption,props] = varargin{1:5};
    showPropOverlay = 1;
else
    [elConnec,nodalCoords,numToggle,caption] = varargin{1:4};
    showPropOverlay = 0;
end

% Clear the current figure and hold the graph
hold on

% Determine the total number of elements and nodes per element
[rows,cols] = size(elConnec);
nel = rows;
nne = cols-1;

% Draw the elements by connecting the nodes using straight lines
for e = 1:nel
    x = [];
    y = [];
    for n = 1:nne 
       node = elConnec(e,1+n);
       x = [x, nodalCoords(nodalCoords(:,1) == node,2)];
       y = [y, nodalCoords(nodalCoords(:,1) == node,3)];    
    end
    xmid = mean(x);
    ymid = mean(y);
    x = [x,x(1)];
    y = [y,y(1)];
    H = line(x,y);
    set(H,'LineWidth',0.5);        % Set line thickness
    set(H,'Color',[0.3 0.75 0.93]);             % Draw in grey color
%     set(H,'Color',[1 1 1]);
    
    % Insert a dot at each node
    H = plot(x,y,'k.');
%     H = plot(x,y,'w');
    set(H,'MarkerSize',20);             % Increase the size of the dots
    %set(H,'MarkerSize',12);             % Increase the size of the dots    
    
    if numToggle == 1
        % Insert the element number near the center of each element
        %str = ['\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {',num2str(e),'}}}'];
        if showPropOverlay == 1
            str = num2str(props(e),3);
        else
            str = num2str(e);
        end
        H = text(xmid,ymid,str,'Interpreter', 'latex');
        set(H,'FontSize',12);            %set font size
        set(H,'FontWeight','bold')       %set to bold font
        set(H,'Color',[180 0 0]/255);    % set font color to red
    end
end

% Enlarge axes limits to include some additional empty space on the sides
axis tight;
H = gca;
XLim =get(H,'XLim');
YLim =get(H,'YLim');
XRange = XLim(2)-XLim(1);
YRange = YLim(2)-YLim(1);
set(H,'XLim',[XLim(1)-0.1*XRange, XLim(2)+0.1*XRange]);
set(H,'YLim',[YLim(1)-0.1*YRange, YLim(2)+0.1*YRange]);

% Insert figure caption
title(caption);

if numToggle == 1
    % Sliglty offset the text from the selected coordinates by an amount Delta
    Delta = max([XRange, YRange])/200;
    
    % Insert node numbers
    for e = 1: nel
        for n = 1: nne
            node = elConnec(e,1+n);
            H = text(nodalCoords(nodalCoords(:,1) == node,2) + Delta, ...
                nodalCoords(nodalCoords(:,1) == node,3) + Delta, ...
                num2str(node));
            set(H,'FontSize',12);            %set font size
            set(H,'FontWeight','bold')       %set to bold font
            set(H,'Color',[0 140 0]/255);    %set font color to green
        end
    end
end
% Set aspect ratio
H=gca;
set(H,'DataAspectRatio',[1 1 1]);
set(H,'PlotBoxAspectRatioMode','manual')
set(H,'PlotBoxAspectRatio',[1 1 1])

end
