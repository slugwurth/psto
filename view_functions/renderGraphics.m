function [s1,s2,s3,s4,f1,f2,f3] = renderGraphics(elDist,node_disp,nodalCoords,scan,prop,E)
%% Render Graphics
%
%	Matt Ireland 2019
%
% Accepts 
%   the global element connectivity matrix, elConnec
%   the global nodal displacement matrix, node_disp
%   the global nodal coordinate matrix, nodalCoords
%   the particle graph representation, graph
%   a logical scan value, scan
%   the global element property matrix, prop
%   the nominal element modulus, E
%   
% Returns 
%   the deflection in x render, s1
%   the deflection in y render, s2
%   the nodal rotation render, s3
%   the graph/network representation, s4
%   the deformed structure render, f2
%   the undeformed structure render, f3
%       with   element numbering overlay
%       or     normalized element stiffness overlay

% Calculate magnitude of nodal displacement
mag = (sum([node_disp(:,1).^2 node_disp(:,2).^2 node_disp(:,3).^2],2)).^(1/2);
% Output variables for colorbar formatting
max1 = max(node_disp(:,1)); min1 = min(node_disp(:,1));
max2 = max(node_disp(:,2)); min2 = min(node_disp(:,2));
max3 = max(node_disp(:,3)); min3 = min(node_disp(:,3));
max4 = max(mag); min4 = min(mag);
span1 = max1 - min1; 
span2 = max2 - min2; 
span3 = max3 - min3;
span4 = max4 - min4;

f1 = figure('Units','normalized','OuterPosition',[0 0 1 1]);% Nodal displacement color plots
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
s1 = subplot(2,2,1); 
    % Plot undeformed truss structure
    Plot2DGeometryUndeformed(elDist,nodalCoords,0,'X Displacement');

    % Overlay colored nodes by deflection
    hold on
    % Deflection in X
    scatter(nodalCoords(:,2),nodalCoords(:,3),150,node_disp(:,1),'filled');
    % Coloring
    colormap(s1,winter);
    c = colorbar('Ticks',(min1:span1/4:max1));
    % Colorbar label
    c.Label.String = 'Displacement in X';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s2 = subplot(2,2,2);
    % Plot undeformed truss structure
    Plot2DGeometryUndeformed(elDist,nodalCoords,0,'Y Displacement');

    % Overlay colored nodes by deflection
    hold on
    % Deflection in Y
    scatter(nodalCoords(:,2),nodalCoords(:,3),150,node_disp(:,2),'filled');
    % Coloring
    colormap(s2,summer);
    c = colorbar('Ticks',(min2:span2/4:max2));
    % Colorbar label
    c.Label.String = 'Displacement in Y';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s3 = subplot(2,2,3);
    % Plot undeformed truss structure
    Plot2DGeometryUndeformed(elDist,nodalCoords,0,'CCW Rotation');

    % Overlay colored nodes by deflection
    hold on
    % Rotation CCW in rad
    scatter(nodalCoords(:,2),nodalCoords(:,3),150,node_disp(:,3),'filled');
    % Coloring
    colormap(s3,autumn);
    c = colorbar('Ticks',(min3:span3/4:max3));
    % Colorbar label
    c.Label.String = 'Rotation CCW in rad';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
s4 = subplot(2,2,4);
    % Plot undeformed truss structure
    Plot2DGeometryUndeformed(elDist,nodalCoords,0,'Displacement Magnitude');

    % Overlay colored nodes by deflection
    hold on
    % Total magnitude of displacement
    scatter(nodalCoords(:,2),nodalCoords(:,3),150,mag,'filled');
    % Coloring
    colormap(s4,jet);
    c = colorbar('Ticks',(min4:span4/4:max4));
    % Colorbar label
    c.Label.String = 'Displacement Magnitude';

f2 = figure('Units','normalized','OuterPosition',[0 0 1 1]);% Plot deformed structure
% Normalize scale factor to maximum displacement
SFactor = 10 / (1 - max(max(node_disp)));
% Plot deformed structure
Plot2DGeometryDeformed(elDist,nodalCoords,0,node_disp,SFactor);

f3 = figure('Units','normalized','OuterPosition',[0 0 1 1]);% Plot undeformed structure
% Overlay element normalized stiffnesses if using scan modulation
if scan == 1
    % Normalize element stiffnesses
    prop = prop(:,1)./E;
    % Plot undeformed structure with normalized stiffness superimposed
    Plot2DGeometryUndeformed(elDist,nodalCoords,1,'Normalized Stiffness',prop);
end
% Regular element number overlay
if scan == 0
    % Plot undeformed structure with element number superimposed
    Plot2DGeometryUndeformed(elDist,nodalCoords,0,'Undeformed');
end

end

