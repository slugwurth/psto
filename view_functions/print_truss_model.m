function print_truss_model(fid,nnd,nel,nne,nodof,eldof,nodalCoords,elConnec,prop,nf,load,n)
%% Write Output File: Truss Model Data
% 		also works for frame elements
%
% written by Amar Khennane 2013
% as published in ISBN 1466580208 
% modified by Matt Ireland 2019
%
% Accepts 
%   the file id of an open file
%   the number of nodes
%   the number of elements
%   the number of nodes per element
%   the number of nodal degrees of freedom
%   the number of element degrees of freedom
%   the global nodal coordinate matrix
%   the global element connectivity matrix
%   the global element property matrix
%   the global nodal freedom matrix
%   the global imposed load matrix
%   the number of free DOF
%
% Returns 
%   Nothing. Writes to file.

fprintf(fid, ' ******* PRINTING MODEL DATA **************\n\n\n');

% Print Nodal coordinates
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Number of nodes:                                 %g\n', nnd );
fprintf(fid, 'Number of elements:                              %g\n', nel );
fprintf(fid, 'Number of nodes per element:                     %g\n', nne );
fprintf(fid, 'Number of degrees of freedom per node:           %g\n', nodof);
fprintf(fid, 'Number of degrees of freedom per element:        %g\n\n\n', eldof);

fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node        X            Y \n');
for i=1:nnd
    fprintf(fid,' %g,      %07.2f,      %07.2f\n',i, nodalCoords(i,1), nodalCoords(i,2));
end
fprintf(fid,'\n');

% Print element connectivity
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Element      Node_1      Node_2 \n');
for i=1:nel
    fprintf(fid,'    %g,         %g,          %g\n',i, elConnec(i,1), elConnec(i,2));
end
fprintf(fid,'\n');

% Print element property
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Element        E             A \n');
for i=1:nel
    fprintf(fid,'    %g,       %g,       %g\n',i, prop(i,1), prop(i,2));
end
fprintf(fid,'\n');

% Print Nodal freedom
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node      disp_U     disp_V\n');
for i=1:nnd
    fprintf(fid,'  %g,        %g,          %g\n',i, nf(i,1), nf(i,2));
end
fprintf(fid,'\n');

% Print Nodal loads
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, 'Node      load_X         load_Y\n');
for i=1:nnd
    fprintf(fid,'  %g,      %07.2f,        %07.2f\n',i, load(i,1), load(i,2));
end

fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid,'\n');
fprintf(fid,'Total number of active degrees of freedom, n = %g\n',n);
fprintf(fid,'\n');

end
