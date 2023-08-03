function [vertices,faces] = readOff(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Name: readOffFiles
% Description: The function loads .off files into MATLAB
%
% Inputs: 
%  - mesh object
%  - filename
% 
% Outputs: 
%  - vertices - matrix of [N vertices, 3]
%  - faces - matrix of [N faces, Number of vertices in a face]
%       Number of vertices in a face should be 3 or 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fid = fopen(filename,'r');
if( fid==-1 )
	error('Can''t open the file.');
	return;
end

% Reading the header:
firstLine = fgets(fid);
if ~strcmp(firstLine(1:3), 'OFF')
    error('The file is not a valid OFF file.');
    return;
end

% Reading the fisrt line: 
% N is a vector containing the number of vertices, ...
%       number of faces, number of edges 
[N,cnt] = fscanf(fid,'%d %d %d', 3);
if cnt~=3
    warning('Problem reading file.');
end
nv = N(1);
nf = N(2);


% Reading the vertices 
Dim = 3;
[vertices,cnt] = fscanf(fid,'%f %f %f', Dim*nv);
if cnt~=Dim*nv
    warning('Problem reading file.');
end
vertices = reshape(vertices, Dim,  nv)';


% Reading the faces:
% Reading the first line first to check if the mesh is triangles of quad
[faces,cnt] = fscanf(fid,'%d %d %d %d', 4*nf);
if cnt~=4*nf
    warning('Problem reading file.');
end

faces = reshape(faces, 4, nf)';
faces = faces(:,2:end)+1;

fclose(fid);

end


