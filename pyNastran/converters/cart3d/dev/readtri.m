function [t, p, id] = readtri( fname )
%READTRI reads a Cart3D *.tri file into Matlab triangle data
%
%   [t, p, id] = readtri( fname ) Reads in a Cart3D *.tri file into the
%   connectivity matrix 't', the point matrix 'p', and the surface id
%   vector 'id'.
%
%   If the file does not contain surface id's, NaN is returned in id.
%
%   See also TRISURF, TRIMESH, TRIMESH.


%   Rob McDonald 
%   ramcdona@calpoly.edu  
%   8 October 2012 v. 1.0

fp = fopen(fname,'r');

% Read in number of points and triangles.
[n] = fscanf(fp, '%d',2);
npt = n(1);
ntri = n(2);

% Read in the point coordinates.
p = fscanf(fp, '%f', [3 npt]);
% plot3(p(1,:),p(2,:),p(3,:),'x')
% axis equal

% Read in the triangle connectivity.
t = fscanf(fp, '%d', [3 ntri])';
% trisurf(t, p(1,:), p(2,:), p(3,:));
% axis equal

% Read in the surface id's if they exist.
id = fscanf(fp, '%d', ntri);

% Check if eof occurred during previous read.
if(feof(fp))
  id = nan;
end

% trisurf(t, p(1,:), p(2,:), p(3,:), id);
% axis equal

fclose(fp);