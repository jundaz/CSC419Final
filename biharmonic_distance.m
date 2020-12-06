function dists = biharmonic_distance(off_file, dist_file, idx_file)
%   Computes pairwise biharmonic distances between points on the surface.
%   This is unoptimized code that demonstrates how the formalism of the paper
%   can be implemented. Will not run on large meshes.
%   
%   dists = biharmonic_distance(off_file, dist_file, idx_file);
%
%   off_file: input mesh file in off format
%
%   dist_file: the output distance file. The distances are saved in the file 
%   in the following order (1,2), (1,3), ..., (1,n), (2,3), ..., (2,n), ...,
%   ..., (n-1,n)
%
%   idx_file: the file with vertex indices between which the distance is to be
%   computed. This argument is optional; if not supplied distances between all
%   mesh vertices are computed.
%
%   Raif Rustamov, Drew University

%read mesh from the off file
fid=fopen(off_file);
fgetl(fid);
nos = fscanf(fid, '%d %d  %d', [3 1]);
nopts = nos(1);
notrg = nos(2);
coord = fscanf(fid, '%g %g  %g', [3 nopts]);
coord = coord';
triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1; %%we have added 1 because the vertex indices start from 0 in off format
fclose(fid);

tic;
bir=[1 1 1];
trgarea=zeros(1,notrg);
for i=1:notrg
    trgx=[coord(triang(i,1),1) coord(triang(i,2),1) coord(triang(i,3),1)];
    trgy=[coord(triang(i,1),2) coord(triang(i,2),2) coord(triang(i,3),2)];
    trgz=[coord(triang(i,1),3) coord(triang(i,2),3) coord(triang(i,3),3)];
    aa=[trgx; trgy; bir];
    bb=[trgx; trgz; bir];
    cc=[trgy; trgz; bir];
    area=sqrt(det(aa)^2+det(bb)^2+det(cc)^2)/2;
    trgarea(i)=area;
end


%find the approximate voronoi area of each vertex
AM = zeros(nopts, 1);
for i=1:notrg
    AM(triang(i,1:3)) = AM(triang(i,1:3)) + trgarea(i)/3;
end

A = sparse(nopts, nopts);
%now construct the cotan laplacian
for i=1:notrg
    for ii=1:3
        for jj=(ii+1):3
            kk = 6 - ii - jj; % third vertex no
            v1 = triang(i,ii);
            v2 = triang(i,jj);
            v3 = triang(i,kk);
            e1 = [coord(v1,1) coord(v1,2) coord(v1,3)] - [coord(v2,1) coord(v2,2) coord(v2,3)];
            e2 = [coord(v2,1) coord(v2,2) coord(v2,3)] - [coord(v3,1) coord(v3,2) coord(v3,3)];
            e3 = [coord(v1,1) coord(v1,2) coord(v1,3)] - [coord(v3,1) coord(v3,2) coord(v3,3)];
            cosa = e2* e3'/sqrt(sum(e2.^2)*sum(e3.^2));
            sina = sqrt(1 - cosa^2);
            cota = cosa/sina;
            w = 0.5*cota;
            A(v1, v1) = A(v1, v1) - w;
            A(v1, v2) = A(v1, v2) + w;
            A(v2, v2) = A(v2, v2) - w;
            A(v2, v1) = A(v2, v1) + w;
        end
    end
end
A = -A;
T = sparse([1:nopts], [1:nopts], (AM), nopts, nopts, nopts);
A = A*inv(T)*A;

%make the first row&col of A zeros, this forces the first entry of
%the solution vector to be zero; takes care of inf many solns
A(1, :) = zeros(1,nopts);
A(:, 1) = zeros(nopts,1);
A(1, 1) = 1;


%read the file with indices of points to find pairwise distances between
if nargin < 3
    sampleidx = (1:nopts)'; %%if no argument for index file, compute all pairs
else
    fid=fopen(idx_file);
    sampleidx = fscanf(fid, '%d', [1 inf]);
    sampleidx=sampleidx + 1; %%we have added 1 because the vertex indices start from 0
    fclose(fid);
end
samplesize = length(sampleidx);

%make the right hand side vector(s) for the system
h = eye(nopts) - (1/nopts)* ones(nopts,nopts);
h(1,:) = zeros(1, nopts);
h = h(1:nopts,sampleidx);

%solve
g = full(A\h);
g = g - repmat(sum(g)/nopts,nopts,1); %%shift the solution so that add up to zero, but weighted
g = g(sampleidx,:);

%compute dists; this is a vectorization of d^2 = g(i,i)+g(j,j)-2g(i,j)
v1 = diag(g)';
v2 = ones(samplesize,1);
m = v2*v1;
dists = sqrt(m + m' - 2*g);

%save the distances in to file
%distances are saved in the following order (1,2), (1,3), ..., (1,n), (2,3), ..., (2,n), ..., ..., (n-1,n)
savedists = squareform(dists); %convert to a vector
fid = fopen(dist_file,'w');
fprintf(fid,'%d\n',savedists);
fclose(fid);

%save in vtk format for visual inspection
%use Paraview to view this format
source_id = 1;
ofid = fopen('control_biharmonic.vtk','w');
fprintf(ofid, '# vtk DataFile Version 3.0\n');
fprintf(ofid,'vtk output\n');
fprintf(ofid,'ASCII\n');
fprintf(ofid,'DATASET POLYDATA\n');
fprintf(ofid,'POINTS %d float\n', nopts);
fprintf(ofid,'%g %g %g\n', coord');
fprintf(ofid,'POLYGONS %d %d\n', notrg, 4*notrg)
fprintf(ofid,'3 %d %d %d\n', triang'-1);
fprintf(ofid,'\n');
fprintf(ofid,'POINT_DATA %d\n', nopts);
fprintf(ofid,'SCALARS distance_from float\n');
fprintf(ofid,'LOOKUP_TABLE default\n');
fprintf(ofid,'%g\n', dists(:,source_id));
fclose(ofid);

