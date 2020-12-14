function [ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)
clc
% coords: global coordinates of the nodes, x, y, and z
% elemConn: element connectivities
% nelem: total number of elements
% nnode: total number of nodes
% nperelem: number of nodes per element for all elements
% ndof: number of degrees of per node

fid=fopen(fname,'r');

% ndim

line=fgets(fid);
linestr = strsplit(line, ",");
ndim = int32(str2num(linestr{1,2}))

% ndof

line=fgets(fid);
linestr = strsplit(line, ",");
ndof = int32(str2num(linestr{1,2}))

% nodes

line=fgets(fid);
linestr = strsplit(line, ",");
nnode = int32(str2num(linestr{1,2}))

nperelem = 2;
nsize = nperelem*ndof;
neq = nnode*ndof;
%if(arclen)
%  neq = neq+1;
%end
%neq

coords = zeros(nnode,ndim);
for i=1:nnode
    line = fgets(fid);
    linestr = strsplit(line, ",");

    coords(i,1) = double(str2num(linestr{1,2}));
    coords(i,2) = double(str2num(linestr{1,3}));
    if(ndim == 3)
      coords(i,3) = double(str2num(linestr{1,4}));
    end
end    


% element data

line=fgets(fid);
linestr = strsplit(line, ",");
nelemData = int32(str2num(linestr{1,2}))

elemData = zeros(nelemData, 10);
for i=1:nelemData
    line = fgets(fid);
    linestr = strsplit(line, ",");

    for j=1:10
      elemData(i,j) = double(str2num(linestr{1,j+1}));
    end
end

% elements

line=fgets(fid);
linestr = strsplit(line, ",");
nelem = int32(str2num(linestr{1,2}))

elemConn = zeros(nelem, 4, "int32");
for i=1:nelem
    line = fgets(fid);
    linestr = strsplit(line, ",");

    elemConn(i,1) = int32(str2num(linestr{1,2}));
    elemConn(i,2) = int32(str2num(linestr{1,3}));
    elemConn(i,3) = int32(str2num(linestr{1,4}));
    elemConn(i,4) = int32(str2num(linestr{1,5}));
end

% Dirichlet boundary conditions

line=fgets(fid);
linestr = strsplit(line, ",");
nDBC    = int32(str2num(linestr{1,2}))

%dbclist = zeros(nDBC, 3);
dbcnodes = zeros(nDBC, 1, "int32");
for i=1:nDBC
    line = fgets(fid);
    linestr = strsplit(line, ",");

    n1 = int32(str2num(linestr{1,1}));
    n2 = int32(str2num(linestr{1,2}));
%    dbclist(i,3) = double(str2num(linestr{1,3}));
    dbcnodes(i) = (n1-1)*ndof+n2;
end

assy4r = setdiff([1:neq], dbcnodes)';

% Force boundary conditions

line=fgets(fid);
linestr = strsplit(line, ",");
nFBC    = int32(str2num(linestr{1,2}))
fbclist = zeros(nFBC, 3);

dof_force = zeros(nFBC, 1, "int32");
Fext = zeros(neq, 1);
for i=1:nFBC
    line = fgets(fid);
    linestr = strsplit(line, ",");

    n1 = int32(str2num(linestr{1,1}));
    n2 = int32(str2num(linestr{1,2}));
    ind = (n1-1)*ndof + n2;

    dof_force(i) = ind;
    Fext(ind) = double(str2num(linestr{1,3}));
end

% for output

line=fgets(fid);
linestr = strsplit(line, ",");
nOutput = int32(str2num(linestr{1,2}))
outputlist = zeros(nOutput, 1, "int32");

for i=1:nOutput
    line = fgets(fid);
    linestr = strsplit(line, ",");

    n1 = int32(str2num(linestr{1,1}));
    n2 = int32(str2num(linestr{1,2}));

    outputlist(i,1) = (n1-1)*ndof+n2;
end


% Arclength parameters

line=fgets(fid);
linestr = strsplit(line, ",");
arclen  = (int32(linestr{1,2}) == 1);

line=fgets(fid)
linestr = strsplit(line, ",");
maxloadSteps = int32(str2num(linestr{1,1}));

line=fgets(fid)
linestr = strsplit(line, ",");
loadincr = double(str2num(linestr{1,1}));

fclose(fid);

% data structures

LM  = zeros(nelem, nsize);

for e=1:nelem
    count = 1;
    for jj=1:nperelem
        ind = ndof*(elemConn(e,jj+2)-1)+1;
        for kk=1:ndof
            LM(e,count) = ind;
            ind = ind + 1;
            count = count + 1;
        end
    end
    count = count - 1;
end



