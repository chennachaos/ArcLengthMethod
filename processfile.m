function [nnode, nelem, coords, elemConn, elemData, ndof, nDBC, dbclist, nFBC, fbclist, LM, neq, assy4r, dof_force, Fext] = processfile(fname)
clc
% coords: global coordinates of the nodes, x, y, and z
% elemConn: element connectivities
% nelem: total number of elements
% nnode: total number of nodes
% nperelem: number of nodes per element for all elements
% ndof: number of degrees of per node


fid=fopen(fname,'r');

line=fgets(fid);
linestr = strsplit(line, ",");
ndof = int32(str2num(linestr(2){:}))


% nodes

line=fgets(fid);
linestr = strsplit(line, ",");
nnode = int32(str2num(linestr(2){:}))
coords = zeros(nnode,2);

for i=1:nnode
    line = fgets(fid);
    linestr = strsplit(line, ",");

    coords(i,1) = double(str2num(linestr(2){:}));
    coords(i,2) = double(str2num(linestr(3){:}));
end    
coords

% elements

line=fgets(fid);
linestr = strsplit(line, ",");
nelemData = int32(str2num(linestr(2){:}))
elemData = zeros(nelemData, 10);

for i=1:nelemData
    line = fgets(fid);
    linestr = strsplit(line, ",");

    for j=1:10
      elemData(i,j) = double(str2num(linestr(j+1){:}));
    endfor
end
elemData

% elements

line=fgets(fid);
linestr = strsplit(line, ",");
nelem = int32(str2num(linestr(2){:}))
elemConn = zeros(nelem, 4, "int32");

for i=1:nelem
    line = fgets(fid);
    linestr = strsplit(line, ",");

    elemConn(i,1) = int32(str2num(linestr(2){:}));
    elemConn(i,2) = int32(str2num(linestr(3){:}));
    elemConn(i,3) = int32(str2num(linestr(4){:}));
    elemConn(i,4) = int32(str2num(linestr(5){:}));
end
elemConn

% Dirichlet boundary conditions

line=fgets(fid);
linestr = strsplit(line, ",");
nDBC    = int32(str2num(linestr(2){:}))
dbclist = zeros(nDBC, 3);

for i=1:nDBC
    line = fgets(fid);
    linestr = strsplit(line, ",");

    dbclist(i,1) = double(str2num(linestr(1){:}));
    dbclist(i,2) = double(str2num(linestr(2){:}));
    dbclist(i,3) = double(str2num(linestr(3){:}));
end
dbclist

% Force boundary conditions

line=fgets(fid);
linestr = strsplit(line, ",");
nFBC    = int32(str2num(linestr(2){:}))
fbclist = zeros(nFBC, 3);

for i=1:nFBC
    line = fgets(fid);
    linestr = strsplit(line, ",");

    fbclist(i,1) = double(str2num(linestr(1){:}));
    fbclist(i,2) = double(str2num(linestr(2){:}));
    fbclist(i,3) = double(str2num(linestr(3){:}));
end
fbclist

% Whether arclength is active or not
line=fgets(fid);
linestr = strsplit(line, ",");
arclen  = (linestr(2){:} == "ON")


fclose(fid);

% data structures

nperelem = 2;
nsize = nperelem*ndof;

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


dbcnodes = []
for i=1:nDBC
  n1 = dbclist(i,1);
  n2 = dbclist(i,2);
  ind = (n1-1)*ndof;
  dbcnodes = [dbcnodes ind+n2];
endfor

neq = nnode*ndof;
if(arclen)
  neq = neq+1;
endif
neq

assy4r = setdiff([1:neq], dbcnodes)';

dof_force = zeros(nFBC, 1, "int32");
Fext = zeros(neq, 1);
for i=1:nFBC
  n1 = fbclist(i,1);
  n2 = fbclist(i,2);
  ind = (n1-1)*ndof+n2;

  dof_force(i) = (fbclist(i,1)-1)*ndof + fbclist(i,2);
  Fext(ind) = fbclist(i,3);
endfor



