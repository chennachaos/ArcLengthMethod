% create and assemble element matrices

function [Kglobal, Fglobal] = get_global_matrix_vector(elemData, elemConn, coords, LM, td, disp, veloCur, acceCur, bf)

neq = max(size(disp))
nelem = size(elemConn)(1)

Kglobal = zeros(neq,neq);
Fglobal = zeros(neq,1);

for e = 1:nelem
    [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, td, disp, veloCur, acceCur, bf);

    Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
    Fglobal = Assembly_Vector(Fglobal,Flocal,LM,e);
end

endfunction

