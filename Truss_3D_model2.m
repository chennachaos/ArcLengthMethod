function [Klocal, Flocal]=Truss_3D_model2(elmDat, IEN, e, XX, td, soln, solnPrev, veloCur, acceCur, bf)

ndof = 3;

af = td(2);
am = td(1);
d1 = td(5);


finite = (int32(elmDat(1)) == 1) ;
rho0 = elmDat(2);
A0   = elmDat(3);
E    = elmDat(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp  = zeros(6,1);
dispPrev  = zeros(6,1);
accC  = zeros(6,1);
velC  = zeros(6,1);
B0    = zeros(6,1);

H  = zeros(6,6);


% rotate nodal displacements and compute nodal positions on element axis   

X1 = XX(IEN(e,3),1);
Y1 = XX(IEN(e,3),2);
Z1 = XX(IEN(e,3),3);
X2 = XX(IEN(e,4),1);
Y2 = XX(IEN(e,4),2);
Z2 = XX(IEN(e,4),3);

disp(1) = soln(ndof*(IEN(e,3)-1)+1);
disp(2) = soln(ndof*(IEN(e,3)-1)+2);
disp(3) = soln(ndof*(IEN(e,3)-1)+3);
disp(4) = soln(ndof*(IEN(e,4)-1)+1);
disp(5) = soln(ndof*(IEN(e,4)-1)+2);
disp(6) = soln(ndof*(IEN(e,4)-1)+3);

dispPrev(1) = solnPrev(ndof*(IEN(e,3)-1)+1);
dispPrev(2) = solnPrev(ndof*(IEN(e,3)-1)+2);
dispPrev(3) = solnPrev(ndof*(IEN(e,3)-1)+3);
dispPrev(4) = solnPrev(ndof*(IEN(e,4)-1)+1);
dispPrev(5) = solnPrev(ndof*(IEN(e,4)-1)+2);
dispPrev(6) = solnPrev(ndof*(IEN(e,4)-1)+3);

dispC = af*disp + (1.0-af)*dispPrev;


% compute the orientation of the element

    dx = X2 - X1;
    dy = Y2 - Y1;
    dz = Z2 - Z1;
    
    L0 = sqrt(dx*dx+dy*dy+dz*dz);

    x1 = X1 + dispC(1);
    y1 = Y1 + dispC(2);
    z1 = Z1 + dispC(3);
    
    x2 = X2 + dispC(4);
    y2 = Y2 + dispC(5);
    z2 = Z2 + dispC(6);
    
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;
    
    L = sqrt(dx*dx+dy*dy+dz*dz);

    B(1) = -dx;   B(2) = -dy;   B(3) = -dz;
    B(4) =  dx;   B(5) =  dy;   B(6) =  dz;
    B = B/L0/L0;

    H(1,1) =  1.0;    H(1,4) = -1.0;
    H(2,2) =  1.0;    H(2,5) = -1.0;
    H(3,3) =  1.0;    H(3,6) = -1.0;

    H(4,1) = -1.0;    H(4,4) =  1.0;
    H(5,2) = -1.0;    H(5,5) =  1.0;
    H(6,3) = -1.0;    H(6,6) =  1.0;

    type = 1;

if type == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% type 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Klocal = zeros(6,6);
    Mlocal = zeros(6,6);
    Flocal = zeros(6,1);

    strain = (L*L/L0/L0 - 1.0)/2.0 ;

    % axial force
    N = (A0*E)*strain;

    % residual
    Flocal = Flocal - N*L0*B';

    % stiffness
    Klocal = Klocal + ((E*A0*L0*B')*B) ;
    Klocal = Klocal + ((N/L0)*H);

    Klocal = Klocal*af;
end



