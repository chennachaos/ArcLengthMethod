function [Klocal, Flocal]=Truss_2D_model1(elmDat, IEN, e, XX, soln, bf)

ndof = 2;

finite = (int32(elmDat(IEN(e,1), 1)) == 1) ;
rho0 = elmDat(IEN(e,1),2);
A0   = elmDat(IEN(e,1),3);
E    = elmDat(IEN(e,1),4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp  = zeros(4,1);

X1 = XX(IEN(e,3),1);
Y1 = XX(IEN(e,3),2);
X2 = XX(IEN(e,4),1);
Y2 = XX(IEN(e,4),2);

disp(1) = soln(ndof*(IEN(e,3)-1)+1);
disp(2) = soln(ndof*(IEN(e,3)-1)+2);
disp(3) = soln(ndof*(IEN(e,4)-1)+1);
disp(4) = soln(ndof*(IEN(e,4)-1)+2);


% compute the orientation of the element

    dx = X2 - X1;
    dy = Y2 - Y1;
    
    L0 = sqrt(dx*dx+dy*dy);

    x1 = X1 + disp(1);
    y1 = Y1 + disp(2);
    x2 = X2 + disp(3);
    y2 = Y2 + disp(4);

    dx = x2 - x1;
    dy = y2 - y1;
    
    L = sqrt(dx*dx+dy*dy);

    B0 = zeros(4,1);
    B(1) = -dx;   B(2) = -dy;   B(3) =  dx;   B(4) =  dy;

    H = zeros(4,4);
    H(1,1) =  1.0;    H(1,3) = -1.0;
    H(2,2) =  1.0;    H(2,4) = -1.0;
    H(3,1) = -1.0;    H(3,3) =  1.0;
    H(4,2) = -1.0;    H(4,4) =  1.0;

    Klocal = zeros(4,4);
    Flocal = zeros(4,1);

    strain = (L/L0 - 1.0) ;

    % axial force
    N = (A0*E)*strain;

    % residual
    Flocal = Flocal - N*B'/L;

    % stiffness
    Klocal = Klocal + ((E*A0*B')*B)/L^3 ;
    Klocal = Klocal + ((N/L)*H);
end



