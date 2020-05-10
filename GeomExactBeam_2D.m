function [Klocal, Flocal]=GeomExactBeam_2D(elmDat, IEN, e, XX, soln, bf)
%%% Shape Function Routine for a 1D Lagrange polynomials
af = 1.0;

p = 1;
nlocal = 2;
ndof   = 3;
nsize  = 6;

rho  = elmDat(2);
A    = elmDat(3);
I    = elmDat(4);
E    = elmDat(5);
nu   = elmDat(6);
kappa= elmDat(7);

G  = E/2.0/(1.0+nu);
EA = E*A;
EI = E*I;
GA = G*A*kappa;


Klocal=zeros(nsize,nsize); % Local stiffness matrix
Flocal=zeros(nsize,1);   % Local load vector

x0(1) = XX(IEN(e,3),1);
y0(1) = XX(IEN(e,3),2);
x0(2) = XX(IEN(e,4),1);
y0(2) = XX(IEN(e,4),2);

dx = x0(2) - x0(1);
dy = y0(2) - y0(1);
h0 = sqrt(dx*dx+dy*dy);

cth0 = dx/h0;
sth0 = dy/h0;
    
RotMat=zeros(6,6);

RotMat(1,1) =  cth0; RotMat(1,2) = -sth0;
RotMat(2,1) =  sth0; RotMat(2,2) =  cth0;
RotMat(3,3) =  1.0;
RotMat(4,4) =  cth0; RotMat(4,5) = -sth0;
RotMat(5,4) =  sth0; RotMat(5,5) =  cth0;
RotMat(6,6) =  1.0;


uxn=[0.0,0.0];
uzn=[0.0,0.0];
btn=[0.0,0.0];

res = zeros(3,1);
B=zeros(3,6);
D=zeros(3,3);

uxn(1) = soln(ndof*(IEN(e,3)-1)+1);
uzn(1) = soln(ndof*(IEN(e,3)-1)+2);
btn(1) = soln(ndof*(IEN(e,3)-1)+3);

uxn(2) = soln(ndof*(IEN(e,4)-1)+1);
uzn(2) = soln(ndof*(IEN(e,4)-1)+2);
btn(2) = soln(ndof*(IEN(e,4)-1)+3);

dummy = RotMat'*[uxn(1); uzn(1); btn(1); uxn(2); uzn(2); btn(2)];
uxn(1) = dummy(1);
uzn(1) = dummy(2);
btn(1) = dummy(3);
uxn(2) = dummy(4);
uzn(2) = dummy(5);
btn(2) = dummy(6);


nGP = 1;
[gpvec, gwvec] = get_Gauss_points(nGP);

for gp = 1:nGP
    [N,dN_dx,d2N_dx2,J,xcoord]=shape_functions_Lagrange_1D(IEN(e,3:end), XX, p, gpvec(gp));

    ux = 0.0;uz =0.0; bt = 0.0;
    dux = 0.0;duz = 0.0; dbt = 0.0;

    for ii=1:nlocal
        ux  = ux  + uxn(ii) * N(ii);
        uz  = uz  + uzn(ii) * N(ii);
        bt  = bt  + btn(ii) * N(ii);
        dux = dux + uxn(ii) * dN_dx(ii);
        duz = duz + uzn(ii) * dN_dx(ii);
        dbt = dbt + btn(ii) * dN_dx(ii);
    end

    sbt = sin(bt);
    cbt = cos(bt);

    %compute average normal strain, shear strain and curvature

    fact = (1.0+dux)*cbt - duz*sbt;

    E = dux + 0.5*(dux*dux + duz*duz);
    G = (1.0+dux)*sbt + duz*cbt;
    K = dbt * fact;

    % compute material response (elastic)

    NF = EA * E; % normal force
    SF = GA * G; % shear force
    BM = EI * K; % bending moment
            
    % multiply with volume element
            
    dvol  = J*gwvec(gp);
    fact1 = dvol * af;
    NF    = NF * dvol;
    SF    = SF * dvol;
    BM    = BM * dvol;
    EAdv  = EA * fact1;
    GAdv  = GA * fact1;
    EIdv  = EI * fact1;


    B(1,1) = (1.0+dux) * dN_dx(1);
    B(1,2) = duz * dN_dx(1);
    B(1,3) = 0.0;

    B(2,1) = sbt * dN_dx(1);
    B(2,2) = cbt * dN_dx(1);
    B(2,3) = fact * N(1);

    B(3,1) = dbt*cbt * dN_dx(1);
    B(3,2) = - dbt*sbt * dN_dx(1);
    B(3,3) = fact * dN_dx(1) - G*dbt * N(1);

    B(1,4) = (1.0+dux) * dN_dx(2);
    B(1,5) = duz * dN_dx(2);
    B(1,6) = 0.0;

    B(2,4) = sbt * dN_dx(2);
    B(2,5) = cbt * dN_dx(2);
    B(2,6) = fact * N(2);

    B(3,4) = dbt*cbt * dN_dx(2);
    B(3,5) = - dbt*sbt * dN_dx(2);
    B(3,6) = fact * dN_dx(2) - G*dbt * N(2);

    D(1,1) = EAdv;
    D(2,2) = GAdv;
    D(3,3) = EIdv;

    res(1) = NF;
    res(2) = SF;
    res(3) = BM;

    Klocal = Klocal + ( B'*D*B );
    Flocal = Flocal - ( B'*res );

    fact1 = (+ SF * cbt - BM * dbt * sbt) * af;
    fact2 = (- SF * sbt - BM * dbt * cbt) * af;

    for ii=1:nlocal
        TI   =  3*(ii-1)+1;
        TIp1 =  TI+1;
        TIp2 =  TI+2;

        for jj=1:nlocal
            TJ   = 3*(jj-1)+1;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            Klocal(TI,TJ)     =  Klocal(TI,TJ)     + dN_dx(ii)*NF*dN_dx(jj) * af;
            Klocal(TIp1,TJp1) =  Klocal(TIp1,TJp1) + dN_dx(ii)*NF*dN_dx(jj) * af;

            fact3 =  dN_dx(ii)*BM*cbt*dN_dx(jj) * af;
            fact4 = -dN_dx(ii)*BM*sbt*dN_dx(jj) * af;

            Klocal(TI,TJp2)   = Klocal(TI,TJp2)    + (fact3 + dN_dx(ii)*fact1*N(jj) );
            Klocal(TIp1,TJp2) = Klocal(TIp1,TJp2)  + (fact4 + dN_dx(ii)*fact2*N(jj) );
            Klocal(TIp2,TJ)   = Klocal(TIp2,TJ)    + (fact3 + N(ii)*fact1*dN_dx(jj) );
            Klocal(TIp2,TJp1) = Klocal(TIp2,TJp1)  + (fact4 + N(ii)*fact2*dN_dx(jj) );

            Klocal(TIp2,TJp2) = Klocal(TIp2,TJp2)  + (N(ii)*(-SF*G-BM*dbt*fact)*N(jj) - dN_dx(ii)*BM*G*N(jj) - N(ii)*BM*G*dN_dx(jj)) * af;
        end
    end
end

h = XX(IEN(e,4)) - XX(IEN(e,3));

Flocal(1) = Flocal(1) + 0.5*h*bf(1);
Flocal(4) = Flocal(4) + 0.5*h*bf(1);

Flocal(2) = Flocal(2) + 0.5*h*bf(2);
Flocal(5) = Flocal(5) + 0.5*h*bf(2);

Flocal = RotMat*Flocal;
Klocal = RotMat*Klocal*RotMat';

