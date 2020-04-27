function [N,dN_dx,d2N_dx2,J,x]=shape_functions_Lagrange_1D(NodeNums, XX, p, xi)
%%% Shape Function Routine for a 1D Lagrange polynomials

nen = p+1;

N =zeros(nen,1); 
dN_dxi   = zeros(nen,1);
d2N_dxi2 = zeros(nen,1);

if(p == 1)
    N(1) = 0.5*(1.0 - xi);
    N(2) = 0.5*(1.0 + xi);
    
    dN_dxi(1) = -0.5;
    dN_dxi(2) =  0.5;
    
    d2N_dxi2(1) = 0.0;
    d2N_dxi2(2) = 0.0;
    
elseif(p == 2)
    val = xi*xi;
    
    N(1) = -0.5 * (xi - val);
    N(2) = 1.0 - val;
    N(3) = 0.5 *(xi + val);
    
    val = 2.0*xi;
    
    dN_dxi(1) = -0.5*(1.0 - val);
    dN_dxi(2) = -val;
    dN_dxi(3) =  0.5*(1.0 + val);
    
    d2N_dxi2(1) =  1.0;
    d2N_dxi2(2) = -2.0;
    d2N_dxi2(3) =  1.0;
elseif(p == 3)
    fact1 = 9/16;
    fact2 = 27/16;
    val1 = xi*xi;
    
    N(1) = -fact1 * (1 - xi)   * (1/9 - val1);
    N(2) =  fact2 * (1 - val1) * (1/3 - xi);
    N(3) =  fact2 * (1 - val1) * (1/3 + xi);
    N(4) = -fact1 * (1 + xi)   * ( 1/9 - val1);
    
    val2 = 3.0*val1;

    dN_dxi(1) = -fact1*(-1/9 - 2.0*xi   +  val2);
    dN_dxi(2) =  fact2*(-1   - 2.0/3*xi +  val2);
    dN_dxi(3) =  fact2*(1    - 2.0/3*xi -  val2);
    dN_dxi(4) = -fact1*(1/9  - 2.0*xi   -  val2);
    
    val2 = 6.0*xi;
    
    d2N_dxi2(1) = -fact1 * (-2   + val2);
    d2N_dxi2(2) =  fact2 * (-2/3 + val2);
    d2N_dxi2(3) =  fact2 * (-2/3 - val2);
    d2N_dxi2(4) = -fact1 * (-2   - val2);
elseif(p == 4)
    fact1 = 2/3;
    fact2 = 8/3;
    val1 = xi*xi;
    val2 = val1*xi;
    val3 = val2*xi;
    
    N(1) =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
    N(2) = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
    N(3) =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
    N(4) =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
    N(5) = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);
    
    val4 = 4.0*val2;

    dN_dxi(1) =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
    dN_dxi(2) = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
    dN_dxi(3) =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
    dN_dxi(4) =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
    dN_dxi(5) = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);
    
    val4 = 12.0*val1;
    
    d2N_dxi2(1) =  fact1 * (-0.5  -  6.0*xi  +  val4);
    d2N_dxi2(2) = -fact2 * (-2.0  -  3.0*xi  +  val4);
    d2N_dxi2(3) =    4.0 * (-2.5  -  0.0     +  val4);
    d2N_dxi2(4) =  fact2 * ( 2.0  -  3.0*xi  -  val4);
    d2N_dxi2(5) = -fact1 * ( 0.5  -  6.0*xi  -  val4);
else
    sprintf('no basis functions defined for this degree');
end

Jx=0.0;
Jy=0.0;
x=0.0;
for kk=1:nen
    Jx = Jx + dN_dxi(kk) * XX(NodeNums(kk),1);
    Jy = Jy + dN_dxi(kk) * XX(NodeNums(kk),2);
    x = x + N(kk) * XX(NodeNums(kk));
end
J = sqrt(Jx*Jx+Jy*Jy);

dN_dx   = (1/J) * dN_dxi;
d2N_dx2 = (1/J^2) * d2N_dxi2;
