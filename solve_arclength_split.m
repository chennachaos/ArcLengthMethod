function  [converged, du, dl, du1] = solve_arclength_split(timeStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, ds, du1)

    psi = 1.0;

    FextReduced = Fext(assy4r);
    FtF = Fext'*Fext;
    if(timeStep > 1)
        A = Du'*Du + psi*Dl*Dl*FtF - ds*ds;

        a = 2.0*Du(assy4r)';
        b = 2.0*psi*Dl*FtF;
    else
        A = 0.0;
        a = 0.0*Du(assy4r)';
        b = 1.0;
    end

    %%% Applying Boundary Conditions

    R = Rglobal(assy4r);

    rNorm = norm(R,2);
    rNorm = sqrt(rNorm*rNorm + A*A);

    fprintf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
    du = R*0.0;
    dl = 0.0;
    converged = false;

    if(rNorm < 1.0e-6)
       converged = true;
       return;
    end

    K1 = Kglobal(assy4r,assy4r);
    [L, U, P] = lu(sparse(K1));

    %% solve the matrix system
    duu = L\(P*FextReduced);
    du1 = U\duu;

    duu = L\(P*R);
    du2 = U\duu;
    du2 = -du2; % this is because the Residual is added to the RHS

    dl = (a*du2 - A)/(b+a*du1);

    du = -du2 + dl*du1;
end