function  [converged, du, dl] = solve_arclength_split(timeStep, neq, iter, Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds)

    if(timeStep > 2)
        A = Du'*Du + Dl*Dl - ds*ds;

        a = 2.0*Du(assy4r(1:end-1))';
        b = 2.0*Dl;
    else
        A = 0.0 ;
        a = 0.0*Du(assy4r(1:end-1))';
        b = 1.0;
    endif

    %%% Applying Boundary Conditions
        
    F1 = Fglobal(assy4r(1:end-1));

    rNorm = norm(F1,2);
    rNorm = sqrt(rNorm*rNorm + A*A);

    printf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
    du = F1*0.0;
    dl = 0.0;
    converged = false;

    if(rNorm < 1.0e-6)
       converged = true;
       return;
    end
        
    K1 = Kglobal(assy4r(1:end-1),assy4r(1:end-1));
    [L, U, P] = lu(sparse(K1));

    %% solve the matrix system
    duu = L\(P*Fext(assy4r(1:end-1)));
    du1 = U\duu;
 
    duu = L\(P*F1);
    du2 = U\duu;

    du2 = -du2; % this is because the Residual is added to the RHS

    dl = (a*du2 - A)/(b+a*du1);

    du = -du2 + dl*du1;
endfunction