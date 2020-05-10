function  [converged, du, dl] = solve_arclength(timeStep, neq, iter, Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds)

    if(timeStep > 1)
        dummy = Du'*Du + Dl*Dl - ds*ds;
        Fglobal(neq) = Fglobal(neq) - dummy;

        Kglobal(:, neq)   = Kglobal(:, neq) - Fext;

        Kglobal(neq, 1:neq-1) = Kglobal(neq, 1:neq-1) + 2.0*Du';

        Kglobal(neq, neq)  = Kglobal(neq, neq)  + 2.0*Dl;
    else
        Fglobal(neq) = 0.0 ;
        Kglobal(neq, neq) = 1.0;
    endif

    %%% Applying Boundary Conditions
        
    K1 = Kglobal(assy4r,assy4r);
    F1 = Fglobal(assy4r);

    rNorm = norm(F1,2);

    printf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);
    du = F1*0.0;
    dl = 0.0;
    converged = false;

    if(rNorm < 1.0e-8)
       converged = true;
       return;
    end
        
    %% solve the matrix system
    dutemp = K1\F1;

    du = dutemp(1:end-1);
    dl = dutemp(end);
endfunction