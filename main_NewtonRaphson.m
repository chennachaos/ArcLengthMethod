
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
more off;
format long g;

%fname = "input_Truss_2D_3members_model1.txt";
%fname = "input_Truss_3D_2members.txt";
%fname = "input_Truss_3D_12members.txt";

%fname = "input_LeeFrame-nelem10.txt";
%fname = "input_LeeFrame-nelem20.txt";
%fname = "input_arch-215deg.txt";
%fname = "input_Arch_semicircle-nelem50-sym.txt";
%fname = "input_Arch_semicircle-nelem50-unsym.txt";
fname = "input-beamEndMoment-nelem20.txt";

%fname = "lattice2.txt"

[ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, dofs_free, dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)

disp = zeros(neq,1);

dispPrev  = disp;
dispPrev2 = disp;
dispPrev3 = disp;
dispPrev4 = disp;


Kglobal = zeros(neq,neq);
Rglobal = zeros(neq,1);


bf=[0.0 0.0];

loadincr = 1.0/double(maxloadSteps);

loadfactor      = 0.0;
loadfactorPrev  = 0.0;
loadfactorPrev2 = 0.0;

converged = false;
convergedPrev = false;

loadStepConverged = 0;
output = [disp(outputlist)];

dispFull = [disp];

for  loadStep=1:maxloadSteps
    loadfactor = loadfactor + loadincr;

    fprintf("load step = %d \t load factor = %f \n", loadStep, loadfactor);

    DsFactor1 = 1.0;
    disp     = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2;

    convergedPrev = converged;
    converged = false;

    for iter = 1:10
        Kglobal(1:end,1:end) = 0.0;
        Rglobal(1:end) = 0.0;

        if(ndim == 2)
          if(ndof == 2) % Truss element
            for e = 1:nelem
                [Klocal, Flocal] = Truss_2D_model1(elemData, elemConn, e, coords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          else % Beam element
            for e = 1:nelem
                [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          end
        else
          if(ndof == 3) % Truss element
            for e = 1:nelem
                [Klocal, Flocal] = Truss_3D_model2(elemData, elemConn, e, coords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          end
        end

        Rglobal = Rglobal + loadfactor*Fext;

        rNorm = norm(Rglobal(dofs_free),2);

        fprintf(' rNorm : %5d ...  %12.6E \n', iter, rNorm);

        if(rNorm < 1.0e-8)
          converged = true;
          break;
        end

        du = Kglobal(dofs_free,dofs_free)\Rglobal(dofs_free);
        disp(dofs_free) = disp(dofs_free) + du;
    end

    if (converged)
%      disp
      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;
      dispPrev2 = dispPrev;
      dispPrev  = disp;

      dispFull = [dispFull; disp];
      output = [output disp(outputlist)];

     figure(1);
     for e=1:nelem
       n1 = elemConn(e,3);
       n2 = elemConn(e,4);
       xx = [coords(n1,1)+disp(ndof*(n1-1)+1) coords(n2,1)+disp(ndof*(n2-1)+1)];
       yy = [coords(n1,2)+disp(ndof*(n1-1)+2) coords(n2,2)+disp(ndof*(n2-1)+2)];
       plot(xx, yy, 'ko-')
       hold on
     end
     %axis([-5.0 12.0 -2 6])
     %axis([-20. 150.0 -20 150])
     hold off

      loadStepConverged = loadStepConverged + 1;
    else
      loadfactor = loadfactorPrev;
      loadincr = loadincr*0.5;
    end

%    waitforbuttonpress
end


