
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
more off;
format long;

fname = "input_LeeFrame.txt";
%fname = "input_Arch_semicircle.txt";
%fname = "input_ArchZienkiewicz.txt";
%fname = "input_Arch_model1.txt"

%fname = "input_Truss_2D_2members.txt";
%fname = "input_Truss_3D_2members.txt";
%fname = "input_Truss_3D_12members.txt";

[ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)

disp = zeros(neq,1);

dispPrev  = disp;
dispPrev2 = disp;
dispPrev3 = disp;
dispPrev4 = disp;


Kglobal = zeros(neq,neq);
Fglobal = zeros(neq,1);


bf=[0.0 0.0];

ds = loadincr;
dsPrev = ds;
dsPrev2 = ds;
ds_max = ds;
ds_min = ds;

loadfactor      = loadincr;
loadfactorPrev2 = 0.0;
loadfactorPrev  = 0.0;

converged = false;
convergedPrev = false;

loadStepConverged = 2;
output = [disp(outputlist)];
llist = [0.0];


for  loadStep=1:maxloadSteps
    printf("load step = %d \n", loadStep);

    if(loadStep > 1)
      ds
      dsPrev
      dsFactor1 = ds/dsPrev
      disp     = (1.0+dsFactor1)*dispPrev - dsFactor1*dispPrev2;
      loadfactor = (1.0+dsFactor1)*loadfactorPrev - dsFactor1*loadfactorPrev2;
    endif

    Du = disp - dispPrev;
    Dl = loadfactor - loadfactorPrev;

    convergedPrev = converged;
    converged = false;

    for iter = 1:10
        Kglobal(1:end,1:end) = 0.0;
        Fglobal(1:end) = 0.0;

%        loadfactor
        
        if(ndim == 2)
          if(ndof == 2) % Truss element
            for e = 1:nelem
                [Klocal, Flocal] = Truss_2D_model2(elemData, elemConn, e, coords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Fglobal = Assembly_Vector(Fglobal,Flocal,LM,e);
            end
          else % Beam element
            for e = 1:nelem
                [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Fglobal = Assembly_Vector(Fglobal,Flocal,LM,e);
            end
          endif
        else
          if(ndof == 3) % Truss element
            for e = 1:nelem
                [Klocal, Flocal] = Truss_3D_model2(elemData, elemConn, e, coords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Fglobal = Assembly_Vector(Fglobal,Flocal,LM,e);
            end
          endif
        endif

        Fglobal = Fglobal + loadfactor*Fext;

%        [converged, du, dl] = solve_arclength(loadStep, neq, iter, Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds);
        [converged, du, dl] = solve_arclength_split(loadStep, neq, iter, Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds);

        if(converged)
          break;
        endif

        disp(assy4r) = disp(assy4r) + du;
        loadfactor = loadfactor + dl;

        Du(assy4r) = Du(assy4r) + du;
        Dl = Dl + dl;
    end

    if (converged)
      if(loadStep == 1)
         ds = norm(disp-dispPrev);
         ds = sqrt(ds*ds + loadfactor*loadfactor)

         dsPrev = ds;
         dsPrev2 = ds;
         ds_max = ds;
         ds_min = ds*1.0e-3;
      end

      dsPrev2 = dsPrev;
      dsPrev = ds;

      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;

      if(convergedPrev)
        ds = min(max(2.0*ds, ds_min), ds_max);
      endif

      dispPrev4 = dispPrev3;
      dispPrev3 = dispPrev2;
      dispPrev2 = dispPrev;
      dispPrev  = disp;

      output = [output disp(outputlist)];
      llist = [llist; loadfactor];

%      plot(abs(output(1,:)), llist, 'bs-'); hold on;
      plot(abs(output(end,:)), llist,'ko-');
%      plot(abs(dy)/R, (E*I/R/R)*llist.^(-1),'ko-')
%      plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
%      axis([-150 150 -150 150])
      hold on
      loadStepConverged = loadStepConverged + 1;
    else
      if(convergedPrev)
        ds = max(ds*0.5, ds_min);
      else
        ds = max(ds*0.1, ds_min);
      endif
    endif

%    waitforbuttonpress
end

%plot(t, dy,'k-')
%figure(1)
%plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
%figure(2)
%plot(dx, llist, 'bs-'); hold on; plot(abs(dy),llist,'ko-')
% axis([0,10,-10,10])

%fileID = fopen('solution.dat','w');
%
%for ii=1:Nt
%    fprintf(fileID,'%12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n', t(ii), dx(ii), dy(ii), vx(ii), vy(ii), ax(ii), ay(ii));
%end





