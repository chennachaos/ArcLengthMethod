
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
more off;
format long;

%fname = "input_LeeFrame-nelem20.txt";
%fname = "input_Arch_semicircle-nelem50-sym.txt";
%fname = "input_Arch_semicircle-nelem50-unsym.txt";
fname = "input_arch-215deg.txt";
%fname = "input_Arch_model1.txt"

%fname = "input_Truss_2D_2members_model1.txt";
%fname = "input_Truss_2D_3members_model1.txt";
%fname = "input_Truss_2D_15deg.txt";
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
dsMax = ds;
dsMin = ds;

loadfactor      = loadincr;
loadfactorPrev2 = 0.0;
loadfactorPrev  = 0.0;

converged = false;
convergedPrev = false;

loadStepConverged = 0;
output = [disp(outputlist)];
llist = [0.0];

dispFull = [disp];

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

        if(ndim == 2)
          if(ndof == 2) % Truss element
            for e = 1:nelem
                [Klocal, Flocal] = Truss_2D_model1(elemData, elemConn, e, coords, disp, bf);

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
         ds = sqrt(Du'*Du + loadfactor*loadfactor*Fext'*Fext);

         dsMax = ds;
         dsMin = ds/1024.0;
      end

      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;
      dispPrev2 = dispPrev;
      dispPrev  = disp;

      dsPrev = ds;
      if(convergedPrev)
        ds = min(max(2.0*ds, dsMin), dsMax);
      endif

      dispFull = [dispFull; disp];
      output = [output disp(outputlist)];
      llist = [llist; loadfactor];
      
%      plot(abs(output(1,:)), llist,'bx-');
%      hold on
      plot(abs(output(2,:)), llist,'bs-');
      hold on

%      linestr = strsplit(fname, ".");
%      figname = strcat(linestr(1){:}, "_", num2str(loadStepConverged), ".pdf")
%      plot_semicircular_arch(coords, disp, figname);
%      plot(abs(output(1,:)), llist, 'bs-'); hold on;
%      plot(abs(output(end,:)), llist,'bx-');
%      plot(abs(dy)/R, (E*I/R/R)*llist.^(-1),'ko-')
%      plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
%      axis([-150 150 -150 150])
%      hold on

%      for e=1:nelem
%        n1 = elemConn(e,3);
%        n2 = elemConn(e,4);
%        xx = [coords(n1,1)+disp(ndof*(n1-1)+1) coords(n2,1)+disp(ndof*(n2-1)+1)];
%        yy = [coords(n1,2)+disp(ndof*(n1-1)+2) coords(n2,2)+disp(ndof*(n2-1)+2)];
%        plot(xx, yy, 'ko-')
%        hold on
%      endfor
%      axis([-0.6 0.6 -3 3])
%      hold off

      loadStepConverged = loadStepConverged + 1;
    else
      if(convergedPrev)
        ds = max(ds*0.5, dsMin);
      else
        ds = max(ds*0.25, dsMin);
      endif
    endif

%    waitforbuttonpress
end

%plot(abs(output(1,:)), llist,'bx-');
%hold on
%plot(abs(output(2,:)), llist,'ko-');
      
%plot(t, dy,'k-')
%figure(1)
%plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
%figure(2)
%plot(dx, llist, 'bs-'); hold on; plot(abs(dy),llist,'ko-')
% axis([0,10,-10,10])

fileID = fopen('solution.dat','w');

for ii=1:size(dispFull)(1)
    fprintf(fileID,'%12.8f \n', dispFull(ii));
end

fclose(fileID)


fileID = fopen('path.dat','w');

for ii=1:size(llist)(1)
    fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \n', llist(ii), output(1,ii), output(2,ii));
%    fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \t %12.8f \n', llist(ii), output(1,ii), output(2,ii), output(3,ii));
end

fclose(fileID)


