
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
more off;
format long;

%fname = "input_Truss_2D_3members_model1.txt";
%fname = "input_Truss_3D_2members.txt";
fname = "input_Truss_3D_12members.txt";

%fname = "input_LeeFrame-nelem20.txt";
%fname = "input_arch-215deg.txt";
%fname = "input_Arch_semicircle-nelem50-sym.txt";
%fname = "input_Arch_semicircle-nelem50-unsym.txt";
%fname = "input-beamEndMoment-nelem10.txt";

[ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist] = processfile(fname)

disp = zeros(neq,1);

dispPrev  = disp;
dispPrev2 = disp;
dispPrev3 = disp;
dispPrev4 = disp;


Kglobal = zeros(neq,neq);
Rglobal = zeros(neq,1);


bf=[0.0 0.0];

Ds = loadincr;
DsPrev = Ds;
DsMax = Ds;
DsMin = Ds;

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
    fprintf("load step = %d \n", loadStep);

    if(loadStep > 1)
      Ds
      DsPrev
      DsFactor1 = Ds/DsPrev
      disp     = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2;
      loadfactor = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2;
    end

    Du = disp - dispPrev;
    Dl = loadfactor - loadfactorPrev;

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

%        [converged, du, dl] = solve_arclength(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);
        [converged, du, dl] = solve_arclength_split(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);

        if(converged)
          break;
        end

        disp(assy4r) = disp(assy4r) + du;
        loadfactor = loadfactor + dl;

        Du(assy4r) = Du(assy4r) + du;
        Dl = Dl + dl;
    end

    if (converged)
%      disp
      if(loadStep == 1)
         Ds = sqrt(Du'*Du + loadfactor*loadfactor*Fext'*Fext);

         DsMax = Ds;
         DsMin = Ds/1024.0;
      end

      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;
      dispPrev2 = dispPrev;
      dispPrev  = disp;

      DsPrev = Ds;
      if(convergedPrev)
        Ds = min(max(2.0*Ds, DsMin), DsMax);
      end

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
        Ds = max(Ds*0.5, DsMin);
      else
        Ds = max(Ds*0.25, DsMin);
      end
    end

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

for ii=1:size(dispFull,1)
    fprintf(fileID,'%12.8f \n', dispFull(ii));
end

fclose(fileID)


fileID = fopen('path.dat','w');

for ii=1:size(llist,1)
    fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \n', llist(ii), output(1,ii), output(2,ii));
%    fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \t %12.8f \n', llist(ii), output(1,ii), output(2,ii), output(3,ii));
end

fclose(fileID)


