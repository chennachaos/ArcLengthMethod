
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
more off;
format long;

fname = "LeeFrame.txt";
%fname = "ArchZienkiewicz.txt";

[nnode, nelem, coords, elemConn, elemData, ndof, nDBC, dbclist, nFBC, fbclist, LM, neq, assy4r, dof_force, Fext] = processfile(fname);

%
soln = zeros(neq,1);

disp = soln;
velo = soln;
acce = soln;

dispPrev2 = soln;
dispPrev3 = soln;
dispPrev4 = soln;

dispPrev = soln;
veloPrev = soln;
accePrev = soln;
dispDotPrev = soln;

dispCur = soln;
veloCur = soln;
acceCur = soln;


Kglobal = zeros(neq,neq);
Fglobal = zeros(neq,1);

dt=1.0;
t=0.0:dt:4.0;

Nt=max(size(t));

td = timeSteppingParameters_Solid(0, 0.0, dt);

%dx = zeros(1,Nt);
%dy = zeros(1,Nt);
vx = zeros(1,Nt);
vy = zeros(1,Nt);
ax = zeros(1,Nt);
ay = zeros(1,Nt);
%llist = dx;

bf=[0.0 0.0];

ds = 0.5;
dsPrev = ds;
dsPrev2 = ds;
ds_max = ds;
ds_min = 0.1;
loadfact = ds;

converged = false;
convergedPrev = false;

loadStep = 2;


for  timeStep=2:50
    printf("time step = %d \n", timeStep);

    if(abs(dispPrev(dof_force(end))) > 90.0)
      break;
    endif

    ds
    dsPrev
    dsFactor1 = ds/dsPrev
    disp = (1.0+dsFactor1)*dispPrev - dsFactor1*dispPrev2;

%    dsFactor1 = (ds^2+ds*(dsPrev+dsPrev2))/dsPrev/dsPrev2;
%    dsFactor2 = -ds*(dsPrev+dsPrev2)/(dsPrev+dsPrev2)/dsPrev2;
%    disp = (1.0+dsFactor1+dsFactor2)*dispPrev - dsFactor1*dispPrev2 - dsFactor2*dispPrev3;

%    Du = (disp(1:neq-1) - dispPrev(1:neq-1)).*sign((1.0+ds/dsPrev)*dispPrev(1:neq-1) - (ds/dsPrev)*dispPrev2(1:neq-1));

    Du = disp(1:neq-1) - dispPrev(1:neq-1);
    Dl = disp(neq)  - dispPrev(neq);
    
    convergedPrev = converged;
    converged = false;

    for iter = 1:10
        Kglobal(1:end,1:end) = 0.0;
        Fglobal(1:end) = 0.0;

        if(timeStep > 2)
          loadfact = disp(neq);
        endif
        loadfact

        for e = 1:nelem
            [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, td, disp, veloCur, acceCur, bf);

            Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
            Fglobal = Assembly_Vector(Fglobal,Flocal,LM,e);
        end

        Fglobal = Fglobal + loadfact*Fext;

%        [converged, du, dl] = solve_arclength(timeStep, neq, iter, Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds);
        [converged, du, dl] = solve_arclength_split(timeStep, neq, iter, Kglobal, Fglobal, dof_force, Fext, assy4r, Du, Dl, ds);

        if(converged)
          break;
        endif

        disp(assy4r(1:end-1)) = disp(assy4r(1:end-1)) + du;
        disp(end) = disp(end) + dl;

        Du(assy4r(1:end-1)) = Du(assy4r(1:end-1)) + du;
        Dl = Dl + dl;
    end

    dsPrev2 = dsPrev;
    dsPrev = ds;
    
    if (converged)
      if(timeStep == 2)
         disp(end)  = ds;
         ds = norm(disp-dispPrev);
         dsPrev = ds;
         dsPrev2 = ds;
         ds_max = ds;
         ds_min = 0.1;
      end

      if(convergedPrev)
        ds = min(max(2.0*ds, ds_min), ds_max);
      endif

      dispPrev4 = dispPrev3;
      dispPrev3 = dispPrev2;
      dispPrev2 = dispPrev;
      dispPrev  = disp;

      dx(loadStep) = disp(dof_force-1);
      dy(loadStep) = disp(dof_force);

%      disp(end)
      llist(loadStep) = disp(end);

%      plot(abs(dx), llist, 'bs-'); hold on;
      plot(abs(dy), llist,'ko-')
%      plot(abs(dy)/R, (E*I/R/R)*llist.^(-1),'ko-')
%      plot(coords(:,1)+disp(1:3:neq-1), coords(:,2)+disp(2:3:neq-1), 'ko-')
%      axis([-150 150 -150 150])
      hold on
      loadStep = loadStep + 1;
    else
%      disp      = dispPrev;
%      dispPrev  = dispPrev2;
%      dispPrev2 = dispPrev3;
%      dispPrev3 = dispPrev4;

      if(convergedPrev)
        ds = max(ds*0.5, ds_min);
      else
%        beta = 10.0;
        ds = max(ds*0.25, ds_min);
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





