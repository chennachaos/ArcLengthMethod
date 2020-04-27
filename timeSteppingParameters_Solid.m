function td = timeSteppingParameters_Solid(tis, rho, dt)

td = zeros(100,1);

switch tis
   case 0
     alpf = 1.0;
     alpm = 1.0;
     gamm = 1.0;

     td(1)  = alpm;
     td(2)  = alpf;
     td(3)  = alpm;
     td(4)  = gamm;
     td(7)  = alpf;

            % velocity is used as the primary variable for 
            % solid dynamics problem 
            % td[10] is the multiplication when converting from
            % displacement based formulation to velocity based formulation
            % It is set to ONE for static problem

      td(10) = 1.0; 

    case 2 % Backward-Euler

       alpf = 1.0;
       alpm = 1.0;
       gamm = 1.0;
       
       td(1)  = alpm;
       td(2)  = alpf;
       td(3)  = alpm;
       td(4)  = gamm;

       td(5)  = 1.0/dt/dt;
       td(6)  = 1.0/dt;
       td(7)  = alpf;
       
       %displacement as the primary variable
       %v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
       %a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;

       td(10) = 1.0/dt;    % d_{n+1}
       td(11) = -td(10);   % d_n
       td(12) = 0.0;       % v_n
       td(13) = 0.0;       % a_n
       td(14) = 0.0;       % ddot_n
       
       td(15) = 1.0/dt/dt; % d_{n+1}
       td(16) = -td(15);   % d_n
       td(17) = -1.0/dt;   % vn
       td(18) = 0.0;       % an
       td(19) = 0.0;       % ddot_n

       %velocity as the primary variable
       %d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
       %a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;

       td(40) = dt;        % v_{n+1}
       td(41) = 1.0;       % d_n
       td(42) = 0.0;       % v_n
       td(43) = 0.0;       % a_n
       td(44) = 0.0;       % ddot_n
       
       td(45) = 1.0/dt;    % v_{n+1}
       td(46) = 0.0;       % d_n
       td(47) = -td(45);   % v_n
       td(48) = 0.0;       % a_n
       td(49) = 0.0;       % ddot_n

   case 3 % CH-aplha

       alpm = (2.0-rho)/(rho+1.0);
       alpf = 1.0/(rho+1.0);
       
       gamm = 0.5 + alpm - alpf;
       beta = 0.25*(1.0+alpm-alpf)*(1.0+alpm-alpf);
       
       td(1)  = alpm;
       td(2)  = alpf;
       td(3)  = alpm;
       td(4)  = gamm;
       
       td(5)  = alpm/beta/dt/dt;
       td(6)  = alpf*gamm/beta/dt;
       td(7)  = alpf;

       %displacement as the primary variable
       %v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
       %a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;

       td(10) = gamm/beta/dt;           % d_{n+1}
       td(11) = -td(10);                % d_n
       td(12) = 1.0-gamm/beta;          % v_n
       td(13) = dt*(1.0-gamm/2.0/beta); % a_n
       td(14) = 0.0;                    % ddot_n
       
       td(15) = 1.0/beta/dt/dt;         % d_{n+1}
       td(16) = -td(15);                % d_n
       td(17) = -1.0/beta/dt;           % v_n
       td(18) = 1.0-1.0/2.0/beta;       % a_n
       td(19) = 0.0;                    % ddot_n

       %velocity as the primary variable
       %d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
       %a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;

       td(40) = dt*beta/gamm;                     % v_{n+1}
       td(41) = 1.0;                              % d_n
       td(42) = dt*(gamm-beta)/gamm;              % v_n
       td(43) = dt*dt*(gamm-2.0*beta)/(2.0*gamm); % a_n
       td(44) = 0.0;                              % ddot_n

       td(45) = 1.0/(gamm*dt);                    % v_{n+1}
       td(46) = 0.0;                              % d_n
       td(47) = -td(45);                          % v_n
       td(48) = (gamm-1.0)/gamm;                  % a_n
       td(49) = 0.0;                              % ddot_n

   case 4 % JHW-alpha or KDP-alpha

       alpf = 1.0/(1.0 + rho);
       alpm = 0.5*(3.0 - rho)/(1.0 + rho);

       gamm = 0.5 + alpm - alpf;

       td(1)  = alpm;
       td(2)  = alpf;
       td(3)  = alpm;
       td(4)  = gamm;

       td(5)  = (alpm*alpm)/(alpf*gamm*gamm*dt*dt);
       td(6)  = alpm/gamm/dt;
       td(7)  = alpf;
       
       %displacement as the primary variable
       %v_{n+1}    = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
       %a_{n+1}    = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;
       %ddot_{n+1} = td(20)*d_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n;

       td(10) = alpm/(alpf*gamm*dt);             % d_{n+1}
       td(11) = -td(10);                         % d_n
       td(12) = (alpf-1.0)/alpf;                 % v_n
       td(13) = 0.0;                             % a_n
       td(14) = (gamm-alpm)/alpf/gamm;           % ddot_n
       
       td(15) = alpm/(alpf*gamm*gamm*dt*dt);     % d_{n+1}
       td(16) = -td(15);                         % d_n
       td(17) = -1.0/(alpf*gamm*dt);             % v_n
       td(18) = (gamm-1.0)/gamm;                 % a_n
       td(19) = (gamm-alpm)/(alpf*gamm*gamm*dt); % ddot_n

       td(20) = 1.0/(gamm*dt);                   % d_{n+1}
       td(21) = -td(20);                         % d_n
       td(22) = 0.0;                             % v_n
       td(23) = 0.0;                             % a_n
       td(24) = (gamm-1.0)/gamm;                 % ddot_n
       
       %velocity as the primary variable
       %d_{n+1}    = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
       %a_{n+1}    = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;
       %ddot_{n+1} = td(30)*v_{n+1} + td(31)*d_n + td(32)*v_n + td(33)*a_n + td(34)*ddot_n;

       td(40) = alpf*gamm*dt/alpm;           % v_{n+1}
       td(41) = 1.0;                         % d_n
       td(42) = (1.0-alpf)*gamm*dt/alpm;     % v_n
       td(43) = 0.0;                         % a_n
       td(44) = (alpm-gamm)*dt/alpm;        % ddot_n
       
       td(45) = 1.0/gamm/dt;                 % v_{n+1}
       td(46) = 0.0;                         % d_n
       td(47) = -td(45);                     % v_n
       td(48) = 1.0-1.0/gamm;                % v_n
       td(49) = 0.0;                         % ddot_n
       
       td(50) = alpf/alpm;                   % v_{n+1}
       td(51) = 0.0;                         % d_n
       td(52) = (1.0-alpf)/alpm;             % v_n
       td(53) = 0.0;                         % a_n
       td(54) = (alpm-1.0)/alpm;             % ddot_n

   case 5 %% Newmark-beta method

       alpm = 1.0;
       alpf = 1.0;

       td(1)  = alpm;
       td(2)  = alpf;
       td(3)  = alpm;
       td(4)  = gamm;
       
       td(5)  = alpm/beta/dt/dt;
       td(6)  = alpf*gamm/beta/dt;
       td(7)  = alpf;

       %displacement as the primary variable
       %v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
       %a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;

       td(10) = gamm/beta/dt;           % d_{n+1}
       td(11) = -td(10);                % d_n
       td(12) = 1.0-gamm/beta;          % v_n
       td(13) = dt*(1.0-gamm/2.0/beta); % a_n
       td(14) = 0.0;                    % ddot_n

       td(15) = 1.0/beta/dt/dt;         % d_{n+1}
       td(16) = -td(15);                % d_n
       td(17) = -1.0/beta/dt;           % v_n
       td(18) = 1.0-1.0/2.0/beta;       % a_n
       td(19) = 0.0;                    % ddot_n

       %velocity as the primary variable
       %d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
       %a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;

       td(40) = dt*beta/gamm;                     % v_{n+1}
       td(41) = 1.0;                              % d_n
       td(42) = dt*(gamm-beta)/gamm;              % v_n
       td(43) = dt*dt*(gamm-2.0*beta)/(2.0*gamm); % a_n
       td(44) = 0.0;                              % ddot_n

       td(45) = 1.0/(gamm*dt);                    % v_{n+1}
       td(46) = 0.0;                              % d_n
       td(47) = -td(45);                          % v_n
       td(48) = (gamm-1.0)/gamm;                  % a_n
       td(49) = 0.0;                              % ddot_n
    
    otherwise
       printf('tis error\n');
       printf('This time integration scheme is not available yet!');
end


