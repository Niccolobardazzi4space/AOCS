%% S/C 12U : detumbling, slew manoeuvre & Earth pointing
% 
%  CREATED BY:             Person Code:           Matriculation Number:
%    Niccolò Bardazzi          10800456                         963039
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  ATTITUDE DYNAMICS PROJECT  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear  
clc
load('sc_data.mat')

%% Orbital parameters
mu = 398600;
a = 7197;
e = 2.0841148e-04;
i = 98.7116;
aop = 95.683;
RAAN = 271.7368;
T = 2*pi/sqrt(mu/a^3);
Greenwich_lon = 180;

%% Components 

I = [11.05,18.25,18.72]*1.e-2;       % inertia [kg*m^2]

m_dipole = [0.01,0.05,0.01];  % due to electronic devices

l = [  0  0.2 0.2    % surfaces' dimensions
      0.3  0  0.2
      0.3 0.2  0
       0  0.2 0.2  
      0.3  0  0.2
      0.3 0.2  0
      0.3 0.2  0     % panels' dimensions (2 panels with 2 normals)
      0.3 0.2  0
      0.3 0.2  0
      0.3 0.2  0 ]; 
  
n = [ 1 0 0          % normal to each surface with respect to body frame
      0 1 0
      0 0 1
     -1 0 0
     0 -1 0
     0 0 -1
      0 0 1
     0 0 -1
      0 0 1
     0 0 -1 ];


% Variable thrusters parameters:
Fmax = 10e-3;             % Thrust of a thruster   [N]
Fmin = 100e-6;            

% trapezoidal approximation on the thrusters' torque 
t_raise  = 10.e-3;          % rise time              [s]
t_fall  = 50.e-3;           % fall time              [s]
t_delay = 5.e-3;            % delay time             [s]
alfaxy = 30;      % [°] inclination of thrusters in xy plane 
p = l(1,3);
x = 0.75*p;       % [m] distance from center of thrusters

conf = [ 1  1 -1 -1
        -1  1  1 -1
         1 -1  1 -1 ];    % piramid configuration (square base)

dist = [ p*sind(alfaxy)       0                  0
              0        p*cosd(alfaxy)             0
              0                0       x*sind(alfaxy)-x*cosd(alfaxy) ];
          
R_mat = dist*conf;
nR = null(R_mat);

Tmax_x = 2*p*Fmax*cosd(alfaxy);
Tmax_y = 2*p*Fmax*sind(alfaxy);
Tmax_z = 2*(x*Fmax*cosd(alfaxy)-x*Fmax*sind(alfaxy));
Tmax = diag([Tmax_x;Tmax_y;Tmax_z]);

Tmin_x = 2*p*Fmin*cosd(alfaxy);
Tmin_y = 2*p*Fmin*sind(alfaxy);
Tmin_z = 2*(x*Fmin*cosd(alfaxy)-x*Fmin*sind(alfaxy));
Tmin = diag([Tmin_x;Tmin_y;Tmin_z]);

% Analog sun sensor
f_ss = 5;               % Update frequency
ss_accuracy  = 0.100;   % [°] 

% EH sensor
f_eh = 2;                % Update frequency [Hz]
eh_accuracy = 0.250;      % [°]  

% Solar panles
roMBd = 0.1;
roMBs = 0.5;
roSPd = 0.1;
roSPs = 0.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DETUMBLING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETUMBLING for determination:

disp('Detumbling has begun')

t_det_det = 2; % [s]
% Detumbling initial conditions:
w0=[0.22 0.26 0.22];
Abn0=[0.5335 0.808 0.25
     -0.808 0.3995 0.433
      0.25 -0.433 0.866];
  
% Target
w_target = [0,0,0];

%% Run simulation
sim('detumb_det') % gain the attitude only (no control applied yet)

%% DETUMBLING WITH EKF
Abn0 = ans.Abn_determination.Data(:,:,end); 
w0 = ans.omega.Data(end,:);

%% EKF
Q = diag([100 100 100 100 0.1 0.1 0.1]);
P0 = diag([1e3 1e3 1e3 1e3 100 100 100])/10;
q0(4) = 0.5*sqrt(1+Abn0(1,1)+Abn0(2,2)+Abn0(3,3));
q0(1) = 1/(4*q0(4))*(Abn0(2,3)-Abn0(3,2));
q0(2) = 1/(4*q0(4))*(Abn0(3,1)-Abn0(1,3));
q0(3) = 1/(4*q0(4))*(Abn0(1,2)-Abn0(2,1));
decimation = 0.1; % [s] 
R = [(ss_accuracy*pi/180)^2
     (eh_accuracy*pi/180)^2];

%% Run simulation

t_det_ekf = 150;
t_before = t_det_det;
sim('detumb_ekf')

%%
% Data
wx_filt = ans.omega_filt.Data(:,1); 
wy_filt = ans.omega_filt.Data(:,2);
wz_filt = ans.omega_filt.Data(:,3);
wx = ans.omega.Data(:,1); 
wy = ans.omega.Data(:,2);
wz = ans.omega.Data(:,3);
Meffx_t = ans.Meff_t.Data(:,1);
Meffy_t = ans.Meff_t.Data(:,2);
Meffz_t = ans.Meff_t.Data(:,3);
F1 = ans.F.Data(:,1);
F2 = ans.F.Data(:,2);
F3 = ans.F.Data(:,3);
F4 = ans.F.Data(:,4);
t = linspace(0,t_det_ekf,length(wx_filt));
Abn_end_detumb = ans.Abn.Data(:,:,end);
w_end_detumb = [wx(end) wy(end) wz(end)];

%% Detumbling Plots

% angular velocities
figure(1)
settings

subplot(2,1,1);
plot(t,wx,t,wy,t,wz)
title('Kinematics');
grid minor
axis tight
ylabel('$\omega$ [rad/s]'), xlabel('$t$ [s]')

subplot(2,1,2);
plot(t,wx_filt,t,wy_filt,t,wz_filt)
title('EKF');
grid minor
axis tight
ylabel('$\omega$ [rad/s]'), xlabel('$t$ [s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$','Interpreter','latex');
Lgnd = legend('show');
Lgnd.Position(1) = 0.87;
Lgnd.Position(2) = 0.80;

% Actuators
figure(2)
settings

plot(t,Meffx_t,t,Meffy_t,t,Meffz_t)
title('Thrusters');
grid minor
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')
legend('$M_{effx}$','$M_{effy}$','$M_{effz}$');
Lgnd = legend('show');
Lgnd.Position(1) = 0.83;
Lgnd.Position(2) = 0.83;
 
figure(3)
subplot(2,2,1)
plot(t,F1,'m')
grid minor 
legend('$F_1$');
ylabel('$F_1$ [N]'), xlabel('$t$ [s]')
subplot(2,2,2)
plot(t,F2,'b')
sgtitle('Force');
grid minor 
legend('$F_2$');
ylabel('$F_2$ [N]'), xlabel('$t$ [s]')
subplot(2,2,3)
plot(t,F3,'c')
sgtitle('Force');
grid minor 
legend('$F_3$');
ylabel('$F_3$ [N]'), xlabel('$t$ [s]')
subplot(2,2,4)
plot(t,F4,'g')
sgtitle('Force');
grid minor 
ylabel('$F_4$ [N]'), xlabel('$t$ [s]')
sgtitle('Forces');
legend('$F_4$');

disp('Detumbling is over') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SLEW MANOEUVRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLEW
t_slew = 200;
t_before = t_det_det+t_det_ekf;

w0 = w_end_detumb;
Abn0 = Abn_end_detumb;

% Reaction wheel 
A = [ 1 0 0 1/sqrt(3)
      0 1 0 1/sqrt(3)
      0 0 1 1/sqrt(3)];    % pyramid configuration
As = pinv(A);
  
Mmax = 2.3e-3;
hmax = 10e-3;

% target conditions:
w_target = [0,0,2*pi/T];

% EKF
Q = diag([100 100 100 100 0.1 0.1 0.1]);
P0 = diag([1e3 1e3 1e3 1e3 100 100 100])/1000;
q0(4) = 0.5*sqrt(1+Abn0(1,1)+Abn0(2,2)+Abn0(3,3));
q0(1) = 1/(4*q0(4))*(Abn0(2,3)-Abn0(3,2));
q0(2) = 1/(4*q0(4))*(Abn0(3,1)-Abn0(1,3));
q0(3) = 1/(4*q0(4))*(Abn0(1,2)-Abn0(2,1));
decimation = 0.1; % [s] 
R = [(ss_accuracy*pi/180)^2
     (eh_accuracy*pi/180)^2];


%% Run simulation
clear ans
sim('slew_ekf');

%%
% Data
wx_filt = ans.omega_filt.Data(:,1); 
wy_filt = ans.omega_filt.Data(:,2);
wz_filt = ans.omega_filt.Data(:,3);
wx = ans.omega.Data(:,1); 
wy = ans.omega.Data(:,2);
wz = ans.omega.Data(:,3);
Meffx_r = ans.Meff_r.Data(:,1);
Meffy_r = ans.Meff_r.Data(:,2);
Meffz_r = ans.Meff_r.Data(:,3);
SRP_x = squeeze(ans.SRP.Data(1,1,:));
SRP_y = squeeze(ans.SRP.Data(1,2,:));
SRP_z = squeeze(ans.SRP.Data(1,3,:));
SRP = [SRP_x, SRP_y, SRP_z];
drag_x = squeeze(ans.drag.Data(1,1,:));
drag_y = squeeze(ans.drag.Data(1,2,:));
drag_z = squeeze(ans.drag.Data(1,3,:));
drag = [drag_x, drag_y, drag_z];
gg_x = ans.gg.Data(:,1);
gg_y = ans.gg.Data(:,2);
gg_z = ans.gg.Data(:,3);
gg = [gg_x, gg_y, gg_z];
mag_x = squeeze(ans.mag.Data(1,1,:));
mag_y = squeeze(ans.mag.Data(1,2,:));
mag_z = squeeze(ans.mag.Data(1,3,:));
h_tot = ans.h.Data;
mag = [mag_x, mag_y, mag_z];
t = linspace(0,t_slew,length(wx));
Abn_end_slew = ans.Abn.Data(:,:,end);
w_end_slew = [wx(end) wy(end) wz(end)];
Ae = ans.Ae.Data;

%% Slew plots

% angular velocities
figure(4)
settings

subplot(1,2,1);
plot(t,wx,t,wy,t,wz)
title('Kinematics');
grid minor
axis tight
ylabel('$\omega$ [rad/s]'), xlabel('$t$ [s]')

subplot(1,2,2);
plot(t,wx_filt,t,wy_filt,t,wz_filt)
title('EKF');
grid minor
axis tight
ylabel('$\omega$ [rad/s]'), xlabel('$t$ [s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$','Interpreter','latex');
Lgnd = legend('show');
Lgnd.Position(1) = 0.84;
Lgnd.Position(2) = 0.84;

% Reaction wheels
figure(5)
settings
plot(t,Meffx_r,t,Meffy_r,t,Meffz_r)
title('Torque ');
grid minor
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')
legend('$M_{effx}$','$M_{effy}$','$M_{effz}$');
Lgnd = legend('show');
Lgnd.Position(1) = 0.83;
Lgnd.Position(2) = 0.83;

% Total momentum exchange 
figure(6)
settings
plot(t,h_tot(:,1),t,h_tot(:,2),t,h_tot(:,3),t,h_tot(:,4))
title('Total momentum exchange');
grid minor
axis tight
ylabel('$h$ [mNm$\cdot$s]'), xlabel('$t$ [s]')
legend('$h_1$','$h_2$','$h_3$','$h_4$');
Lgnd = legend('show');
Lgnd.Position(1) = 0.83;
Lgnd.Position(2) = 0.83;


% Perturbations
figure(7)
settings

subplot(2,2,1)
plot(t,SRP);
title('SRP')
grid minor
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')

subplot(2,2,2)
plot(t,drag)
title('Drag')
grid minor 
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')

subplot(2,2,3)
plot(t,mag)
title('Mag. field');
grid minor 
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')
legend('$M_x$','$M_y$','$M_z$');

subplot(2,2,4)
plot(t,gg)
title('G. gradient');
grid minor 
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')
legend('$M_x$','$M_y$','$M_z$');
Lgnd = legend('show');
Lgnd.Position(1) = 0.86;
Lgnd.Position(2) = 0.83;

% Error matrix
figure(8)
Ae1 = squeeze(Ae(1,:,:));
Ae2 = squeeze(Ae(2,:,:));
Ae3 = squeeze(Ae(3,:,:));
plot(t,Ae1,t,Ae2,t,Ae3);
title('Error on attitude matrix');
grid minor 
axis tight
xlabel('$t$ [s]')

disp('Slew manoeuvre is over')

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  EARTH POINTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EARTH POINTING
disp('Earth pointing has begun')

t_ep = T;
w0 = w_end_slew;
Abn0 = Abn_end_slew;

% target conditions:
w_target = [0,0,2*pi/T];

t_before = t_det_det+t_det_ekf+t_slew;

% %% LQR max value Earth pointing
% % pitch axis 
% Ap = [0 1.e-3 
%       1 0];
% Bp = [1/I(3)
%         0    ];
% Cp = [1 0
%       0 1];
% 
% % roll and yaw axes
% n_E = 2*pi/T;
% Kx = (I(3)-I(2))/I(1);
% Ky = (I(3)-I(1))/I(2);
% Kz = (I(2)-I(1))/I(3);
% Axy = [    0   -n_E*(Kx-1) -n_E^2*Kx    0
%     -n_E*(1-Ky)   0        0   -n_E^2*Ky
%          1        0        0      0
%          0        1        0      0    ];
% Bx = [1/I(1)
%         0
%         0
%         0    ];
% Cx = [1 0 0 0
%       0 0 1 0];
% By = [  0
%       1/I(2)
%         0
%         0    ];
% Cy = [0 1 0 0
%       0 0 0 1];
% 
% % LQR weights
% % performance on pitch and yaw,roll
% tmax = [0.0001 0.0001 0.00008 0.0001 0.0001 0.00007]; % dy,dr,dp,y,r,p
% z1 = 1/tmax(1)^2;
% z2 = 1/tmax(2)^2;
% z3 = 1/tmax(3)^2;
% z4 = 1/tmax(4)^2;
% z5 = 1/tmax(5)^2;
% z6 = 1/tmax(6)^2;
% Wzz_p = diag([z3 z5]);
% Wzz = diag([z1 z2 z4 z5]);
% % input
% umax = [2 2 2]*1.e-4;
% u1 = 1/umax(1)^2;
% u2 = 1/umax(2)^2;
% u3 = 1/umax(3)^2;
% Wuu_p = u3;
% Wuu = u1;  % it's the same for yaw abd roll so I picked the one on yaw
% % LQR solutions on pitch 
% Cz_p = [1 0
%         0 1];
% Q = Cz_p'*Wzz_p*Cz_p;
% R = Wuu_p;
% K_p = lqr(Ap,Bp,Q,R,[]);
% % LQR solutions on yaw,roll 
% Cz = [1 0 0 0
%       0 1 0 0
%       0 0 1 0
%       0 0 0 1];
% Q = Cz'*Wzz*Cz;
% R = Wuu;
% K_yaw = lqr(Axy,Bx,Q,R,[]);
% K_roll = lqr(Axy,By,Q,R,[]);
% K = K_roll;

K = [ 0.1 0.1 0.01 0.01];
K_p = [1.1 0.007];

% EKF
Q = diag([100 100 100 100 0.1 0.1 0.1])*1.e-1;
P0 = diag([1e3 1e3 1e3 1e3 100 100 100])*1.e-4;
q0(4) = 0.5*sqrt(1+Abn0(1,1)+Abn0(2,2)+Abn0(3,3));
q0(1) = 1/(4*q0(4))*(Abn0(2,3)-Abn0(3,2));
q0(2) = 1/(4*q0(4))*(Abn0(3,1)-Abn0(1,3));
q0(3) = 1/(4*q0(4))*(Abn0(1,2)-Abn0(2,1));
R = [(ss_accuracy*pi/180)^2
     (eh_accuracy*pi/180)^2];
 
%% Run simulation
clear ans
sim('earth_pointing_ekf')

%%
wx_filt = squeeze(ans.omega_filt.Data(1,1,:)); 
wy_filt = squeeze(ans.omega_filt.Data(2,1,:));
wz_filt = squeeze(ans.omega_filt.Data(3,1,:));
wx = ans.omega.Data(:,1); 
wy = ans.omega.Data(:,2);
wz = ans.omega.Data(:,3);
Meffx_r = ans.Meff_r.Data(:,1);
Meffy_r = ans.Meff_r.Data(:,2);
Meffz_r = ans.Meff_r.Data(:,3);
SRP_x = squeeze(ans.SRP.Data(1,1,:));
SRP_y = squeeze(ans.SRP.Data(1,2,:));
SRP_z = squeeze(ans.SRP.Data(1,3,:));
SRP = [SRP_x, SRP_y, SRP_z];
mag_x = squeeze(ans.mag.Data(1,1,:));
mag_y = squeeze(ans.mag.Data(1,2,:));
mag_z = squeeze(ans.mag.Data(1,3,:));
Abn = ans.Abn.Data;
Aln = ans.Aln.Data;
mag = [mag_x, mag_y, mag_z];
clear t
t = linspace(0,t_ep,length(wx));

%%
% Errors
Ae = ans.Ae.Data;
q = zeros(length(t),4);
for i = 1:length(t)
    q(i,4) = 0.5*sqrt(1+Ae(1,1,i)+Ae(2,2,i)+Ae(3,3,i));
    q(i,1) = 1/(4*q(i,4))*(Ae(2,3,i)-Ae(3,2,i));
    q(i,2) = 1/(4*q(i,4))*(Ae(3,1,i)-Ae(1,3,i));
    q(i,3) = 1/(4*q(i,4))*(Ae(1,2,i)-Ae(2,1,i));
end


% %% From workspace, avoiding the simulation
% wx_filt = squeeze(out.omega_filt.Data(1,1,:)); 
% wy_filt = squeeze(out.omega_filt.Data(2,1,:));
% wz_filt = squeeze(out.omega_filt.Data(3,1,:));
% wx = out.omega.Data(:,1); 
% wy = out.omega.Data(:,2);
% wz = out.omega.Data(:,3);
% Meffx_r = out.Meff_r.Data(:,1);
% Meffy_r = out.Meff_r.Data(:,2);
% Meffz_r = out.Meff_r.Data(:,3);
% SRP_x = squeeze(out.SRP.Data(1,1,:));
% SRP_y = squeeze(out.SRP.Data(1,2,:));
% SRP_z = squeeze(out.SRP.Data(1,3,:));
% SRP = [SRP_x, SRP_y, SRP_z];
% mag_x = squeeze(out.mag.Data(1,1,:));
% mag_y = squeeze(out.mag.Data(1,2,:));
% mag_z = squeeze(out.mag.Data(1,3,:));
% Abn = out.Abn.Data;
% Aln = out.Aln1.Data;
% mag = [mag_x, mag_y, mag_z];
% clear t
% t = linspace(0,t_ep,length(wx));
% 
% Ae = out.Ae.Data;
% q = zeros(length(wx),4);
% for i = 1:length(t)
%     q(i,4) = 0.5*sqrt(1+Ae(1,1,i)+Ae(2,2,i)+Ae(3,3,i));
%     q(i,1) = 1/(4*q(i,4))*(Ae(2,3,i)-Ae(3,2,i));
%     q(i,2) = 1/(4*q(i,4))*(Ae(3,1,i)-Ae(1,3,i));
%     q(i,3) = 1/(4*q(i,4))*(Ae(1,2,i)-Ae(2,1,i));
% end
% 
% we = ([wx wy wz]-w_target)*360/pi;
% 
% %% Earth pointing plots
% 
% % angular velocities
% figure(9)
% settings
% 
% plot(t,we(:,1),t,we(:,2),t,we(:,3))
% hold on
% xline(380)
% hold on
% xline(5084)
% hold off
% title('Error on angular velocity');
% grid minor
% axis tight
% ylabel('$\Delta \omega_e$ $[^\circ/s]$'), xlabel('$t$ [s]')
% legend('$\omega_x$','$\omega_y$','$\omega_z$','Interpreter','latex');
% Lgnd = legend('show');
% Lgnd.Position(1) = 0.87;
% Lgnd.Position(2) = 0.84;

% angular velocities
figure(9)
settings

plot(t,wx(:,1),t,wy(:,2),t,wz(:,3))
title('Angular velocity');
grid minor
axis tight
ylabel('$\omega$ $[rad/s]$'), xlabel('$t$ [s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$','Interpreter','latex');
Lgnd = legend('show');
Lgnd.Position(1) = 0.87;
Lgnd.Position(2) = 0.84;

% Reaction wheels
figure(10)
settings
title('Reaction wheels');
plot(t,Meffx_r,t,Meffy_r,t,Meffz_r)
grid minor
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')
legend('$M_{effx}$','$M_{effy}$','$M_{effz}$');
Lgnd = legend('show');
Lgnd.Position(1) = 0.83;
Lgnd.Position(2) = 0.83;

% Perturbations
figure(11)
settings

subplot(2,1,1)
plot(t,SRP);
title('SRP')
grid minor
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')

subplot(2,1,2)
plot(t,mag)
title('Mag.field')
grid minor 
axis tight
ylabel('$M$ [Nmm]'), xlabel('$t$ [s]')
Lgnd.Position(1) = 0.86;
Lgnd.Position(2) = 0.83;

% Aln-Abn
figure(12)
settings

subplot(1,2,1)
Aln1 = squeeze(Aln(1,:,:));
Aln2 = squeeze(Aln(2,:,:));
Aln3 = squeeze(Aln(3,:,:));
plot(t,Aln1,t,Aln2,t,Aln3);
title('Aln')
grid minor
axis tight
xlabel('$t$ [s]')

subplot(1,2,2)
Abn1 = squeeze(Abn(1,:,:));
Abn2 = squeeze(Abn(2,:,:));
Abn3 = squeeze(Abn(3,:,:));
plot(t,Abn1,t,Abn2,t,Abn3);
title('Abn')
grid minor 
axis tight
xlabel('$t$ [s]')

% Quaternion error
figure(13)
settings

subplot(2,2,1)
plot(t,q(:,1));
grid minor
axis tight
ylabel('$q_1$'), xlabel('$t$ [s]')

subplot(2,2,2)
plot(t,q(:,2));
grid minor
axis tight
ylabel('$q_2$'), xlabel('$t$ [s]')

subplot(2,2,3)
plot(t,q(:,3));
grid minor
axis tight
ylabel('$q_3$'), xlabel('$t$ [s]')

subplot(2,2,4)
plot(t,(q(:,4)-1));
grid minor
axis tight
ylabel('$q_4$'), xlabel('$t$ [s]')
sgtitle('Pointing error on quaternions')
Lgnd.Position(1) = 0.86;
Lgnd.Position(2) = 0.83;

disp('Earth pointing is over')

%% Functions

function settings
    % Settings for the plots
    set(0,'defaulttextinterpreter','latex')
    set(0,'defaultlegendinterpreter','latex')
    set(0,'defaultAxesTickLabelInterpreter','latex')
    set(0,'DefaultTextFontSize',18)
    set(0,'DefaultAxesFontSize',14)
    set(0,'DefaultTextFontName','Times')
    set(0,'DefaultAxesFontName','Times')
    set(0,'DefaultLineLinewidth',1.2)
end


