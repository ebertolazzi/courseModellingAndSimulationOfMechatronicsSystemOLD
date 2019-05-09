%
% Matlab code for the Course:
%
%     Modelling and Simulation Mechatronics System
%
% by
% Enrico Bertolazzi
% Dipartimento di Ingegneria Industriale
% Universita` degli Studi di Trento
% email: enrico.bertolazzi@unitn.it
%
addpath('../matlab');
addpath('../matlab/explicit');
addpath('../matlab/implicit');

% initialize to variable solver with the Explicit Euler
% solver object instance
%solver = ExplicitEuler();
%solver = ImplicitEuler();
solver     = Heun();
solver_DAE = Heun_DAE();

% initialize to variable ode with the Pendulum
% ODE object instance
ode = SliderCrank();

%
% tell  solver that the ODE to be used is the Pendulum ode
solver.setODE(ode);
solver_DAE.setODE(ode);

% set the time steps 
tt  = 0:0.000001:0.00001;

% now choose the initial condition (not consistent)
ini = [ 60/180*pi, -30/180*pi, 0, 0, 0 ];
% project initial condition to the hidden constraints
ini = ode.project( 0, ini ) ;
%ini = ode.project2( 0, ini ) ;

fprintf('compute solution, ODE with drift\n');
xy     = solver.advance( tt, ini );

% Extract the solution
theta1     = xy(1,:);
theta2     = xy(2,:);
theta1_dot = xy(3,:);
theta2_dot = xy(4,:);
lambda     = xy(5,:);

subplot( 2, 1, 1);
plot( tt, theta1, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Solution of ODE (not stabilized)');

stop
fprintf('compute solution, ODE with projection\n');
xy_DAE = solver_DAE.advance_DAE( tt, ini );


% Extract the solution
theta1     = xy(1,:);
theta2     = xy(2,:);
theta1_dot = xy(3,:);
theta2_dot = xy(4,:);
lambda     = xy(5,:);

subplot( 2, 1, 2);
plot( tt, theta1, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Solution of ODE with projection');
