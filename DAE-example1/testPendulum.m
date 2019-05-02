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
ell = 2;
g   = 9.8;
ode = Pendulum( ell, g );

%
% tell  solver that the ODE to be used is the Pendulum ode
solver.setODE(ode);
solver_DAE.setODE(ode);

% set the time steps 
tt  = 0:0.05:20;

% now choose the initial condition (not consistent)
ini = [ ell, 0, 0, -1, 0 ];
% project initial condition to the hidden constraints
ini = ode.project( 0, ini ) ;
%ini = ode.project2( 0, ini ) ;

fprintf('compute solution, ODE with drift\n');
xy     = solver.advance( tt, ini );
fprintf('compute solution, ODE with projection\n');
xy_DAE = solver_DAE.advance_DAE( tt, ini );

% Extract the solution
x      = xy(1,:);
y      = xy(2,:);
u      = xy(3,:);
v      = xy(4,:);
lambda = xy(5,:);

subplot( 2, 1, 1);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Solution of ODE (not stabilized)');

% Extract the solution
x      = xy_DAE(1,:);
y      = xy_DAE(2,:);
u      = xy_DAE(3,:);
v      = xy_DAE(4,:);
lambda = xy_DAE(5,:);

subplot( 2, 1, 2);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Solution of ODE with projection');
