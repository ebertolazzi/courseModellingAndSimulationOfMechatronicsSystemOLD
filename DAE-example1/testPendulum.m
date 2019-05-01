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
solver = Heun();

% initialize to variable ode with the Pendulum
% ODE object instance
ell = 2;
g   = 9.8;
ode = Pendulum( ell, g );

%
% tell  solver that the ODE to be used is the Pendulum ode
solver.setODE(ode);

% set the time steps 
tt  = 0:0.01:20;

% now choose the initial condition (not consistent)
ini = [ ell, 0, 0, -1, 0 ];
%ini = ode.project( 0, ini ) ; 
ini = ode.project2( 0, ini ) ; 

xy = solver.advance( tt, ini );

% Extract the solution
x      = xy(1,:);
y      = xy(2,:);
u      = xy(3,:);
v      = xy(4,:);
lambda = xy(5,:);

plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('The solution');
