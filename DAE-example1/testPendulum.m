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
solver = ImplicitEuler();

% initialize to variable ode with the Pendulum
% ODE object instance
ell  = 2;
mass = 10;
ode = Pendulum( ell, mass, 9.8 );

%
% tell  solver that the ODE to be used is the Pendulum ode
solver.setODE(ode);

% set the time steps 
tt  = 0:0.01:10;

% now choose the initial condition (not consistent)
ini = [ ell, 0, 0, -1, mass/4  ]; 

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
