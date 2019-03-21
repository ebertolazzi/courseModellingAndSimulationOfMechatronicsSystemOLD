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

mass    = 10;
gravity = 9.8;
solver  = Pendulum2( mass, gravity );

% set the time steps
tt  = 0:0.001:100;

% now choose the initial condition (not consistent)
P0 = [ 1; 0];
V0 = [ 0; 0];

[P,V] = solver.advance( tt, P0, V0 );

% Extract the solution
x = P(1,:);
y = P(2,:);
u = V(1,:);
v = V(2,:);

plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('The solution');
