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
eta     = 0.5;
omega   = 1000;
solver  = Pendulum2( mass, gravity );
solver2 = PendulumStab( eta, omega, mass, gravity );

% set the time steps
tt  = 0:0.01:10;

% now choose the initial condition (not consistent)
P0 = [ 1; 0];
V0 = [ 0; 0];

[P,V]   = solver.advance( tt, P0, V0 );
[PP,VV] = solver2.advance( tt, P0, V0 );

% Extract the solution
x = P(1,:);
y = P(2,:);
u = V(1,:);
v = V(2,:);

xx = PP(1,:);
yy = PP(2,:);
uu = VV(1,:);
vv = VV(2,:);

subplot( 2, 1, 1);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('The solution');

subplot( 2, 1, 2);
plot( xx, yy, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('The solution');
