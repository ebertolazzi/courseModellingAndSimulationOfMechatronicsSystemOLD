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
addpath('explicit');
addpath('implicit');
addpath('odes');
addpath('daes');

%e1   = ExplicitEuler();
%e1   = LobattoIIIC();
e1  = GaussLegendre4();
ode1 = ParabolicPendulum();
e1.setODE(ode1);

tt  = 0:0.025:500;
% setup initial condition
x0      = 0;
y0      = 1;
u0      = 1;
v0      = 0;
lambda0 = (1/4)*9.81-3/2;
ini = [x0;y0;u0;v0;lambda0];
sol = e1.advance( tt, ini );
x = sol(1,:);
y = sol(2,:);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
hold on
axis equal
xx = -1:0.01:1;
yy = xx.^2+sqrt(1-xx.^2);
plot( xx, yy, '-r', 'Linewidth', 2 );
x = sol(1,:);
title(e1.getName());
