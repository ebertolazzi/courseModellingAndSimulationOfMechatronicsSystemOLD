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

close all;

ee   = Collatz();
ei   = CrankNicolson();
ode1 = ODE1();

ee.setODE(ode1);
ei.setODE(ode1);

tt  = 0:pi/10:2*pi;
ini = [0;1];
xy1 = ee.advance( tt, ini );
xy2 = ei.advance( tt, ini );

plot( xy1(1,:), xy1(2,:), '-o', 'MarkerSize', 12, 'Linewidth', 3 );
hold on
plot( xy2(1,:), xy2(2,:), '--*', 'MarkerSize', 12, 'Linewidth', 3 );
axis equal