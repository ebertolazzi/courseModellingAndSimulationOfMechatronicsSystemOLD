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
%
% Compare with results in:
%
% https://archimede.dm.uniba.it/~testset/report/chemakzo.pdf
%
addpath('explicit');
addpath('implicit');
addpath('odes');
addpath('daes');

close all;

solver     = Heun();
solver_DAE = Heun_DAE();
ell        = 2;
gravity    = 9.81;
dae        = Pendulum( ell, gravity );

solver.setODE(dae);
solver_DAE.setODE(dae);

tt  = 0:0.1:5;
% setup initial condition
x0      = ell;
y0      = 0;
u0      = 0;
v0      = -1;
lambda0 = 0;
ini     = [x0;y0;u0;v0;lambda0];
ini     = dae.project(0,ini);
fprintf('avance with ODE and possible drift\n');
sol = solver.advance( tt, ini );
fprintf('avance with ODE+PROJECTION\n');
sol_DAE = solver_DAE.advance_DAE( tt, ini );
fprintf('done\n');

subplot(2,1,1);
x = sol(1,:);
y = sol(2,:);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
hold on
axis equal
xx = ell*cos(0:pi/100:2*pi);
yy = ell*sin(0:pi/100:2*pi);
plot( xx, yy, '-r', 'Linewidth', 2 );
title('ODE+no stabilization');

subplot(2,1,2);
x = sol_DAE(1,:);
y = sol_DAE(2,:);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
hold on
axis equal
plot( xx, yy, '-r', 'Linewidth', 2 );
title('ODE+projection');
