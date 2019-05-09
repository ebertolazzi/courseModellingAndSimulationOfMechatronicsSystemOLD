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
dae        = ParabolicPendulum();

solver.setODE(dae);
solver_DAE.setODE(dae);

tt  = 0:0.025:25;
% setup initial condition
x0      = 0;
y0      = 1;
u0      = 1;
v0      = 0;
lambda0 = (1/4)*9.81-3/2;
ini = [x0;y0;u0;v0;lambda0];
fprintf('avance with ODE and possible drift\n');
sol= solver.advance( tt, ini );
fprintf('avance with ODE+PROJECTION\n');
sol_DAE = solver_DAE.advance_DAE( tt, ini );
fprintf('done\n');

subplot(2,1,1);
x = sol(1,:);
y = sol(2,:);
plot( x, y, '-o', 'MarkerSize', 6, 'Linewidth', 2 );
hold on
axis equal
xx = -1:0.01:1;
yy = xx.^2+sqrt(1-xx.^2);
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
