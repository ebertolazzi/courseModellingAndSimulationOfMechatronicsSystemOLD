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

e1   = ExplicitEuler();
e2   = ImplicitEuler();
e3   = Collatz();
e4   = Heun();
e5   = CrankNicolson();
ode1 = ODE1();

e1.setODE(ode1);
e2.setODE(ode1);
e3.setODE(ode1);
e4.setODE(ode1);
e5.setODE(ode1);

tt1 = 0:pi/40:2*pi;
tt  = 0:pi/8:2*pi;
ini = [0;1]; 
xy1 = e1.advance( tt1, ini );
xy2 = e2.advance( tt1, ini );
xy3 = e3.advance( tt,  ini );
xy4 = e4.advance( tt,  ini );
xy5 = e5.advance( tt,  ini );

subplot(2,2,1)
plot( xy1(1,:), xy1(2,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
hold on
plot( xy2(1,:), xy2(2,:), '--*', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Eueler');

subplot(2,2,2)
plot( xy3(1,:), xy3(2,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Collatz');

subplot(2,2,3)
plot( xy4(1,:), xy4(2,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Heun');

subplot(2,2,4)
plot( xy5(1,:), xy5(2,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
axis equal
title('Crank-Nicolson');