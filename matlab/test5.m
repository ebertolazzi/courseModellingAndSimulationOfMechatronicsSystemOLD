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

e1   = RadauIA();
e2   = RadauIIA();
e3   = LobattoIIIB();
e4   = LobattoIIIC();

ode1 = VanDerPol(5);

e1.setODE(ode1);
e2.setODE(ode1);
e3.setODE(ode1);
e4.setODE(ode1);
%e5.setODE(ode1);

tt  = 0:0.01:100;
ini = [0;1];
disp('advance RadauIA');
xy1 = e1.advance( tt, ini );
disp('advance RadauIIA');
xy2 = e2.advance( tt, ini );
disp('advance LobattoIIIB');
xy3 = e3.advance( tt, ini );
disp('advance LobattoIIIC');
xy4 = e4.advance( tt, ini );
%xy5 = e5.advance( tt, ini );

subplot(2,2,1)
plot( tt, xy1(1,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
%axis equal
title(e1.getName());

subplot(2,2,2)
plot( tt, xy2(1,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
%axis equal
title(e2.getName());

subplot(2,2,3)
plot( tt, xy3(1,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
%axis equal
title(e3.getName());

subplot(2,2,4)
plot( tt, xy4(1,:), '-o', 'MarkerSize', 6, 'Linewidth', 2 );
%axis equal
title(e4.getName());
