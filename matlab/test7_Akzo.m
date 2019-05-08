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

%e1   = ExplicitEuler();
%e1   = LobattoIIIC();
e1   = GaussLegendre4();
ode1 = Akzo();
e1.setODE(ode1);

tt  = 0:0.1:200;
% setup initial condition
y1 = 0.444;
y4 = 0.007;
Ks = 115.83;
ini = [ y1, 0.00123, 0, y4, 0, Ks*y1*y4];
sol = e1.advance( tt, ini );

for kk=1:6
  subplot(3,2,kk);
  %plot( tt, sol(kk,:), '-o', 'MarkerSize', 2, 'Linewidth', 2 );
  plot( tt, sol(kk,:), 'Linewidth', 2 );
end