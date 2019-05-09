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
classdef Heun_DAE < Heun
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    % Heun tableau
    %
    % 0 | 0   0
    % 1 | 1   0
    % --+-------
    %   | 1/2 1/2
    %
    function self = Heun_DAE( )
      self@Heun();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %  for the explicit RK methods
    %
    %  x(k+1) = x(k) + sum(j) b(j) * K(k)
    %
    %  K(i) = h*f( t(k) + c(i) * dt, x(k) + sum(j) A(i,j) K(j) ), i=1..s
    %
    function x1 = step_DAE( self, t0, x0, dt )
      x1 = self.step( t0, x0, dt );
      x1 = self.odeClass.project( t0+dt, x1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    %  compute approximate solution on a series of point given by the vector t
    %  starting at initial point x0
    %
    function x = advance_DAE( self, t, x0 )
      x      = zeros(length(x0),length(t));
      x(:,1) = x0(:);
      for k=1:length(t)-1
        x(:,k+1) = self.step_DAE( t(k), x(:,k), t(k+1)-t(k) );
      end
    end
  end
end
