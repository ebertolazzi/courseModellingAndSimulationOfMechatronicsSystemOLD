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
classdef Pendulum2 < handle
  properties (SetAccess = protected, Hidden = true)
    mass;
    gravity;
  end
  methods
    function self = Pendulum2( mass, gravity )
      self.mass    = mass;
      self.gravity = gravity;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Z = eval( self, t, P, V )
      x  = P(1); y  = P(2);
      vx = V(1); vy = V(2);
      M = [ self.mass, 0,         2*x ; ...
            0,         self.mass, 2*y ; ...
            2*x,       2*y,         0 ];
      rhs = [ 0; -self.mass *self.gravity; -2*(vx^2+vy^2) ];
      Z   = M\rhs;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [P1,V1] = step( self, t, P0, V0, dt )
      Z  = self.eval( t, P0, V0 );
      P1 = P0 + dt * ( V0 + (dt/2) * Z(1:2) );
      V1 = V0 + dt * Z(1:2);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [P,V] = advance( self, t, P0, V0 )
      P = zeros(2,length(t));
      V = zeros(2,length(t));
      P(:,1) = P0;
      V(:,1) = V0;
      for k=2:length(t)
        dt = t(k)-t(k-1);
        [P(:,k),V(:,k)] = step( self, t(k-1), P(:,k-1), V(:,k-1), dt );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
