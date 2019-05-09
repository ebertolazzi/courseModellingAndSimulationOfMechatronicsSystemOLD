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
classdef RadauIA_DAE < RadauIA
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = RadauIA_DAE( )
      self@RadauIA();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function x1 = step_DAE( self, t0, x0, dt )
      x1 = self.step( t0, x0, dt );
      x1 = self.odeClass.project( t0+dt, x1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function x = advance_DAE( self, t, x0 )
      x      = zeros(length(x0),length(t));
      x(:,1) = x0(:);
      for k=1:length(t)-1
        x(:,k+1) = self.step_DAE( t(k), x(:,k), t(k+1)-t(k) );
      end
    end
  end
end
