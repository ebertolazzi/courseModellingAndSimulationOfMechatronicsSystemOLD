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
classdef ODE1 < ODEbaseClass
  %% MATLAB class wrapper for the underlying C++ class
  methods
    function self = ODE1()
      self@ODEbaseClass('ODE1');
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f( self, t, x )
      x1  = x(1);
      x2  = x(2);
      res = [-x2; x1];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDx( self, t, x )
      x1  = x(1);
      x2  = x(2);
      res = [0, -1; 1, 0];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDt( self, t, x )
      res = [0;0];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t0, x0, t )
      R    = hypot(x0(2),x0(1));
      phi0 = atan2(x0(2),x0(1));
      res  = R*[cos(t(:)-t0+phi0),sin(t(:)-t0+phi0)].';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
