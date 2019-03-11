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
classdef VanDerPol < ODEbaseClass
  properties (SetAccess = protected, Hidden = true)
    mu;
  end
  methods
    function self = VanDerPol( mu )
      self@ODEbaseClass('VanDerPol');
      self.mu = mu;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f( self, t, x )
      y = x(1);
      v = x(2);
      res = [v;self.mu*(1-y^2)*v-y];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDx( self, t, x )
      y = x(1);
      v = x(2);
      res = [0, 1;-1-2*self.mu*y,self.mu*(1-y^2)];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDt( self, t, x )
      res = [0;0];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t0, x0, t )
      res  = [0;0]*t(:).';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
