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
classdef Pendulum < ODEbaseClass
  properties (SetAccess = protected, Hidden = true)
    ell;
    mass;
    gravity;
  end
  methods
    function self = Pendulum( ell, mass, gravity )
      self@ODEbaseClass('Pendulum');
      self.ell     = ell;
      self.mass    = mass;
      self.gravity = gravity;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Zp = f( self, t, Z )
      x      = Z(1);
      y      = Z(2);
      u      = Z(3);
      v      = Z(4);
      lambda = Z(5);
      Zp     = [ u; ...
                 v; ...
                 -lambda*x/self.mass; ...
                 -lambda*y/self.mass-self.gravity; ...
                 (-4*lambda*(x*u+y*v)-3*v*self.mass*self.gravity)/(x^2+y^2)];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDx( self, t, Z )
      x      = Z(1);
      y      = Z(2);
      u      = Z(3);
      v      = Z(4);
      lambda = Z(5);
      t1 = (-self.gravity * self.mass);
      t5 = (u * lambda);
      t6 = (x ^ 2);
      t9 = (y ^ 2);
      t12 = (lambda * v);
      t13 = (x * y);
      t17 = t6 + t9;
      t18 = t17 ^ 2;
      t19 = 1 / t18;
      c1 = t19 * (-6 * v * x * t1 + 8 * t13 * t12 + 4 * t6 * t5 - 4 * t9 * t5);
      t20 = v * y;
      c2 = t19 * (-6 * t20 * t1 - 4 * t6 * t12 + 4 * t9 * t12 + 8 * t13 * t5);
      t31 = 1 / t17;
      c3 = -4 * t31 * x * lambda;
      c4 = t31 * (-4 * lambda * y + 3 * t1);
      c5 = t31 * (-4 * u * x - 4 * t20);

      res = [ 0, 0, 1, 0, 0; ...
              0, 0, 0, 1, 0; ...
              -lambda/self.mass, 0, 0, 0, -x/self.mass; ...
              0, -lambda/self.mass, 0, 0, -y/self.mass; ...
              c1, c2, c3, c4, c5];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDt( self, t, x )
      res = zeros(5,1);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t0, x0, t )
      res  = zeros(5,1)*t(:).';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
