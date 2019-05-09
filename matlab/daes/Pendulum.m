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
classdef Pendulum < DAEbaseClass
  properties (SetAccess = protected, Hidden = true)
    ell;
    gravity;
  end
  methods
    function self = Pendulum( ell, gravity )
      self@DAEbaseClass('Pendulum');
      self.ell     = ell;
      self.gravity = gravity;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ode = f( self, t, Z )
      x      = Z(1);
      y      = Z(2);
      u      = Z(3);
      v      = Z(4);
      mu     = Z(5);
      g      = self.gravity;
      ode    = zeros(5,1);
      ode_1 = u;
      ode_2 = v;
      t1 = (mu * x);
      ode_3 = -t1;
      ode_4 = -mu * y - g;
      t11 = (x ^ 2);
      t12 = (y ^ 2);
      ode_5 = 1 / (t11 + t12) * (-4 * y * v * mu - 3 * v * g - 4 * u * t1);
      ode(1) = ode_1;
      ode(2) = ode_2;
      ode(3) = ode_3;
      ode(4) = ode_4;
      ode(5) = ode_5;      
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function jac = DfDx( self, t, Z )
      x       = Z(1);
      y       = Z(2);
      u       = Z(3);
      v       = Z(4);
      mu      = Z(5);
      g       = self.gravity;
      jac     = zeros(5,5);
      jac_1_3 = 1;
      jac_2_4 = 1;
      jac_3_1 = -mu;
      jac_3_5 = -x;
      jac_4_2 = jac_3_1;
      jac_4_5 = -y;
      t1 = (u * mu);
      t2 = (x ^ 2);
      t5 = mu * y;
      t11 = (y ^ 2);
      t15 = t2 + t11;
      t16 = t15 ^ 2;
      t17 = 0.1e1 / t16;
      jac_5_1 = t17 * ((4 * t2 * t1) + 0.6e1 * x * (0.4e1 / 0.3e1 * t5 + g) * v - (4 * t11 * t1));
      t21 = v * mu;
      jac_5_2 = t17 * (0.6e1 * g * v * y + 0.8e1 * x * y * t1 + 0.4e1 * t11 * t21 - 0.4e1 * t2 * t21);
      t31 = 0.1e1 / t15;
      jac_5_3 = -0.4e1 * t31 * x * mu;
      jac_5_4 = t31 * (-0.4e1 * t5 - 0.3e1 * g);
      jac_5_5 = t31 * (-0.4e1 * u * x - 0.4e1 * v * y);
      jac(1,3) = jac_1_3;
      jac(2,4) = jac_2_4;
      jac(3,1) = jac_3_1;
      jac(3,5) = jac_3_5;
      jac(4,2) = jac_4_2;
      jac(4,5) = jac_4_5;
      jac(5,1) = jac_5_1;
      jac(5,2) = jac_5_2;
      jac(5,3) = jac_5_3;
      jac(5,4) = jac_5_4;
      jac(5,5) = jac_5_5;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDt( self, t, x )
      res = zeros(5,1);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t0, x0, t )
      res = zeros(5,1)*t(:).';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = algebraic( self, t, Z )
      x    = Z(1);
      y    = Z(2);
      u    = Z(3);
      v    = Z(4);
      mu   = Z(5);
      g    = self.gravity;
      L    = self.ell;
      t1 = L ^ 2;
      t2 = (x ^ 2);
      t3 = (y ^ 2);
      hidden_1 = -t1 + t2 + t3;
      hidden_2 = 2 * u * x + 2 * v * y;
      t10 = u ^ 2;
      t11 = v ^ 2;
      hidden_3 = -2 * y * g - 2 * t2 * mu - 2 * mu * t3 + 2 * t10 + 2 * t11;
      res = zeros(3,1);
      res(1) = hidden_1;
      res(2) = hidden_2;
      res(3) = hidden_3;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function hjac = DalgebraicDx( self, t, Z )
      x    = Z(1);
      y    = Z(2);
      u    = Z(3);
      v    = Z(4);
      mu   = Z(5);
      g    = self.gravity;
      L    = self.ell;
      hjac = zeros(3,5);
      hjac_1_1 = 2 * x;
      hjac_1_2 = 2 * y;
      hjac_2_1 = 2 * u;
      hjac_2_2 = 2 * v;
      hjac_2_3 = hjac_1_1;
      hjac_2_4 = hjac_1_2;
      hjac_3_1 = -4 * mu * x;
      hjac_3_2 = -4 * mu * y - 2 * g;
      hjac_3_3 = 4 * u;
      hjac_3_4 = 4 * v;
      t6 = x ^ 2;
      t7 = y ^ 2;
      hjac_3_5 = -2 * t6 - 2 * t7;
      hjac(1,1) = hjac_1_1;
      hjac(1,2) = hjac_1_2;
      hjac(2,1) = hjac_2_1;
      hjac(2,2) = hjac_2_2;
      hjac(2,3) = hjac_2_3;
      hjac(2,4) = hjac_2_4;
      hjac(3,1) = hjac_3_1;
      hjac(3,2) = hjac_3_2;
      hjac(3,3) = hjac_3_3;
      hjac(3,4) = hjac_3_4;
      hjac(3,5) = hjac_3_5;      
    end
  end
end
