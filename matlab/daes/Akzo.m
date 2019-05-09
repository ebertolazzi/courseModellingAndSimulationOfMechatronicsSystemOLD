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
classdef Akzo < DAEbaseClass
  properties (SetAccess = protected, Hidden = true)
    k1
    k2
    k3
    k4
    K
    klA
    Ks
    pCO2
    H
  end
  methods
    function self = Akzo()
      self@DAEbaseClass('Akzo');
      self.k1     = 18.7;
      self.k2     = 0.58;
      self.k3     = 0.09;
      self.k4     = 0.42;
      self.K      = 34.4;
      self.klA    = 3.3;
      self.Ks     = 115.83;
      self.pCO2   = 0.9;
      self.H      = 737;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ode = f( self, t, Z )
      x1   = Z(1);
      x2   = Z(2);
      x3   = Z(3);
      x4   = Z(4);
      x5   = Z(5);
      x6   = Z(6);
      k1   = self.k1;
      k2   = self.k2;
      k3   = self.k3;
      k4   = self.k4;
      K    = self.K;
      klA  = self.klA;
      Ks   = self.Ks;
      pCO2 = self.pCO2;
      H    = self.H;
      % automatically generated by MAPLE
      t1 = x1 ^ 2;
      t2 = t1 ^ 2;
      t3 = t2 * k1;
      t4 = sqrt(x2);
      t5 = K * t4;
      t6 = t5 * t3;
      t8 = k3 * x1;
      t9 = x4 ^ 2;
      t10 = K * t9;
      t11 = t10 * t8;
      t12 = k2 * x3;
      t14 = x4 * K * t12;
      t16 = k2 * x1 * x5;
      t18 = 0.1e1 / K;
      ode_1 = t18 * (-0.2e1 * t6 - t11 + t14 - t16);
      t19 = x6 ^ 2;
      t20 = t19 * k4;
      ode_2 = 0.1e1 / H * (-t4 * (t3 + t20) * H + H * (-0.2e1 * klA * x2 - 0.2e1 * t9 * t8) + 0.2e1 * klA * pCO2) / 0.2e1;
      ode_3 = t18 * (t6 - t14 + t16);
      ode_4 = t18 * (-0.2e1 * t11 - t14 + t16);
      ode_5 = t18 * (t5 * t20 + t14 - t16);
      t48 = K * x1;
      ode_6 = -t18 * (0.2e1 * k1 * x4 * t4 * t2 * K + 0.2e1 * k3 * t9 * t1 * K + x3 * x4 * k2 * t48 + x1 * x4 * x5 * k2 + k3 * t9 * x4 * t48 - k2 * x5 * t1 - t12 * t10) * Ks;
      ode    = zeros(6,1);
      ode(1) = ode_1;
      ode(2) = ode_2;
      ode(3) = ode_3;
      ode(4) = ode_4;
      ode(5) = ode_5;
      ode(6) = ode_6;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function jac = DfDx( self, t, Z )
      x1   = Z(1);
      x2   = Z(2);
      x3   = Z(3);
      x4   = Z(4);
      x5   = Z(5);
      x6   = Z(6);
      k1   = self.k1;
      k2   = self.k2;
      k3   = self.k3;
      k4   = self.k4;
      K    = self.K;
      klA  = self.klA;
      Ks   = self.Ks;
      pCO2 = self.pCO2;
      H    = self.H;
      % automatically generated by MAPLE
      jac = zeros(6,6);
      t1 = x1 ^ 2;
      t2 = t1 * x1;
      t4 = sqrt(x2);
      t5 = K * t4;
      t6 = t5 * t2 * k1;
      t8 = x4 ^ 2;
      t9 = t8 * k3;
      t10 = K * t9;
      t11 = k2 * x5;
      t13 = 0.1e1 / K;
      jac_1_1 = t13 * (-0.8e1 * t6 - t10 - t11);
      t14 = t1 ^ 2;
      t15 = t14 * k1;
      t16 = 0.1e1 / t4;
      t17 = t16 * t15;
      jac_1_2 = -t17;
      jac_1_3 = k2 * x4;
      t18 = k3 * x1;
      t19 = (x4 * t18);
      t20 = 2 * t19;
      t21 = (k2 * x3);
      jac_1_4 = -t20 + t21;
      t22 = k2 * x1;
      t23 = t13 * t22;
      jac_1_5 = -t23;
      jac_2_1 = -0.2e1 * t2 * k1 * t4 - t9;
      t27 = x6 ^ 2;
      t28 = t27 * k4;
      jac_2_2 = -t16 * (0.4e1 * t4 * klA + t15 + t28) / 0.4e1;
      jac_2_4 = -t20;
      t35 = t4 * x6 * k4;
      jac_2_6 = -t35;
      jac_3_1 = t13 * (0.4e1 * t6 + t11);
      jac_3_2 = t17 / 0.2e1;
      jac_3_3 = -jac_1_3;
      jac_3_4 = -t21;
      jac_3_5 = t23;
      jac_4_1 = t13 * (-0.2e1 * t10 + t11);
      jac_4_3 = jac_3_3;
      jac_4_4 = -4 * t19 + jac_3_4;
      jac_4_5 = jac_3_5;
      jac_5_1 = -t13 * t11;
      jac_5_2 = t16 * t28 / 0.2e1;
      jac_5_3 = jac_1_3;
      jac_5_4 = t21;
      jac_5_5 = jac_1_5;
      jac_5_6 = 0.2e1 * t35;
      t52 = K * t8 * t18;
      t59 = x5 * t22;
      jac_6_1 = -0.8e1 * t13 * Ks * (k1 * x4 * t4 * t2 * K + k3 * t8 * x4 * K / 0.8e1 + t52 / 0.2e1 + k2 * (K * x3 + x5) * x4 / 0.8e1 - t59 / 0.4e1);
      jac_6_2 = -k1 * x4 * t16 * t14 * Ks;
      t69 = Ks * k2;
      t70 = x1 - x4;
      jac_6_3 = -t70 * x4 * t69;
      jac_6_4 = -t13 * (0.4e1 * x4 * k3 * t1 * K + jac_5_4 * K * x1 - 0.2e1 * x4 * K * jac_5_4 + 0.2e1 * t5 * t15 + 0.3e1 * t52 + t59) * Ks;
      jac_6_5 = t13 * t70 * x1 * t69;
      jac(1,1) = jac_1_1;
      jac(1,2) = jac_1_2;
      jac(1,3) = jac_1_3;
      jac(1,4) = jac_1_4;
      jac(1,5) = jac_1_5;
      jac(2,1) = jac_2_1;
      jac(2,2) = jac_2_2;
      jac(2,4) = jac_2_4;
      jac(2,6) = jac_2_6;
      jac(3,1) = jac_3_1;
      jac(3,2) = jac_3_2;
      jac(3,3) = jac_3_3;
      jac(3,4) = jac_3_4;
      jac(3,5) = jac_3_5;
      jac(4,1) = jac_4_1;
      jac(4,3) = jac_4_3;
      jac(4,4) = jac_4_4;
      jac(4,5) = jac_4_5;
      jac(5,1) = jac_5_1;
      jac(5,2) = jac_5_2;
      jac(5,3) = jac_5_3;
      jac(5,4) = jac_5_4;
      jac(5,5) = jac_5_5;
      jac(5,6) = jac_5_6;
      jac(6,1) = jac_6_1;
      jac(6,2) = jac_6_2;
      jac(6,3) = jac_6_3;
      jac(6,4) = jac_6_4;
      jac(6,5) = jac_6_5;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = DfDt( self, t, x )
      res = zeros(6,1);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t0, x0, t )
      res  = [0;0;0;0;0;0]*t(:).';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = algebraic( self, t, Z )
      x1   = Z(1);
      x2   = Z(2);
      x3   = Z(3);
      x4   = Z(4);
      x5   = Z(5);
      x6   = Z(6);
      k1   = self.k1;
      k2   = self.k2;
      k3   = self.k3;
      k4   = self.k4;
      K    = self.K;
      klA  = self.klA;
      Ks   = self.Ks;
      pCO2 = self.pCO2;
      H    = self.H;
      % automatically generated by MAPLE
      res = Ks * x1 * x4 - x6;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function hjac = DalgebraicDx( self, t, Z )
      x1   = Z(1);
      x2   = Z(2);
      x3   = Z(3);
      x4   = Z(4);
      x5   = Z(5);
      x6   = Z(6);
      k1   = self.k1;
      k2   = self.k2;
      k3   = self.k3;
      k4   = self.k4;
      K    = self.K;
      klA  = self.klA;
      Ks   = self.Ks;
      pCO2 = self.pCO2;
      H    = self.H;
      % automatically generated by MAPLE
      hjac     = zeros(1,6);
      hjac_1_1 = Ks * x4;
      hjac_1_4 = Ks * x1;
      hjac_1_6 = -1;
      hjac(1,1) = hjac_1_1;
      hjac(1,4) = hjac_1_4;
      hjac(1,6) = hjac_1_6;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
