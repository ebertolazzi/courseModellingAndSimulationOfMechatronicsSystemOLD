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
    gravity;
  end
  methods
    function self = Pendulum( ell, gravity )
      self@ODEbaseClass('Pendulum');
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
      ode(1) = u;
      ode(2) = v;
      t1     = (mu * x);
      ode(3) = -t1;
      t2     = (mu * y);
      ode(4) = -t2 - g;
      t10    = (x ^ 2);
      t11    = (y ^ 2);
      ode(5) = 1 / (t10 + t11) * (-3 * v * g - 4 * u * t1 - 4 * v * t2);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function jac = DfDx( self, t, Z )
      x        = Z(1);
      y        = Z(2);
      u        = Z(3);
      v        = Z(4);
      mu       = Z(5);
      g        = self.gravity;
      jac      = zeros(5,5);
      jac(1,3) = 1;
      jac(2,4) = 1;
      jac(3,1) = -mu;
      jac(3,5) = -x;
      jac(4,2) = jac(3,1);
      jac(4,5) = -y;
      t1       = (u * mu);
      t2       = (x ^ 2);
      t5       = mu * y;
      t11      = (y ^ 2);
      t15      = t2 + t11;
      t16      = t15 ^ 2;
      t17      = 0.1e1 / t16;
      jac(5,1) = t17 * ((4 * t2 * t1) + 0.6e1 * x * (0.4e1 / 0.3e1 * t5 + g) * v - (4 * t11 * t1));
      t21      = v * mu;
      jac(5,2) = t17 * (0.6e1 * g * v * y + 0.8e1 * x * y * t1 + 0.4e1 * t11 * t21 - 0.4e1 * t2 * t21);
      t31      = 0.1e1 / t15;
      jac(5,3) = -0.4e1 * t31 * x * mu;
      jac(5,4) = t31 * (-0.4e1 * t5 - 0.3e1 * g);
      jac(5,5) = t31 * (-0.4e1 * u * x - 0.4e1 * v * y);
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
    function H = hidden( self, t, Z )
      x    = Z(1);
      y    = Z(2);
      u    = Z(3);
      v    = Z(4);
      mu   = Z(5);
      g    = self.gravity;
      L    = self.ell;
      t1   = L ^ 2;
      t2   = (x ^ 2);
      t3   = (y ^ 2);
      H(1) = -t1 + t2 + t3;
      H(2) = 2 * u * x + 2 * v * y;
      t10  = u ^ 2;
      t11  = v ^ 2;
      H(3) = -2 * y * g - 2 * t2 * mu - 2 * t3 * mu + 2 * t10 + 2 * t11;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Hjac = DhiddenDx( self, t, Z )
      x    = Z(1);
      y    = Z(2);
      u    = Z(3);
      v    = Z(4);
      mu   = Z(5);
      g    = self.gravity;
      L    = self.ell;
      Hjac = zeros(3,5);
      Hjac(1,1) = 2 * x;
      Hjac(1,2) = 2 * y;
      Hjac(2,1) = 2 * u;
      Hjac(2,2) = 2 * v;
      Hjac(2,3) = Hjac(1,1);
      Hjac(2,4) = Hjac(1,2);
      Hjac(3,1) = -4 * mu * x;
      Hjac(3,2) = -4 * mu * y - 2 * g;
      Hjac(3,3) = 4 * u;
      Hjac(3,4) = 4 * v;
      t6        = x ^ 2;
      t7        = y ^ 2;
      Hjac(3,5) = -2 * t6 - 2 * t7;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [c,ceq,gc,gceq] = const( self, t, Z )
      c   = [];
      ceq = self.hidden( t, Z );
      if nargout > 2
        gc   = [];
        gceq = self.DhiddenDx( t, Z ).';
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [c,ceq,gc,gceq] = Pconst( self, t, Z )
      H   = self.hidden( t, Z );
      c   = [];
      ceq = H(1);
      if nargout > 2
        DHDx = self.DhiddenDx( t, Z );
        gc   = [];
        gceq = DHDx(1,1:2).';
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [c,ceq,gc,gceq] = Vconst( self, t, Z )
      H   = self.hidden( t, Z );
      c   = [];
      ceq = H(2);
      if nargout > 2
        DHDx = self.DhiddenDx( t, Z );
        gc   = [];
        gceq = DHDx(2,3:4).';
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [c,ceq,gc,gceq] = Lconst( self, t, Z )
      H   = self.hidden( t, Z );
      c   = [];
      ceq = H(3);
      if nargout > 2
        DHDx = self.DhiddenDx( t, Z );
        gc   = [];
        gceq = DHDx(3,5);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = Qtarget( self, X, X0 )
      f = 0.5*dot(X(:)-X0(:),X(:)-X0(:));
      if nargout > 1
        g = X(:)-X0(:);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Zproj = project( self, t, Z )
      Zproj = Z(:);
      x0    = Z(1);
      y0    = Z(2);
      u0    = Z(3);
      v0    = Z(4);
      mu0   = Z(5);
      options = optimoptions('fmincon',...
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true, ...
        'CheckGradients', false, ...
        'Display','off' ...
      );
      % project position (x,y)
      fun     = @(P) self.Qtarget( P,[x0;y0]);
      nonlcon = @(P) self.Pconst( t, [P(:);Zproj(3:5)] );
      res     = fmincon(...
        fun, Z(1:2), [], [], [], [], [-Inf;-Inf], [Inf;Inf], ...
        nonlcon, options ...
      );
      Zproj(1) = res(1);
      Zproj(2) = res(2);
      % project velocity (u,v)
      fun     = @(V) self.Qtarget(V,[u0;v0]);
      nonlcon = @(V) self.Vconst( t, [Zproj(1:2);V(:);mu0] );
      res     = fmincon(...
        fun, Z(3:4), [], [], [], [], [-Inf;-Inf], [Inf;Inf], ...
        nonlcon, options ...
      );
      Zproj(3) = res(1);
      Zproj(4) = res(2);
      % project multiplier
      if true
        g  = self.gravity;
        L  = self.ell;
        t2 = Zproj(3)^ 2;
        t3 = Zproj(4) ^ 2;
        t5 = Zproj(1)^ 2;
        t6 = Zproj(2)^ 2;
        Zproj(5) = -0.1e1 / (t5 + t6) * (g * Zproj(2) - t2 - t3);
      else
        fun     = @(mu) self.Qtarget(mu,mu0);
        nonlcon = @(mu) self.Lconst( t, [Zproj(1:4),mu] );
        Zproj(5) = fmincon( fun, mu0, [], [], [], [], [-Inf], [Inf], nonlcon, options );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Zproj = project2( self, t, Z )
      options = optimoptions('fmincon',...
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true, ...
        'CheckGradients', false, ...
        'Display','iter' ...
      );
      % project position (x,y)
      fun     = @(ZZ) self.Qtarget( ZZ, Z );
      nonlcon = @(ZZ) self.const( t, ZZ );
      Zproj   = fmincon(...
        fun, Z, [], [], [], [], -Inf*ones(1,5), Inf*ones(1,5), ...
        nonlcon, options ...
      );
    end
  end
end
