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
classdef DAEbaseClass < ODEbaseClass
  methods (Abstract)
    %
    %  Abstract functions defining the invariants of the DAE
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    algebraic( self, t, Z )
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DalgebraicDx( self, t, Z )
  end

  methods
    function self = DAEbaseClass( name )
      self@ODEbaseClass( name );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = Qtarget( self, X, X0 )
      dx = X(:)-X0(:);
      f  = 0.5*dot(dx,dx);
      if nargout > 1
        g = dx;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [c,ceq,Jc,Jceq] = Constraints( self, t, Z  )
      c   = [];
      ceq = self.algebraic( t, Z );
      if nargout > 2
        Jc   = [];
        Jceq = self.DalgebraicDx( t, Z ).';
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Zproj = project( self, t, Z )
      options = optimoptions('fmincon',...
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true, ...
        'CheckGradients', false, ...
        'Display','off' ...
      );
      % project position (x,y)
      fun     = @(ZZ) self.Qtarget( ZZ, Z );
      nonlcon = @(ZZ) self.Constraints( t, ZZ );
      LUbound = Inf*ones(size(Z));
      Zproj   = fmincon(...
        fun, Z, [], [], [], [], -LUbound, LUbound, nonlcon, options ...
      );
    end
  end
end
