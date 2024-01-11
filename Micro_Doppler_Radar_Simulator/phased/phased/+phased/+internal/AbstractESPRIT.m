classdef (Hidden) AbstractESPRIT < phased.internal.AbstractULASubspaceDOA
%This class is for internal use only. It may be removed in the future.

%AbstractESPRIT Abstract class for ESPRIT DOA estimation

%   Copyright 2009-2018 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Dependent, Nontunable)
    %SpatialSmoothing Spatial smoothing
    %   Specify the effective reduction in the size of the sensor array due
    %   to spatial smoothing as a nonnegative integer. If the array
    %   consists of M elements and the value of this property is L,
    %   maximally overlapped subarrays of M-L elements are formed.
    %
    %   Covariance matrices are estimated for each subarray of M-L elements
    %   and averaged together to produce a covariance matrix of size
    %   (M-L)x(M-L). The maximum value of L is M-2 resulting in subarrays
    %   consisting of two elements.
    %   
    %   Note that each additional increment in this property (decorrelates)
    %   handles one additional coherent source, but reduces the effective
    %   size of the array aperture. The default value of this property is 0
    %   resulting in no spatial smoothing.
    SpatialSmoothing 
end

properties (Nontunable)
    %Method Type of least square method
    %   Specify the least square method used for ESPRIT as one of 'TLS' |
    %   'LS'. 'TLS' refers to total least squares and 'LS' refers to least
    %   squares. The default is 'TLS'.
    Method = 'TLS';
end
    
properties(Constant, Hidden)
    MethodSet = matlab.system.StringSet({'TLS','LS'});
end

methods
    function set.SpatialSmoothing(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar', 'integer', 'nonnegative', 'finite' }, '', 'SpatialSmoothing');
        obj.pSpatialSmoothing = value;
    end
    
    function value = get.SpatialSmoothing(obj)
        value = obj.pSpatialSmoothing;
    end
    
end

methods (Access = protected)
    function obj = AbstractESPRIT(varargin)
        obj = obj@phased.internal.AbstractULASubspaceDOA(varargin{:});
    end
end

methods (Access = protected)
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractULASubspaceDOA(obj);
    end
end

methods (Access = protected)

    function psieig = privCoreESPRIT(obj,eigenvects,D,Js1,Js2,method) %#ok<MANU>
    % privCoreESPRIT  Core ESPRIT algorithm

        % Selecting subarray signal subspaces
        Us1 = Js1*eigenvects(:,1:D);
        Us2 = Js2*eigenvects(:,1:D);

        switch lower(method)
            case 'ls'
                % LS-ESPRIT
                psi = (Us1'*Us1)\Us1'*Us2;    % Eq. (9.121) in [1]
            case 'tls'
                % TLS-ESPRIT
                C = [Us1';Us2']*[Us1 Us2];    % Eq. (9.123) in [1]
                [U,~,~] = svd(C);             % C is 2*D x 2*D  
                V12 = U(1:D,D+1:2*D);         % D x D
                V22 = U(D+1:2*D,D+1:2*D);     % D x D
                psi = -V12/V22;               % Eq. (9.122) in [1]
        end

        psieig = eig(psi);

    end

end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
    groups = getPropertyGroupsImpl@phased.internal.AbstractULASubspaceDOA;
    props = {...
      'SpatialSmoothing',...
      'Method'};
    groups(1).PropertyList = [groups(1).PropertyList props];
    % SpatialSmoothing as dependent on private-only
    groups(1).DependOnPrivatePropertyList =  {'SpatialSmoothing'};
  end
end
   
end

