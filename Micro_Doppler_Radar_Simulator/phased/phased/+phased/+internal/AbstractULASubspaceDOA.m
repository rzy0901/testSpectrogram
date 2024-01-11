classdef (Hidden) AbstractULASubspaceDOA < phased.internal.AbstractULADOA & ...
     matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%ABSTRACTULASUBSPACEDOA Summary of this class goes here
%   Detailed explanation goes here

%   Copyright 2010-2018 The MathWorks, Inc.
%    


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %NumSignalsSource Source of number of signals
    %   Specify the source of the number of signals as one of 'Auto' | 
    %   'Property'. The default is 'Auto'. If you set this property to
    %   'Auto', the number of signals is estimated by the method specified
    %   by the NumSignalsMethod property.
    NumSignalsSource = 'Auto';
    %NumSignalsMethod Method to estimate number of signals
    %   Specify the method to estimate the number of signals as one of 
    %   'AIC' | 'MDL'. The default is 'AIC'. The 'AIC' uses the Akaike
    %   Information Criterion and the 'MDL' uses Minimum Description Length
    %   Criterion. This property applies when you set the NumSignalsSource
    %   property to 'Auto'.
    NumSignalsMethod = 'AIC';
end

properties (Nontunable, PositiveInteger) 
    %NumSignals Number of signals
    %   Specify the number of signals as a positive integer scalar. The
    %   default value is 1. This property applies when you set the
    %   NumSignalsSource property to 'Property'.
    NumSignals = 1;
end

properties(Constant, Hidden)
    NumSignalsSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    NumSignalsMethodSet = matlab.system.StringSet({ 'AIC','MDL' });
end

methods (Access = protected)
    
    function obj = AbstractULASubspaceDOA(varargin)
       obj = obj@phased.internal.AbstractULADOA(varargin{:});
    end
    
    function flag = isOutputComplexityLockedImpl(obj,index) %#ok<MANU>
        flag = false;
        if index == 1
            flag = true;
        end
    end
    
    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractULADOA(obj);
        if strcmp(obj.NumSignalsSource,'Property')
            maxNumSig = getMaxNumSignal(obj);
            cond = obj.NumSignals > maxNumSig;
            if cond
                coder.internal.errorIf(cond,'phased:phased:internal:AbstractULASubspaceDOA:InvalidNumSources', maxNumSig);
            end
        end
    end

    function maxNumSig = getMaxNumSignal(obj)
    % NumSignals should be no greater than possible signal subspace dim - 1
        maxNumSig = getEffectiveChannel(obj) - 1;
    end

    function D = getNumSignals(obj,eigenvals,K,fb)
    % Obtain signal subspace dimension
        if strcmp(obj.NumSignalsSource,'Auto')
            if strcmp(obj.NumSignalsMethod,'AIC')
                D = phased.internal.aictest(eigenvals,K,fb);
            else
                D = phased.internal.mdltest(eigenvals,K,fb);
            end
        else
            D = obj.NumSignals;
        end        
    end
    
    function [eigenvals, eigenvects] = privEig(obj,Sx)  %#ok<INUSL>
        % Calculate eigenvectors and eigenvalues assuming Sx is a positive
        % semidefinite matrix
        [eigenvectsC,eigenvalsC] = svd((Sx+Sx')/2);
        eigenvals = diag(eigenvalsC);
        eigenvals(eigenvals<0|eigenvals<=eps(eigenvals(1))) = 0;
        eigenvects = complex(eigenvectsC);
        
    end

    function doasRow = privConvertAngles(obj,doas)
    %privConvertAngles(Hdoa,DOAs) Convert angles for ULA.
    %    ANG = privConvertAngles(Hdoa,DOAs) converts estimated angle in
    %    sin-space to degrees. This method is valid for ULA only.

        wavelength = obj.PropagationSpeed/obj.OperatingFrequency;
        d_lambda = obj.SensorArray.ElementSpacing/wavelength;
        d_lambda=cast(d_lambda,class(doas));
        u = doas/(2*pi*d_lambda);
        % check whether all elements of u are within [-1,1] 
        idx = find(abs(u)<=1);
        %assert(size(idx,1) <= size(doas,1));
        if (obj.NumSignalsSource(1) == 'P') %Property
            numValidSignals = length(idx);
            numSignals = obj.NumSignals;
            if numValidSignals < numSignals && isempty(coder.target)
                warning(message('phased:phased:internal:AbstractULASubspaceDOA:InvalidPsi', obj.NumSignals));
            end
            if ~isempty(idx)
                doas = asin(u(idx));
                % In degrees
                doas = doas*180/pi;
                % convert to row vector
                doasCol = doas(:);
                doasRow = NaN(1,numSignals,class(doas));%invalid angles
                doasRow(1:numValidSignals) = doasCol.';
            else
                doasRow = NaN(1,numSignals,class(doas));%invalid angles
            end
        else %auto mode, varsize (never from Simulink)
            if isempty(idx)
                doasRow = zeros(1,0,class(doas));
            else
                doas = asin(u(idx));
                % In degrees
                doas = doas*180/pi;
                % convert to row vector
                doasCol = doas(:);
                doasRow = doasCol.';
            end
        end
    end

    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractULADOA(obj, prop);
        SourceAuto = strcmp(obj.NumSignalsSource,'Auto');
        if SourceAuto && strcmp(prop,'NumSignals')
            flag = true;
        end
        if ~SourceAuto && strcmp(prop,'NumSignalsMethod')
            flag = true;
        end
    end
    
end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl(sensorType)
      if nargin == 1
          groups = getPropertyGroupsImpl@phased.internal.AbstractULADOA(sensorType);
      else
          groups = getPropertyGroupsImpl@phased.internal.AbstractULADOA;
      end
    pNumSignalsSource = matlab.system.display.internal.Property('NumSignalsSource', ...
        'IsGraphical', false, ...
        'UseClassDefault', false,'Default','Property');

    props = {...
      pNumSignalsSource ,...
      'NumSignalsMethod',...
      'NumSignals'};
    groups(1).PropertyList = [groups(1).PropertyList props];
  end
end

methods (Access = protected) %for Simulink
    function varargout = getOutputNamesImpl(~)
        varargout = {'Ang'};
    end
    function varargout = getInputNamesImpl(~)
        varargout = {''};
    end
    function flag = isInputSizeLockedImpl(~,~)
        flag = false;
    end
    function varargout = getOutputSizeImpl(obj)
        varargout{1} = [1,obj.NumSignals];
    end
    function varargout = isOutputFixedSizeImpl(~)
        varargout = {true};
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout = {propagatedInputDataType(obj,1)};
    end
    function varargout = isOutputComplexImpl(~)
        varargout = {false};
    end
end
end

