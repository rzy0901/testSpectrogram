classdef (Sealed, StrictDefaults) RangeAngle < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%RangeAnge Range and Angle Calculator
%   H = phased.internal.RangeAngle creates a System object, H, which
%   calculates the range and/or azimuth and elevation angles of several 
%   positions with respect to a reference.
%
%   H = phased.internal.RangeAngle(Name,Value) creates a range and angle
%   calculator System object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair 
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   ANG = step(H,POS) returns the azimuth and elevation angles, ANG 
%   (in degrees), of the input positions POS with respect to the reference  
%   position specified in the ReferencePosition property.
%
%   RNG = step(H,POS) returns the range, RNG (in meters), when you set
%   the Output property to 'Range'.
%
%   POS must be a 3xN matrix with each column specifying a position in the
%   form of [x; y; z] (in meters) coordinates. RNG is a 1xN vector whose
%   entries are the ranges for the corresponding positions specified in
%   POS. ANG is an 2xN matrix whose columns are the angles, in the form of
%   [azimuth; elevation], for the corresponding positions specified in POS.
%
%   [RNG ANG] = step(H,POS) returns the range and angles when you set
%   the Output property to 'Range and Angle'.
%
%   [...] = step(H,POS,REFPOS) uses REFPOS as the reference position
%   when you set the ReferencePositionSource property to 'Input port'.
%
%   [...] = step(H,POS,REFAXES) uses REFAXES as the reference axes when you
%   set the ReferenceAxesSource property to 'Input port' and the
%   ReferencePositionSource property to 'Property'.
%
%   [...] = step(H,POS,REFPOS,REFAXES) uses REFPOS as the reference
%   position and REFAXES as the reference axes when you set the
%   ReferenceAxesSource property to 'Input port' and the
%   ReferencePositionSource property to 'Input port'.
%
%   % Example:
%   %   A target is located at [500; 0; 0] meters and a receiver is located
%   %   at [100; 100; 100] meters.  Determine the range and angles of the 
%   %   target with respect to the receiver.
%
%   rngang = phased.internal.RangeAngle('Output','Range and Angle', ...
%                                  'ReferencePositionSource','Input port');
%   [tgt_rng,tgt_ang] = step(rngang,[500; 0; 0],[100; 100; 100]);

%   Copyright 2014-2015 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
properties (Nontunable)
    %Model  Propagation model 
    %   Specify the propagation model used for range and angle computation
    %   as one of 'Free space' | 'Two-ray', where the default is 'Free
    %   space'.
    Model = 'Free space'
    %ReferencePositionSource Reference position source
    %   Specify the source of reference position using one of
    %   'Property' | 'Input port', where the default is 'Property'.
    ReferencePositionSource = 'Property';  
    %ReferencePosition Reference position
    %   Specify position to calculate the range or angles from in the 
    %   form of [x; y; z] (in meters) coordinates. The default value 
    %   is [0;0;0]. This property applies when you set the
    %   ReferencePositionSource property to 'Property'. 
    ReferencePosition = [0; 0; 0];
    %ReferenceAxesSource Reference axes source
    %   Specify the source of reference orientation axes using one of
    %   'Property' | 'Input port', where the default is 'Property'.
    ReferenceAxesSource = 'Property';  
    %ReferenceAxes Reference axes
    %   Specify the reference axes to calculate the range or angles from in
    %   the form of a 3x3 matrix. Each column of the matrix specifies the
    %   direction of an axis for the coordinate system in the form of [x;
    %   y; z] coordinates. The default value is the identity matrix. This
    %   property applies when you set the ReferenceAxesSource property
    %   to 'Property'.
    ReferenceAxes = eye(3);
    %Output  Output(s)
    %   Specify the output ports to enable as one of 'Angle' | 'Range'
    %    | 'Range and Angle', where the default is 'Angle'. 
    Output = 'Angle'
end

properties (Constant, Hidden)
    ReferencePositionSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    ReferenceAxesSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    OutputSet = matlab.system.StringSet({'Angle','Range','Range and Angle'});
    ModelSet = matlab.system.StringSet({'Free space','Two-ray'});
end

properties (Access = private, Nontunable)
    % output 
    pOutput
end

properties (Access = private, Nontunable, Logical) 
    %pReferencePositionViaProp - private flag whether the reference is
    %specified via property
    pReferencePositionViaProp;
    %pReferenceAxesViaProp - private flag whether the reference is
    %specified via property
    pReferenceAxesViaProp;
    
end


methods
    function obj = RangeAngle(varargin)
        setProperties(obj, nargin, varargin{:});
    end

    function set.ReferencePosition(obj,value)
        sigdatatypes.validate3DCartCoord(value,'','ReferencePosition',{'column'});
        obj.ReferencePosition = value;
    end
    
    function set.ReferenceAxes(obj,value)
        sigdatatypes.validate3DCartCoord(value,'','ReferenceAxes',{'size',[3 3]});
        obj.ReferenceAxes = value;
    end

end
methods (Access=protected)

    function setupImpl(obj)
        obj.pReferencePositionViaProp = ...
            (obj.ReferencePositionSource(1) == 'P'); %Property
        obj.pReferenceAxesViaProp = ...
            (obj.ReferenceAxesSource(1) == 'P'); %Property
        if strcmp(obj.Output, 'Angle')
            obj.pOutput = 1;
        elseif strcmp(obj.Output, 'Range')
            obj.pOutput = 2;
        else % Range and Angle
            obj.pOutput = 3;
        end
    end
    
    function varargout = stepImpl(obj, pos, refPosArg, refAxesArg)

        if obj.pReferencePositionViaProp
            refPos = obj.ReferencePosition;
        else
            refPos = refPosArg;
        end
        if obj.pReferenceAxesViaProp
            refAxes = obj.ReferenceAxes;
        else
            if obj.pReferencePositionViaProp
                refAxes = refPosArg;
            else
                refAxes = refAxesArg;
            end
        end
        if strcmp(obj.Model,'Free space')
            lclcoord = phased.internal.global2localcoord(pos,'rs',refPos,refAxes);
        else
            cond = any(pos(3,:)*refPos(3) < 0);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:rangeangle:invalidHeight','Pos','RefPos');
            end
            pos_dr = reshape([pos;pos(1:2,:);-pos(3,:)],3,[]);
            lclcoord = phased.internal.global2localcoord(pos_dr,'rs',refPos,refAxes);
        end
        if obj.pOutput == 1; %Angle
            varargout{1} = lclcoord(1:2,:);
        elseif obj.pOutput == 2; %Range
            varargout{1} = lclcoord(3,:);
        else
            varargout{1} = lclcoord(3,:);
            varargout{2} = lclcoord(1:2,:);
        end
    end


    function num = getNumInputsImpl(obj)
        num = 1;
        if strncmpi(obj.ReferencePositionSource,'Input port',1)
            num = num + 1;            
        end
        if strncmpi(obj.ReferenceAxesSource,'Input port',1)
            num = num + 1;            
        end

    end

    function num = getNumOutputsImpl(obj)
        if strcmp(obj.Output, 'Range and Angle')
            num = 2;
        else
            num = 1;
        end
    end


    function validateInputsImpl(obj, pos, refpos, refaxes)

        sigdatatypes.validate3DCartCoord(pos,'','Pos');
        if strcmp(obj.ReferencePositionSource,'Input port')
            sigdatatypes.validate3DCartCoord(refpos,'','RefPos',{'column'});
        end
        if strcmp(obj.ReferenceAxesSource,'Input port')
            if strcmp(obj.ReferencePositionSource,'Input port')
                sigdatatypes.validate3DCartCoord(refaxes,'','RefAxes',{'size',[3 3]});
            else
                sigdatatypes.validate3DCartCoord(refpos,'','RefAxes',{'size',[3 3]});
            end
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        if isLocked(obj)
            s.pReferencePositionViaProp = obj.pReferencePositionViaProp;
            s.pReferenceAxesViaProp = obj.pReferenceAxesViaProp;
            s.pOutput = obj.pOutput;
        end
    end
    
    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if strcmp(prop,'ReferencePosition') && ...
                strcmp(obj.ReferencePositionSource,'Input port')
            flag = true;
        end
        if strcmp(prop,'ReferenceAxes') && ...
                strcmp(obj.ReferenceAxesSource,'Input port')
            flag = true;
        end
    end
end

methods (Access = protected) %for Simulink
    
    function varargout = getInputNamesImpl(obj) 
        varargout = {'Pos'};
        if strcmp(obj.ReferencePositionSource,'Input port')
            varargout{end+1} = 'RefPos'; 
        end
        if strcmp(obj.ReferenceAxesSource,'Input port')
            varargout{end+1} = 'RefAxes'; 
        end
    end

    function varargout = getOutputNamesImpl(obj)
        if strcmp(obj.Output, 'Angle')
            varargout{1} = 'Ang';
        elseif strcmp(obj.Output, 'Range')
            varargout{1} = 'Range';
        else % Range and Angle
            varargout{1} = 'Range';
            varargout{2} = 'Ang';
        end
    end
    function varargout = getOutputSizeImpl(obj)
        posSz = propagatedInputSize(obj,1);
        if strcmp(obj.Model,'Free space')
            N = posSz(2);
        else % 'Two-ray'
            N = 2*posSz(2);
        end
        if strcmp(obj.Output, 'Angle')
            varargout{1} = [2 N];
        elseif strcmp(obj.Output, 'Range')
            varargout{1} = [1 N];
        else % Range and Angle
            varargout{1} = [1 N];
            varargout{2} = [2 N];
        end
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
        if strcmp(obj.Output, 'Range and Angle')
            varargout{2} = propagatedInputFixedSize(obj, 1);
        end
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj, 1);
        if strcmp(obj.Output, 'Range and Angle')
            varargout{2} = propagatedInputDataType(obj, 1);
        end
    end
    function varargout = isOutputComplexImpl(obj)
        varargout{1} = false;
        if strcmp(obj.Output, 'Range and Angle')
            varargout{2} = false;
        end
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Range\nAngle');
    end        
end

methods (Static,Hidden,Access=protected)
    function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:RangeAngleTitle')),...
          'Text',getString(message('phased:library:block:RangeAngleDesc')));
    end
end
end
