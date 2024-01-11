function [freq,v,plotArgs] = parsePlotResponseInputs(UseArray,H,varargin)
%This function is for internal use only. It may be removed in the future.

%PARSEPLOTRESPONSEINPUTS parse inputs to plotResponse

%   Copyright 2014-2018 The MathWorks, Inc.

freq = varargin{1};

if UseArray
    v = varargin{2};
    if ischar(v) || isstring(v)
        error(message('phased:system:array:InvalidArgumentForObject'));
    end
    sigdatatypes.validateSpeed(v,'plotResponse',...
        'V',{'scalar','positive'});
    N = getDOF(H);
    weights = ones(N,1);
    pvpair = varargin(3:end);
else
    v = physconst('lightspeed');
    pvpair = varargin(2:end);
    if ~isempty(pvpair)&&~ischar(pvpair{1})
        error(message('phased:system:element:InvalidArgumentForElement'));
    end
    weights = 1;
end


format = 'line';
respcut = '';
cutangle = [];
unit = 'db';
normalizeresponse = [];
overlayfreq = [];
steerangle = [];
polarization = '';
azimuthangles = [];
elevationangles = [];
ugrid = [];
vgrid = [];
elementweights = [];

sigutils.pvparse(pvpair{:});

format = validatestring(format,{'line','polar','uv'},...
    'plotResponse','Format');

if ~strncmp(format,'uv',2)
    if isempty(respcut)
        respcut = 'az';
    end
    respcut = validatestring(respcut,{'Az','El','3D'},...
        'plotResponse','RespCut');

else
    if isempty(respcut)
        respcut = 'u';
    end
    respcut = validatestring(respcut,{'U','3D'},...
        'plotResponse','RespCut');
   
end

unit = validatestring(unit,{'mag','pow','db','dbi'},...
    'plotResponse','Unit');
directivityflag = strcmp(unit,'dbi');

if ~directivityflag
    if isempty(normalizeresponse)
        normalizeresponse = true;
    end
    validateattributes(normalizeresponse,{'logical'},{'scalar'},...
        'plotResponse','NormalizeResponse');
else
    if ~isempty(normalizeresponse)
        error(message('phased:system:array:IrrelevantSetting','NormalizeResponse','Unit','dbi'));
    end
    normalizeresponse = false;  % directivity is not normalized
end


%check polarization capability
if ~directivityflag
    if isempty(polarization)
        polarization = 'none';
    end
    polarization = validatestring(polarization,{'None','Combined','H','V'},...
        'plotResponse','Polarization');
    %polarization can only be specified with polarization capable elements
    if ~strcmp(polarization,'None') && ~isPolarizationCapable(H)
        error(message('phased:polarization:invalidElementPolarizationSetting'));
    end
else
    if ~isempty(polarization)
        error(message('phased:system:array:IrrelevantSetting','Polarization','Unit','dbi'));
    end
end

%ensure FREQ is a scalar for 3D and row vector for 2D plots and ensure that
%WEIGHTS is a vector for 3D plots
if ~strncmpi(format,'uv',2)
    if(~isempty(ugrid))
        error(message('phased:system:array:IrrelevantSettingWhenNot','UGrid','Format','UV'));
    end
    if(~isempty(vgrid))
        error(message('phased:system:array:IrrelevantSettingWhenNot','VGrid','Format','UV'));
    end
    if strcmpi(respcut,'3D')
        sigdatatypes.validateFrequency(freq,'plotResponse',...
            'FREQ',{'scalar'});
        
    else
        sigdatatypes.validateFrequency(freq,'plotResponse',...
            'FREQ',{'row'});
        if (strncmpi(respcut,'el',2)&&(~isempty(azimuthangles)))
            error(message('phased:system:array:IrrelevantSetting','AzimuthAngles','RespCut','El'));
        elseif (strncmpi(respcut,'az',2)&&(~isempty(elevationangles)))
            error(message('phased:system:array:IrrelevantSetting','ElevationAngles','RespCut','Az'));
        end
        
        if isempty(cutangle)
            cutangle = 0;
        end
        if strncmpi(respcut,'e',1)
            sigdatatypes.validateAngle(cutangle,'plotResponse',...
                'CutAngle',{'scalar','<=',180,'>=',-180});
        else % strncmpi(respcut,'a',1)
            sigdatatypes.validateAngle(cutangle,'plotResponse',...
                'CutAngle',{'scalar','<=',90,'>=',-90});
        end
    end
    
    if isempty(azimuthangles)
        azimuthangles = -180:180;
    end
    
    if isempty(elevationangles)
        elevationangles = -90:90;
    end

    validateattributes(azimuthangles,{'numeric'},{'nondecreasing','>=',-180,'<=',180},...
        'plotResponse','AzimuthAngles');
    validateattributes(elevationangles,{'numeric'},{'nondecreasing','>=',-90,'<=',90},...
        'plotResponse','ElevationAngles');
  
    
   
else
    if ~isempty(cutangle)
        error(message('phased:system:array:IrrelevantSetting','CutAngle','Format','UV'));
    end
    if ~isempty(elevationangles)
        error(message('phased:system:array:IrrelevantSetting','ElevationAngles','Format','UV'));
    end
    if ~isempty(azimuthangles)
        error(message('phased:system:array:IrrelevantSetting','AzimuthAngles','Format','UV'));
    end
    
    if strcmpi(respcut,'3D')
        sigdatatypes.validateFrequency(freq,'plotResponse',...
            'FREQ',{'scalar'});
    else
        if (strncmpi(respcut,'u',1)&&(~isempty(vgrid)))
            error(message('phased:system:array:IrrelevantSetting','VGrid','RespCut','U'));
        end
        sigdatatypes.validateFrequency(freq,'plotResponse',...
            'FREQ',{'row'});
    end
    
    if(isempty(ugrid))
        ugrid = -1:0.01:1;
    end
     if(isempty(vgrid))
        vgrid = -1:0.01:1;
     end
      validateattributes(ugrid,{'numeric'},{'nondecreasing','>=',-1,'<=',1},...
        'plotResponse','UGrid');
      validateattributes(vgrid,{'numeric'},{'nondecreasing','>=',-1,'<=',1},...
        'plotResponse','VGrid');
end

if strcmpi(respcut,'3D')
    if ~isempty(overlayfreq)
        error(message('phased:system:array:IrrelevantSetting','OverlayFreq','RespCut','3D'));
    end
    if ~isempty(cutangle)
        error(message('phased:system:array:IrrelevantSetting','CutAngle','RespCut','3D'));
    end
end

if strncmpi(format,'polar',1)
    if ~isempty(overlayfreq)
        error(message('phased:system:array:IrrelevantSetting','OverlayFreq','Format','''Polar'''));
    end
end

if isempty(overlayfreq)
    overlayfreq = true;
end
validateattributes(overlayfreq,{'logical'},{'scalar'},...
    'plotResponse','OverlayFreq');
    
if ~overlayfreq
    if isscalar(freq)
        error(message('phased:system:array:expectRow','FREQ'));
    end
end

if UseArray && isa(H,'phased.internal.AbstractSubarray') && ...
        ~strncmp(H.SubarraySteering,'None',1)
    SteerSubarray = true;
else
    SteerSubarray = false;
end

if SteerSubarray
    if ~strncmp(H.SubarraySteering,'Custom',1)
        if isempty(steerangle)
            steerangle = [0;0];
        end
        if isscalar(steerangle)
            steerangle = [steerangle;0];
        end
        sigdatatypes.validateAzElAngle(steerangle,'plotResponse','SteerAngle',{'size',[2 1]});
    else % Custom
        % does not support multiple weights for multiple
        % frequency yet because this is still an analog
        % behavior so at any moment, there is only one set of
        % weights can be applied.

        ws = elementweights;  % weights
        Ns = getNumSubarrays(H);
        if (~iscell(ws) && ~ismatrix(ws)) || isempty(ws)
            error(message('phased:phased:expectedCellOrMatrix','ElementWeights'));
        end
        Nse = zeros(1,Ns);
        for m = 1:Ns
            Nse(m) = getNumElements(H,m);
        end
        if iscell(ws)
            if ~isrow(ws) || (numel(ws)~= Ns)
                error(message('phased:phased:expectedMatrixSize','ElementWeights',1,Ns));
            end
            for m = 1:Ns
                if ~iscolumn(ws{m}) || (numel(ws{m})~=Nse(m))
                    error(message('phased:system:array:SubarrayElementWeightsSizeMismatch',...
                        m,'ElementWeights',Nse(m)));
                end
                if ~isa(ws{m},'double')
                    error(message('phased:system:array:SubarrayElementWeightsInvalidDataType',...
                        m,'ElementWeights','double'));
                end
            end
        else
            sz_ws = size(ws);
            Nsemax = max(Nse);
            if ~isequal(sz_ws,[Nsemax Ns])
                error(message('phased:phased:expectedMatrixSize','ElementWeights',Nsemax,Ns));
            end
            if ~isa(ws,'double')
                error(message('MATLAB:system:invalidInputDataType','ElementWeights','double'));
            end
        end
    end
else
    if ~isempty(steerangle)
        error(message('phased:system:array:InvalidSteeringSetting',...
            'SubarraySteering','Phase','Time'));
    end
end

plotArgs = {'Format',format,'RespCut',respcut,'CutAngle',cutangle,...
    'Unit',unit,'NormalizeResponse',normalizeresponse,...
    'OverlayFreq',overlayfreq,'SteerAngle',steerangle,...
    'Polarization',polarization,'AzimuthAngles',azimuthangles,...
    'ElevationAngles',elevationangles,'UGrid',ugrid,'VGrid',vgrid,...
    'Weights',weights,'ElementWeights',elementweights,...
    'FuncName','plotResponse'};
