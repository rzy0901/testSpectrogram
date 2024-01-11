function [freq,c,plotArgs] = parsePatternInputs(UseArray,H,freq,varargin)
%This function is for internal use only. It may be removed in the future.

%PARSEPATTERNINPUTS parse inputs to pattern

%   Copyright 2014-2017 The MathWorks, Inc.

% check if the right arguments are passed with Array Object and Element
% system object

p = inputParser;
addRequired(p,'freq');
addOptional(p,'az',[]);
addOptional(p,'el',[]);

addParameter(p,'CoordinateSystem','polar');
addParameter(p,'Type','directivity');
addParameter(p,'Normalize',[]);
addParameter(p,'PlotStyle','Overlay');
addParameter(p,'Polarization','');
addParameter(p,'PropagationSpeed',physconst('lightspeed'));
addParameter(p,'Weights',[]);
addParameter(p,'SteerAngle',[]);
addParameter(p,'ElementWeights',[]);

parse(p,freq,varargin{:});

freq = p.Results.freq;
az = p.Results.az;
el = p.Results.el;
CoordinateSystem = p.Results.CoordinateSystem;
Type = p.Results.Type;
Normalize = p.Results.Normalize;
PlotStyle = p.Results.PlotStyle;
Polarization = p.Results.Polarization;
PropagationSpeed = p.Results.PropagationSpeed;
Weights = p.Results.Weights;
SteerAngle = p.Results.SteerAngle;
ElementWeights = p.Results.ElementWeights;

if ~isempty(Weights) && ~UseArray
    error(message('phased:system:element:NotSupportedForElement','Weights','pattern'));
elseif isempty(Weights)
    Weights = ones(getDOF(H),1);
end

CoordinateSystem = validatestring(CoordinateSystem,...
    {'polar','rectangular','uv'},'pattern','CoordinateSystem');
if strcmp(CoordinateSystem,'rectangular')
    format = 'line';
else % polar and uv
    format = CoordinateSystem;
end

Type = validatestring(Type,{'directivity','efield','power','powerdb'},...
    'pattern','Type');
if strcmp(Type,'directivity')
    unit = 'dbi';
elseif strcmp(Type,'efield')
    unit = 'mag';
elseif strcmp(Type,'power')
    unit = 'pow';
else
    unit = 'db';
end

if ~isempty(Polarization) && ~isPolarizationCapable(H)
    error(message('phased:polarization:invalidElementPolarizationSetting'));
elseif isPolarizationCapable(H)
    if ~isempty(Polarization) && strcmp(Type,'directivity')
        error(message('phased:system:array:IrrelevantSetting','Polarization','Type','directivity'));
    end
    if isempty(Polarization)
        Polarization = 'combined';
    end
    Polarization = validatestring(Polarization,{'combined','H','V'},...
        'pattern','Polarization');
    polarization = Polarization;
else
    polarization = '';
end

if isempty(polarization)
    polarization = 'None';
end
       
sigdatatypes.validateFrequency(freq,'pattern','FREQ',{'row','positive'});
sigdatatypes.validateSpeed(PropagationSpeed,'pattern','PropagationSpeed',{'scalar'});
if ~strcmp(CoordinateSystem,'uv')
    ugrid = [];
    vgrid = [];
    if isempty(az) && isempty(el)
        if isElementFromAntenna(H)
            az = -180:5:180;
            el = -90:5:90;
        else
            az= -180:180;
            el = -90:90;
        end
    elseif isempty(az)
        if (isempty(el) || ~isscalar(el)) && isElementFromAntenna(H)
            az = -180:5:180;
        else
            az = -180:180;
        end
    elseif isempty(el)
        if (isempty(az) || ~isscalar(az)) && isElementFromAntenna(H)
            el = -90:5:90;
        else
            el = -90:90;
        end
    end
    sigdatatypes.validateAngle(az,'pattern','AZ',...
        {'nondecreasing','>=',-180,'<=',180});
    sigdatatypes.validateAngle(el,'pattern','EL',...
        {'nondecreasing','>=',-90,'<=',90});
    if isscalar(freq)
        if isscalar(el) 
            respcut = 'az';
            cutangle = el;
            azimuthangles = az;
            elevationangles = [];
        elseif isscalar(az) 
            respcut = 'el';
            cutangle = az;
            azimuthangles = [];
            elevationangles = el;
        else
            respcut = '3d';
            cutangle = [];
            azimuthangles = az;
            elevationangles = el;
        end
    else
        if ~isscalar(az) && ~isscalar(el)
            error(message('phased:system:array:AtLeastOneScalarOutOfTwo','AZ','EL'));
        elseif isscalar(az)
            respcut = 'el';
            cutangle = az;
            azimuthangles = [];
            elevationangles = el;
        else %isscalar(el)
            respcut = 'az';
            cutangle = el;
            azimuthangles = az;
            elevationangles = [];
        end
    end
    
else  % uv
    if isempty(az)
        az = -1:0.01:1;
    end
    if isempty(el)
        el = -1:0.01:1;
    end
    u = az;
    v = el;
    azimuthangles = [];
    elevationangles = [];
    cutangle = [];  % not applicable
    validateattributes(u,{'numeric'},{'real','nondecreasing','>=',-1,'<=',1},...
        'pattern','U');
    validateattributes(v,{'numeric'},{'real','nondecreasing','>=',-1,'<=',1},...
        'pattern','V');
    if isscalar(freq)
        if ~isscalar(u) && ~isscalar(v)
            respcut = '3d';
            ugrid = u;
            vgrid = v;
        else
            if ~isscalar(v) || v~=0
                error(message('phased:system:array:ValueMustBeZero','V'));
            end
            respcut = 'u';
            cutangle = v;
            ugrid = u;
            vgrid = [];
        end
        
    else
        if ~isscalar(u) && ~isscalar(v)
            error(message('phased:system:array:AtLeastOneScalarOutOfTwo','U','V'));
        elseif ~isscalar(v) || v~=0
            error(message('phased:system:array:ValueMustBeZero','V'));
        end
        respcut = 'u';
        ugrid = u;
        vgrid = [];
    end
    
end

PlotStyle = validatestring(PlotStyle,{'Overlay','Waterfall'},'pattern','PlotStyle');
if strcmp(PlotStyle,'Waterfall')
    if isscalar(freq)
        error(message('phased:system:array:expectRow','FREQ'));
    end
    if strcmp(CoordinateSystem,'polar')
        error(message('phased:system:array:IrrelevantSetting','waterfall','CoordinateSystem','polar'));
    end
    overlayfreq = false;
else 
    overlayfreq = true;
end

if ~isempty(Normalize) && strcmp(Type,'directivity')
    error(message('phased:system:array:IrrelevantSetting','Normalize','Type','directivity'));
elseif ~strcmp(Type,'directivity')
    if isempty(Normalize)
        Normalize = true;
    end
    validateattributes(Normalize,{'logical'},{'scalar'},'pattern','Normalize');
    normalizeresponse = Normalize;
else
    normalizeresponse = false;
end


if UseArray && isa(H,'phased.internal.AbstractSubarray') && ...
        ~strncmp(H.SubarraySteering,'None',1)
    SteerSubarray = true;
else
    SteerSubarray = false;
end

if SteerSubarray
    if ~strncmp(H.SubarraySteering,'Custom',1)
        if isempty(SteerAngle)
            SteerAngle = [0;0];
        end
        if isscalar(SteerAngle)
            SteerAngle = [SteerAngle;0];
        end
        sigdatatypes.validateAzElAngle(SteerAngle,'pattern','SteerAngle',{'size',[2 1]});
        steerangle = SteerAngle;
        elementweights = [];
    else
        % does not support multiple weights for multiple
        % frequency yet because this is still an analog
        % behavior so at any moment, there is only one set of
        % weights can be applied.
        Ns = getNumSubarrays(H);
        Nse = zeros(1,Ns);
        for m = 1:Ns
            Nse(m) = getNumElements(H,m);
        end
        if isempty(ElementWeights)
            for m = Ns:-1:1
                ws{m} = ones(Nse(m),1);
            end
        else
            ws = ElementWeights;
        end

        if (~iscell(ws) && ~ismatrix(ws)) || isempty(ws)
            error(message('phased:phased:expectedCellOrMatrix','ElementWeights'));
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
        elementweights = ws;
        steerangle = [];
    end
else
    if ~isempty(SteerAngle)
        error(message('phased:system:array:InvalidSteeringSetting',...
            'SubarraySteering','Phase','Time'));
    else
        steerangle = SteerAngle;
    end
    if ~isempty(ElementWeights)
        error(message('phased:system:array:InvalidCustomSteeringSetting',...
            'SubarraySteering','Custom'));
    else
        elementweights = ElementWeights;
    end
end

c = PropagationSpeed;
plotArgs = {'Format',format,'RespCut',respcut,'CutAngle',cutangle,...
    'Unit',unit,'NormalizeResponse',normalizeresponse,...
    'OverlayFreq',overlayfreq,'SteerAngle',steerangle,...
    'Polarization',polarization,'AzimuthAngles',azimuthangles,...
    'ElevationAngles',elevationangles,'UGrid',ugrid,'VGrid',vgrid,...
    'Weights',Weights,'ElementWeights',elementweights,...
    'FuncName','pattern','TwoSidedElevation',false};

