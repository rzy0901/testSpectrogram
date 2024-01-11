function [freq,c,plotArgs] = parsePatternElevationInputs(UseArray,H,freq,varargin)
%This function is for internal use only. It may be removed in the future.

%PARSEPATTERNELEVATIONINPUTS parse inputs to pattern

%   Copyright 2014-2017 The MathWorks, Inc.

% check if the right arguments are passed with Array Object and Element
% system object

p = inputParser;
addRequired(p,'freq');
addOptional(p,'az',0);

addParameter(p,'Type','directivity');
addParameter(p,'PropagationSpeed',physconst('lightspeed'));
addParameter(p,'Weights',[]);
addParameter(p,'SteerAngle',[]);
addParameter(p,'Elevation',[]);
addParameter(p,'ElementWeights',[]);

parse(p,freq,varargin{:});

freq = p.Results.freq;
az = p.Results.az;
Type = p.Results.Type;
PropagationSpeed = p.Results.PropagationSpeed;
Weights = p.Results.Weights;
SteerAngle = p.Results.SteerAngle;
el = p.Results.Elevation;
ElementWeights = p.Results.ElementWeights;

if isempty(el)
    el = -90:90;
    TwoSided = true;
else
    TwoSided = false;
end

if ~isempty(Weights) && ~UseArray
    error(message('phased:system:element:NotSupportedForElement','Weights','patternElevation'));
elseif isempty(Weights)
    Weights = ones(getDOF(H),1);
end

format = 'polar';

Type = validatestring(Type,{'directivity','efield','power','powerdb'},...
    'patternElevation','Type');
if strcmp(Type,'directivity')
    unit = 'dbi';
elseif strcmp(Type,'efield')
    unit = 'mag';
elseif strcmp(Type,'power')
    unit = 'pow';
else
    unit = 'db';
end

if isPolarizationCapable(H)
    polarization = 'combined';
else
    polarization = 'None';
end

sigdatatypes.validateFrequency(freq,'patternElevation','FREQ',{'scalar','positive'});
sigdatatypes.validateAngle(el,'patternElevation','ELEVATION',{'row','nondecreasing','>=',-90,'<=',90});
sigdatatypes.validateSpeed(PropagationSpeed,'patternElevation','PropagationSpeed',{'scalar'});
sigdatatypes.validateAngle(az,'patternElevation','AZ',{'row','nondecreasing','>=',-180,'<=',180});

respcut = 'el';
cutangle = az;
azimuthangles = [];
elevationangles = el;
ugrid = [];
vgrid = [];

overlayfreq = true;
normalizeresponse = false;

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
        sigdatatypes.validateAzElAngle(SteerAngle,'patternElevation','SteerAngle',{'size',[2 1]});
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
    'FuncName','patternElevation','TwoSidedElevation',TwoSided};

