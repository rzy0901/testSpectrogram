function varargout = rocpfa(Pfa,varargin)
%rocpfa Receiver operating characteristic curves on varying Pfa
%   [Pd, SNR] = rocpfa(Pfa) returns the calculated receiver operating
%   characteristic (ROC) curve, i.e., probability of detection (Pd) vs.
%   single pulse signal to noise ratio (SNR) (in dB), for a given
%   probability of false alarm (Pfa). Pd and SNR are returned in columns. 
%
%   If a vector of Pfa values are specified, Pd will be a matrix with each
%   column corresponding to a given Pfa value. SNR is always a column
%   vector and is shared by all Pfa values.
%
%   [Pd, SNR] = rocpfa(...,'SignalType',TYPE) specifies the type of the
%   received signal when calculating the ROC curve. TYPE can be one of the
%   following: 'Real' | 'NonfluctuatingCoherent' | 
%   'NonfluctuatingNoncoherent' | 'Swerling1' | 'Swerling2' | 'Swerling3' |
%   'Swerling4'. The default TYPE is 'NonfluctuatingCoherent'. The noise is
%   assumed to be white Gaussian.
%
%   [Pd, SNR] = rocpfa(...,'MinSNR',MINSNR) specifies the minimum SNR (in
%   dB) included in the ROC calculation. The default value of MINSNR is 0.
%
%   [Pd, SNR] = rocpfa(...,'MaxSNR',MAXSNR) specifies the maximum SNR (in
%   dB) included in the ROC calculation. The default value of MAXSNR is 20.
%
%   [Pd, SNR] = rocpfa(...,'NumPoints',NUMPTS) specifies the number of
%   points, i.e., steps, used to calculate the ROC curves. The default
%   value of NUMPTS is 101.
%
%   [Pd, SNR] = rocpfa(...,'NumPulses',NUMINT) integrates NUMINT pulses in
%   calculating the ROC curves. The default value of NUMINT is 1, i.e., no
%   pulse integration.
%
%   rocpfa(...) plots the receiver operating characteristic curves.
%
%   % Example:
%   %   Plot the ROC curve for the nonfluctuating-coherent signal type
%   %   using a set of given Pfa values.
%
%   p = [1e-10 1e-8 1e-6 1e-4]; % Pfa 
%   rocpfa(p,'SignalType','NonfluctuatingCoherent');
%
%   See also phased, rocsnr, npwgnthresh.

%   Copyright 2010-2018 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing

%#codegen
%#ok<*EMCA>

cond = nargin > 0;
if ~cond
    coder.internal.assert(cond,'MATLAB:narginchk:notEnoughInputs');
end
cond = nargout < 3;
if ~cond
    coder.internal.assert(cond,'MATLAB:nargoutchk:tooManyOutputs');
end

[SignalType, MinSNR, MaxSNR, NumPoints, NumPulses] = ...
    parseInput(Pfa,varargin{:});

SNR = linspace(MinSNR,MaxSNR,NumPoints);
d = db2pow(SNR);

Pd = privrocpdcalc(Pfa,d,NumPulses,SignalType);
Pd = Pd.';  % match number of rows to number of SNR values

if ~nargout
    privrocplot('rocpfa',SNR,Pd,Pfa,@plot,'SNR (dB)','Pfa','',SignalType);
else
    varargout{1} = Pd;
end
if nargout > 1
    varargout{2} = SNR(:);
end

%--------------------------------------------------------------------------
function [SignalType, MinSNR, MaxSNR, NumPoints, NumPulses] = ...
    parseInput(Pfa,varargin)

eml_assert_no_varsize(1:nargin,Pfa,varargin{:});
% Define default values for optional inputs.
defaultSignalType = 'NonfluctuatingCoherent';
defaultMinSNR = 0;
defaultMaxSNR = 20;
defaultNumPoints = 101;
defaultNumPulses = 1;

 if isempty(coder.target)
     p = inputParser;
     p.addParameter('SignalType',defaultSignalType);   
     p.addParameter('MinSNR',defaultMinSNR); 
     p.addParameter('MaxSNR',defaultMaxSNR); 
     p.addParameter('NumPoints',defaultNumPoints);  
     p.addParameter('NumPulses',defaultNumPulses);  
     p.parse(varargin{:});
     SignalType = p.Results.SignalType;
     MinSNR = p.Results.MinSNR;
     MaxSNR = p.Results.MaxSNR;
     NumPoints = p.Results.NumPoints;     
     NumPulses = p.Results.NumPulses;          
 else

     parms = struct('SignalType',uint32(0), ...
           'MinSNR',uint32(0), ...
           'MaxSNR',uint32(0), ...
           'NumPoints',uint32(0), ...                    
           'NumPulses',uint32(0));

     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     SignalType = ...
         eml_get_parameter_value(pstruct.SignalType,defaultSignalType,varargin{:});
     MinSNR = eml_get_parameter_value(pstruct.MinSNR,defaultMinSNR,varargin{:});
     MaxSNR = eml_get_parameter_value(pstruct.MaxSNR,defaultMaxSNR,varargin{:});
     NumPoints = eml_get_parameter_value(pstruct.NumPoints,defaultNumPoints,varargin{:});     
     NumPulses = eml_get_parameter_value(pstruct.NumPulses,defaultNumPulses,varargin{:});          
 end

SignalType = validatestring(SignalType,{'Real','NonfluctuatingCoherent',...
    'NonfluctuatingNoncoherent','Swerling1','Swerling2','Swerling3',...
    'Swerling4'},'rocpfa','SignalType');

validateattributes(Pfa,{'numeric'},{'vector','positive','>',0,'<=',1},...
    'rocpfa','Pfa');

validateattributes(MinSNR,{'numeric'},{'finite','scalar','real'},...
    'rocpfa','MinSNR');

validateattributes(MaxSNR,{'numeric'},{'finite','scalar','real'},...
    'rocpfa','MaxSNR');

cond = MinSNR < MaxSNR;
if ~cond
    coder.internal.assert(cond,'phased:rocpfa:InvalidInput');
end
 
validateattributes(NumPoints,{'numeric'},{'scalar','positive','integer'},...
    'rocpfa','NumPoints');

validateattributes(NumPulses,{'numeric'},{'scalar','positive','integer'},...
    'rocpfa','NumPulses');



% [EOF]
