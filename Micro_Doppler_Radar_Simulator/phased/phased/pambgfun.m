function varargout = pambgfun(x,fs,varargin)
%pambgfun Periodic ambiguity function
%   PAFMAG = pambgfun(X,Fs) returns the magnitude of the normalized
%   periodic ambiguity function in PAFMAG for a single period of a periodic
%   signal specified in X. X must be a vector. Fs is the sampling frequency
%   in Hertz.
%
%   PAFMAG is an M-by-N matrix, where M is the number of Doppler shifts and
%   N is the number of time delays.
%
%   [PAFMAG,DELAY,DOPPLER] = pambgfun(X,Fs) also returns the time delay
%   vector, DELAY, and the Doppler shift vector, DOPPLER.
%
%   [PAFMAG,DELAY, Doppler] = pambgfun(X,Fs,P) returns the P-period
%   periodic ambiguity function PAFMAG. P is the number of identical
%   reference pulses.
%
%   Note that the syntax pambgfun(X,Fs) is equivalent to the syntax
%   pambgfun(X,Fs,P,'Cut','2D') with P = 1.
%
%   [PAFMAG,DELAY] = pambgfun(...,'Cut','Doppler') returns the periodic
%   ambiguity function PAFMAG along zero Doppler cut. DELAY contains the
%   corresponding time delay vector.
%
%   [PAFMAG,DOPPLER] = pambgfun(...,'Cut','Delay') returns the periodic
%   ambiguity function PAFMAG along zero delay cut. DOPPLER contains the
%   corresponding Doppler shift vector.
%
%   [...] = pambgfun(...,'Cut','Doppler','CutValue',V) or [...] =
%   pambgfun(...,'Cut','Delay','CutValue',V) returns the cut of the
%   periodic ambiguity function at the specified Delay or Doppler value, V.
%   If 'Cut'is 'Delay', V should be between -M and M where M is the
%   duration of X in seconds.  If 'Cut' is 'Doppler', V should be between
%   -Fs/2 and Fs/2. The parameter 'CutValue' is only applicable when 'Cut'
%   is set to 'Delay' or 'Doppler'. The default value of V is 0. The units
%   for V are in Hertz.
%   
%   pambgfun(...) produces a contour plot of the periodic ambiguity
%   function if 'Cut' is '2D' or a line plot of the periodic ambiguity
%   function cut if 'Cut' is 'Delay' or 'Doppler'.
%
%   % Example:
%   %   Plot the periodic ambiguity function of a rectangular pulse. Assume
%   %   the sampling rate is 10 Hz.
%
%   x = [ones(1,10) zeros(1,5)]; fs = 10;
%   pambgfun(x,fs);
%
%   See also ambgfun, phased, phased.LinearFMWaveform, 
%   phased.RectangularWaveform, phased.SteppedFMWaveform, dutycycle.

%   Copyright 2016-2018 The MathWorks, Inc.

%   Reference
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [2] Eli Brookner, Radar Technology, Lex Books, 1996
%   [3] E. Brigham, The Fast Fourier Transform and Its Applications,
%       Prentice Hall, 1988
%   [4] Levanon, N. and E. Mozeson. Radar Signals. Hoboken, NJ: John Wiley
%       & Sons, 2004.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(2,7,nargin);
[cut, cutvalue, P] = parseInput(x,fs,varargin{:});
xlen = numel(x(:));
t = (0:(xlen-1))/fs;
xnorm = x(:)/norm(x);

tau = -(xlen-1):(xlen-1);
taulen = 2*xlen - 1;

nfreq = 2^nextpow2(taulen);

Fd = -fs/2:fs/nfreq:fs/2-fs/nfreq;

if strcmp(cut,'Delay')
    cutValueIdx = val2ind(cutvalue,1/fs,-(xlen-1)/fs);
    v = localshift(xnorm(:),cutValueIdx,xlen);
    amf_delaycut = nfreq.*abs(ifftshift((ifft(conj(v).*xnorm,nfreq))));
    if P~=1
        Tr = xlen/fs;
        sinc_func = zeros(numel(Fd),1);
        sinc_func(:,1) = abs((sin(P*pi*Fd*Tr))./(P*sin(pi*Fd*Tr)));
        sinc_func(isnan(sinc_func),1) = 1 ;
        amf_delaycut = amf_delaycut.*sinc_func;
    end
elseif strcmp(cut,'Doppler')
    x_temp = toeplitz(xnorm,xnorm([1 end:-1:2]));
    x_modulated = xnorm.'.*exp(1i*2*pi*cutvalue*t);
    amf_dopplercut_temp = abs(x_modulated*conj(x_temp));
    amf_dopplercut = amf_dopplercut_temp(:,[2:end 1:end]);
    if P~=1
        Tr = xlen/fs;
        sinc_func = abs((sin(P*pi*cutvalue*Tr))./(P*sin(pi*cutvalue*Tr)));
        sinc_func(isnan(sinc_func)) = 1 ;
        sinc_func_replicated = ones(1,2*xlen-1).*sinc_func;
        amf_dopplercut = amf_dopplercut.*sinc_func_replicated;
    end
else
    amf = zeros(nfreq,taulen);
    for m = 1:taulen
            v = localshift(xnorm(:),tau(m),xlen);
            amf(:,m) = abs(ifftshift((ifft(conj(v).*xnorm,nfreq))));
    end
    amf = nfreq*amf;
    Tr = xlen/fs;
    sinc_func = zeros(numel(Fd),1);
    sinc_func(:,1) = abs((sin(P*pi*Fd*Tr))./(P*sin(pi*Fd*Tr)));
    sinc_func(isnan(sinc_func),1) = 1 ;
    sinc_func_replicated = repmat(sinc_func,1,2*xlen-1);
    amf = amf.*sinc_func_replicated;
end
tau = [-t(end:-1:2) t];

if nargout == 0
    cond = isempty(coder.target);
    if ~cond
        coder.internal.assert(cond,'phased:pambgfun:invalidCodegenOutput');
    end
    [~,tau_scale,tau_label] = engunits(max(abs(tau)));
    [~,fd_scale,fd_label] = engunits(max(abs(Fd)));
    switch cut
        case '2D'
            [~,hc] = contour(tau*tau_scale,Fd*fd_scale,amf); grid on;
            set(get(hc,'Parent'),'Tag','ambg');
            xlabel_str = getString(message('phased:pambgfun:xlabel2D',texlabel('tau'),tau_label));
            xlabel(xlabel_str);
            ylabel_str = getString(message('phased:pambgfun:ylabel2D',texlabel('f_d'),fd_label));
            ylabel(ylabel_str);
            title(getString(message('phased:pambgfun:PeriodicAmbiguityFunction')));
            colorbar;
        case 'Doppler'
            hl = plot(tau*tau_scale,amf_dopplercut); grid on;
            set(get(hl,'Parent'),'Tag','ambg');
            xlabel_str = getString(message('phased:pambgfun:xlabel2D',texlabel('tau'),tau_label));
            xlabel(xlabel_str);
            ylabel(getString(message('phased:pambgfun:PeriodicAmbiguityFunction')));
            title(getString(message('phased:pambgfun:DopplerCutTitle',num2str(cutvalue*fd_scale),fd_label)));
        case 'Delay'
            hl = plot(Fd*fd_scale,amf_delaycut); grid on;
            set(get(hl,'Parent'),'Tag','ambg');
            xlabel_str = getString(message('phased:pambgfun:ylabel2D',texlabel('f_d'),fd_label));
            xlabel(xlabel_str);
            ylabel(getString(message('phased:pambgfun:PeriodicAmbiguityFunction')));
            title(getString(message('phased:pambgfun:DelayCutTitle',num2str(cutvalue*tau_scale),tau_label)));
    end
else
    
    if strcmp(cut,'2D')
      varargout{1} = amf;
      varargout{2} = tau;
      varargout{3} = Fd;
    elseif strcmp(cut,'Doppler')
       varargout{1} = amf_dopplercut; 
       varargout{2} = tau;
    else
       varargout{1} = amf_delaycut; 
       varargout{2} = Fd;
    end
end

%--------------------------------------
function [cut, cutvalue,P] = parseInput(x,fs,varargin)
eml_assert_no_varsize(1:nargin, x,fs,varargin{:});
validateattributes(x,{'double'},{'vector'},'pambgfun','X');
validateattributes(fs,{'double'},{'scalar','positive','finite'},'pambgfun','Fs');
if rem(nargin,2)~=0
    P = varargin{1};
    validateattributes(P,{'double'},{'scalar','positive','finite'},'pambgfun','N');
    numericInput = 2;
else
    P = 1; 
    numericInput = 1;
end


defaultCut = '2D';
defaultCutValue = 0;

if isempty(coder.target)
     p = inputParser;
     p.addParameter('Cut',defaultCut);
     p.addParameter('CutValue',defaultCutValue);
     p.parse(varargin{numericInput:end});
     cut = p.Results.Cut;
     cutvalue = p.Results.CutValue;
 else
     parms = struct('Cut',uint32(0), ...
                    'CutValue',uint32(0));
     pstruct = eml_parse_parameter_inputs(parms,[],varargin{numericInput:end});
     cut = eml_get_parameter_value(pstruct.Cut,defaultCut,varargin{numericInput:end});
     cutvalue = eml_get_parameter_value(pstruct.CutValue,defaultCutValue,varargin{numericInput:end});
 end

cut =  validatestring(cut,{'Delay','Doppler','2D'},'pambgfun','CUT');
xlen = numel(x);
if strcmp(cut, 'Delay')
    validateattributes(cutvalue,{'double'},{'scalar','finite','real','>=',-(xlen-1)/fs,'<=',(xlen-1)/fs},'pambgfun','CUTVALUE');
else
    validateattributes(cutvalue,{'double'},{'scalar','finite','real','>=',-fs/2,'<=',fs/2},'pambgfun','CUTVALUE');
end

function v = localshift(x,tau,xlen)
if tau >= 0
    v = circshift(x,xlen-tau);
else
    v = circshift(x,abs(tau));
end


% [EOF]
