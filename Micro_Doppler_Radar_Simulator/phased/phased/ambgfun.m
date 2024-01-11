function varargout = ambgfun(x,varargin)
%ambgfun Ambiguity function and cross ambiguity function
%   AFMAG = ambgfun(X,Fs,PRF) returns the magnitude of the normalized
%   ambiguity function in AFMAG for the pulse train specified in X. X must
%   be a vector. Fs is the sampling frequency (in Hz). PRF is the pulse
%   repetition frequency (in Hz).
%
%   AFMAG is an M-by-N matrix, where M is the number of Doppler shifts and
%   N is the number of time delays.
%
%   [AFMAG,DELAY,DOPPLER] = ambgfun(X,Fs,PRF) also returns the time delay
%   vector, DELAY, and the Doppler shift vector, DOPPLER.
%
%   Note that the syntax ambgfun(X,Fs,PRF) is equivalent to the syntax
%   ambgfun(X,Fs,PRF,'Cut','2D').
%
%   [AFMAG, DELAY, DOPPLER] = ambgfun(X,Y,Fs,PRF) returns the magnitude of
%   the normalized cross ambiguity function in AFMAG for the pulse train
%   specified in X and Y. X and Y must be vectors. Fs is the sampling
%   frequency (in Hz). PRF is a vector containing the PRF for both X and
%   Y. Note that the syntax ambgfun(X,Y,Fs,PRF) is equivalent to the syntax
%   ambgfun(X,Y,Fs,PRF,'Cut','2D').
%
%   [AFMAG,DELAY] = ambgfun(...,'Cut','Doppler') returns the ambiguity
%   function AFMAG along zero Doppler cut. DELAY contains the corresponding
%   time delay vector.
%
%   [AFMAG,DOPPLER] = ambgfun(...,'Cut','Delay') returns the ambiguity
%   function AFMAG along zero delay cut. DOPPLER contains the corresponding
%   Doppler shift vector.
%
%   [...] = ambgfun(...,'Cut','Doppler','CutValue',V) or [...] =
%   ambgfun(...,'Cut','Delay','CutValue',V) returns the cut of the
%   ambiguity function at the specified Delay or Doppler value, V. If 'Cut'
%   is 'Delay', V should be between -M and M where M is the duration of X
%   in seconds.  If 'Cut' is 'Doppler', V should be between -Fs/2 and Fs/2.
%   The parameter 'CutValue' is only applicable when 'Cut' is set to
%   'Delay' or 'Doppler'. The default value of V is 0.
%   
%   ambgfun(...) produces a contour plot of the ambiguity function if 'Cut'
%   is '2D' or a line plot of the ambiguity function cut if 'Cut' is
%   'Delay' or 'Doppler'.
%
%   % Examples:
%
%   % Example 1:
%   %   Plot the ambiguity function of a rectangular pulse. Assume the
%   %   sampling rate is 10 Hz and the PRF is 1 Hz.
%
%   x = ones(1,10); fs = 10; prf = 1;
%   ambgfun(x,fs,prf);
%
%   % Example 2:
%   %   Plot the cross ambiguity function between a rectangular pulse and
%   %   its distorted version. Compare the cross ambiguity function with
%   %   the theoretical ambiguity function. The sampling rate is 10 Hz and
%   %   the PRF is 1 Hz.
%
%   x = [ones(5,1); zeros(5,1)]; fs = 10; xprf = 1;
%   y = [hamming(5); zeros(5,1)]; fs = 10; yprf = 1;
%   ambgfun(x,y,fs,[xprf yprf],'Cut','Doppler')
%   hold on
%   ambgfun(x,fs,xprf,'Cut','Doppler')
%   legend('Approximate matching signal', 'Known matching signal')
%
%   See also pambgfun, phased, phased.LinearFMWaveform, 
%   phased.RectangularWaveform, phased.SteppedFMWaveform, dutycycle.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [2] Eli Brookner, Radar Technology, Lex Books, 1996
%   [3] E. Brigham, The Fast Fourier Transform and Its Applications,
%       Prentice Hall, 1988

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,8,nargin);
[cut, cutvalue, xsub, xlen, ysub, ylen, fs, isCross] = parseInput(x,varargin{:});

xnorm = xsub/norm(xsub);
ynorm = ysub/norm(ysub);

seqlen = xlen+ylen;



tau = -(seqlen/2-1):(seqlen/2-1);
taulen = seqlen - 1;
nfreq = 2^nextpow2(taulen);
Fd = -fs/2:fs/nfreq:fs/2-fs/nfreq;


if xlen>=ylen
    ynorm = [ynorm;zeros(xlen-ylen,1)];
else
    xnorm = [xnorm;zeros(ylen-xlen,1)];    
end

xlen = numel(xnorm);
t = (0:xlen-1)/fs;



if strcmp(cut,'Delay')
    fftx = (fft(xnorm.',nfreq));
    xshift = ifft(fftx.*exp(1i*2*pi*Fd*cutvalue));
    
    ffty = (fft(ynorm.',nfreq));
    yoriginal = ifft(ffty);
    
    amf_delaycut = nfreq*abs(ifftshift(ifft(yoriginal.*conj(xshift),nfreq))).';
elseif strcmp(cut,'Doppler')
    if isempty(coder.target)
        ffty = fft(ynorm.',taulen);
        fftxfreqshifted = fft(xnorm.'.*exp(1i*2*pi*cutvalue*t),taulen);
        amf_dopplercut = abs(fftshift(ifft(ffty.*conj(fftxfreqshifted))));
    else
        ffty = fft(ynorm.',nfreq);
        fftxfreqshifted = fft(xnorm.'.*exp(1i*2*pi*cutvalue*t),nfreq);
        amf_temp = abs(fftshift(ifft(ffty.*conj(fftxfreqshifted))));
        dc_idx = nfreq/2+1;
        side_len = xlen-1;
        amf_dopplercut = amf_temp(dc_idx-side_len:dc_idx+side_len);
    end
else
    amf = zeros(nfreq,taulen);
    for m = 1:taulen
        v = localshift(xnorm(:),tau(m),xlen);
        amf(:,m) = abs(ifftshift(ifft(ynorm.*conj(v),nfreq)));
    end
    amf = nfreq*amf;

end

tau = tau./fs;

if nargout == 0
    cond = isempty(coder.target);
    if ~cond
        coder.internal.assert(cond,'phased:ambgfun:invalidCodegenOutput');
    end
    [~,tau_scale,tau_label] = engunits(max(abs(tau)));
    [~,fd_scale,fd_label] = engunits(max(abs(Fd)));
    switch cut
        case '2D'
            [~,hc] = contour(tau*tau_scale,Fd*fd_scale,amf); grid on;
            set(get(hc,'Parent'),'Tag','ambg');
            xlabel_str = getString(message('phased:ambgfun:xlabel2D',texlabel('tau'),tau_label));
            xlabel(xlabel_str);
            ylabel_str = getString(message('phased:ambgfun:ylabel2D',texlabel('f_d'),fd_label));
            ylabel(ylabel_str);
            if isCross
                title(getString(message('phased:ambgfun:CrossAmbiguityFunction')));
            else
                title(getString(message('phased:ambgfun:AmbiguityFunction')));
            end
            colorbar;
        case 'Doppler'
            hl = plot(tau*tau_scale,amf_dopplercut); grid on;
            set(get(hl,'Parent'),'Tag','ambg');
            xlabel_str = getString(message('phased:ambgfun:xlabel2D',texlabel('tau'),tau_label));
            xlabel(xlabel_str);
            if isCross
                ylabel(getString(message('phased:ambgfun:CrossAmbiguityFunction')));
                title(getString(message('phased:ambgfun:CrossDopplerCutTitle',num2str(cutvalue*fd_scale),fd_label)));
            else
                ylabel(getString(message('phased:ambgfun:AmbiguityFunction')));
                title(getString(message('phased:ambgfun:DopplerCutTitle',num2str(cutvalue*fd_scale),fd_label)));
            end
        case 'Delay'
            hl = plot(Fd*fd_scale,amf_delaycut); grid on;
            set(get(hl,'Parent'),'Tag','ambg');
            xlabel_str = getString(message('phased:ambgfun:ylabel2D',texlabel('f_d'),fd_label));
            xlabel(xlabel_str);
            if isCross
                ylabel(getString(message('phased:ambgfun:CrossAmbiguityFunction')));
                title(getString(message('phased:ambgfun:CrossDelayCutTitle',num2str(cutvalue*tau_scale),tau_label)));
            else
                ylabel(getString(message('phased:ambgfun:AmbiguityFunction')));
                title(getString(message('phased:ambgfun:DelayCutTitle',num2str(cutvalue*tau_scale),tau_label)));
            end
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
function v = localshift(x,tau,xlen)
vreal = zeros(xlen,1);
if ~isreal(x)
    v = complex(vreal);
else
    v = vreal;
end
if tau >= 0
    v(1:xlen-tau) = x(1+tau:xlen);
else
    v(1-tau:xlen) = x(1:xlen+tau);
end

%--------------------------------------
function [cut, cutvalue, xsub, xlen, ysub, ylen, fs, isCross] = parseInput(x,varargin)
eml_assert_no_varsize(1:nargin, x,varargin{:});
validateattributes(x,{'double'},{'vector'},'ambgfun','X');
if rem(nargin,2)~=0
   y = x;
   fs = varargin{1};
   prf = varargin{2};
   numericInput = 3;
   validateattributes(prf,{'double'},{'scalar','positive','finite'},'ambgfun','PRF');
   xprf = prf;
   yprf = prf;
   isCross = false;
else
    y = varargin{1};
    fs = varargin{2};
    prf = varargin{3};
    numericInput = 4;
    validateattributes(prf,{'double'},{'vector','positive','finite','numel', 2},'ambgfun','PRF');
    prf = prf(:);
    xprf = prf(1,1);
    yprf = prf(2,1);
    isCross = true;
end
validateattributes(y,{'double'},{'vector'},'ambgfun','Y');
validateattributes(fs,{'double'},{'scalar','positive','finite'},'ambgfun','Fs');

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

cut =  validatestring(cut,{'Delay','Doppler','2D'},'ambgfun','CUT');

% validate cutvalue
xcol = x(:);  % to column form
ycol = y(:);
xleninit = numel(xcol);
yleninit = numel(ycol);
xnumperpulse = floor(fs/xprf);
xNpulseInit = xleninit/xnumperpulse;

ynumperpulse = floor(fs/yprf);
yNpulseInit = yleninit/ynumperpulse;

if abs(xNpulseInit-round(xNpulseInit))<=eps(xNpulseInit)  % remove round off error
    xNpulse = round(xNpulseInit);
else
    xNpulse = xNpulseInit;
end

cond = xNpulse < 1;
if cond
    coder.internal.errorIf(cond,'phased:ambgfun:invalidInput','X');
end

if  mod(xNpulse,1) %Npulse ~= floor(Npulse)
    xNpulse = floor(xNpulse);
    if isempty(coder.target)
       warning(message('phased:ambgfun:incompletePulseX', xNpulse));
    end
    xlen = xNpulse*xnumperpulse;
    xsub = xcol(1:xlen);    
else
    xsub = xcol;
    xlen = xleninit;
end

if abs(yNpulseInit-round(yNpulseInit))<=eps(yNpulseInit)  % remove round off error
    yNpulse = round(yNpulseInit);
else
    yNpulse = yNpulseInit;
end

cond = yNpulse < 1;
if cond
    coder.internal.errorIf(cond,'phased:ambgfun:invalidInput','Y');
end

if  mod(yNpulse,1) %Npulse ~= floor(Npulse)
    yNpulse = floor(yNpulse);
    if isempty(coder.target) && rem(nargin,2)==0
        warning(message('phased:ambgfun:incompletePulseY',yNpulse));
    end
    ylen = yNpulse*ynumperpulse;
    ysub = ycol(1:ylen);    
else
    ysub = ycol;
    ylen = yleninit;
end

if strcmp(cut, 'Delay')
    validateattributes(cutvalue,{'double'},{'scalar','finite','real','>=',-(xlen-1)/fs,'<=',(xlen-1)/fs},'ambgfun','CUTVALUE');
else
    validateattributes(cutvalue,{'double'},{'scalar','finite','real','>=',-fs/2,'<=',fs/2},'ambgfun','CUTVALUE');
end

% [EOF]
