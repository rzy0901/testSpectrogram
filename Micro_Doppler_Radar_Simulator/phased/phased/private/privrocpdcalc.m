function Pd = privrocpdcalc(Pfa,SNR,N,type)
%PRIVROCPDCALC Calculate Pd value for ROC curves
%   Pd = PRIVROCPDCALC(Pfa,SNR,N,TYPE) returns the probability of detection
%   Pd corresponding the given probability of false alarm Pfa, signal to
%   noise ratio SNR and number of pulses N used in pulse integration.
%   TYPE specifies the signal type and can be one of the following values:
%   [ 'Real' | 'NonfluctuatingNoncoherent' | {'NonfluctuatingCoherent'} | 
%   'Swerling1' | 'Swerling2' | 'Swerling3' | 'Swerling4' ].
%
%   Both Pfa and SNR can be vectors and Pd's dimension is given by
%   [length(Pfa) length(SNR)];

%   Copyright 2008-2012 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing

%#codegen

% Checking letters was enough and ok here because type was already validated and strcmp
% is auto extrinsic.
if type(1) == 'R' %Real
    Pd =localRealROC(Pfa,SNR,N);
elseif type(1) == 'N' 
    if type(15) == 'C' %NonfluctuatingCoherent
        Pd =localNonfluctuatingCoherentROC(Pfa,SNR,N);
    else               %NonfluctuatingNoncoherent
        Pd =localNonfluctuatingNoncoherentROC(Pfa,SNR,N);
    end
else 
    if type(9) == '1'     %Swerling1
        Pd =localSwerling1ROC(Pfa,SNR,N);
    elseif type(9) == '2' %Swerling2
        Pd =localSwerling2ROC(Pfa,SNR,N);
    elseif type(9) == '3' %Swerling3
        Pd =localSwerling3ROC(Pfa,SNR,N);
    else %Swerling4
        Pd =localSwerling4ROC(Pfa,SNR,N);        
    end
end

%-------------------------------------------------------------------------
function Pd = localNonfluctuatingCoherentROC(Pfa,SNR,N) 

SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate

erfcinv2PFA = erfcinv(2*Pfa);
for k = 1:SNRlen,
    Pd(:,k) = 0.5*erfc( erfcinv2PFA-(sqrt(SNR(k)*N)) );
end

%-------------------------------------------------------------------------
function Pd = localNonfluctuatingNoncoherentROC(Pfa,SNR,N) 

if ~isempty(coder.target)
    coder.internal.assert(false,'phased:rocsnr:invalidCodegenNN');
end
SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate
T = real(gammaincinv(1-Pfa,N));
T = T(:);
SQRT2T = sqrt(2*T);

for k = 1:SNRlen,
    NX = N*SNR(k);
    NXT = NX*T;
    Pd(:,k) = marcumq(sqrt(2*NX),SQRT2T);
    Pdtemp = 0;
    if SNR(k)~=0 && SNR(k)~=inf
        for r = 2:N
            Pdtemp = Pdtemp + (T/NX).^((r-1)/2).*besseli(r-1,2.*sqrt(NXT),1);
        end
    end
    Pd(:,k) = Pd(:,k)+exp(-(T+NX)+2.*sqrt(NXT)).*Pdtemp;
end

%-------------------------------------------------------------------------
function Pd = localRealROC(Pfa,SNR,N)

SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate

erfcinv2PFA = erfcinv(2*Pfa);
SQRT2 = sqrt(2);
for k = 1:SNRlen,
    Pd(:,k) = 0.5*erfc( erfcinv2PFA-(sqrt(SNR(k)*N)/SQRT2) );
end

%-------------------------------------------------------------------------
function Pd = localSwerling1ROC(Pfa,SNR,N)

SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate
T = real(gammaincinv(1-Pfa,N));
T = T(:).';

a = 1 + N*SNR;
b = 1 + 1./(N*SNR);
for k = 1:SNRlen,
    if N == 1
        pdtmp = exp(-T/a(k));
    else
        pdtmp = 1-real(gammainc(T,N-1))+b(k)^(N-1)*real(gammainc(T/b(k),N-1)).*exp(-T/a(k));
    end
    Pd(:,k) = pdtmp;
end

%-------------------------------------------------------------------------
function Pd = localSwerling2ROC(Pfa,SNR,N)

SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate
T = real(gammaincinv(1-Pfa,N));
T = T(:).';

for k = 1:SNRlen,
    Pd(:,k) = 1 - real(gammainc(T/(1+SNR(k)),N));
end

%-------------------------------------------------------------------------
function Pd = localSwerling3ROC(Pfa,SNR,N)

SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate
T = real(gammaincinv(1-Pfa,N));
T = T(:).';

a = 1 + N*SNR/2;
b = 1 + 2./(N*SNR);
if N > 2
    c = (N-1)*log(T)-T-sum(log(1:N-2));
end
for k = 1:SNRlen,
    K0 = 1 + T/a(k) - 2*(N-2)/N/SNR(k);
    pdtmp = b(k)^(N-2)*exp(-T/a(k)).*K0;
    if N > 2
        pdtmp = exp(c-log(a(k)))+1-real(gammainc(T,N-1))+pdtmp.*real(gammainc(T/b(k),N-1));
    end
    Pd(:,k) = pdtmp;
end

%-------------------------------------------------------------------------
function Pd = localSwerling4ROC(Pfa,SNR,N)

SNRlen = numel(SNR);
Pd = zeros(numel(Pfa),SNRlen);  % preallocate
T = real(gammaincinv(1-Pfa,N));
T = T(:).';

for k = 1:SNRlen,
    tmp1 = T/(1+SNR(k)/2);
    gammaSum = real(gammainc(tmp1,N));
    for gammaI = 1:N
        c = gammaI*log(SNR(k)/2)+sum(log(N-gammaI+1:N))-sum(log(1:gammaI));
        gammaSum = gammaSum + exp(c)*real(gammainc(tmp1,N+gammaI));
    end
    pdtmp = 1-gammaSum*(1+SNR(k)/2)^(-N);
    Pd(:,k) = pdtmp;
end


% [EOF]
