function [w_pre,w_comb,Pi,G,C] = diagbfweights(Hchann_in,P_in,NpowSubchan_in,Popt_in)
%diagbfweights  Diagonalization beamforming weights
%   [WP,WC] = diagbfweights(HCHAN) returns the precoding weights, WP, and
%   combining weights, WC, for the channel matrix, HCHAN. These weights
%   together diagonalize the channel into independent subchannels so that
%   the result of WP*HCHAN*WC has all its off-diagonal elements equal to 0.
%
%   HCHAN can be either a matrix or a 3-dimensional array. If HCHAN is a
%   matrix, HCHAN has a size of NtxNr where Nt is number of elements in the
%   transmit array and Nr is the number of elements in the receive array.
%   If HCHAN is a 3-dimensional array, its dimension is LxNtxNr where L is
%   the number of subcarriers.
%
%   If HCHAN is a matrix, WP is an NtxNt matrix and WC has a size of NrxNr.
%   If HCHAN is a 3-dimensional array, WP is an LxNtxNt matrix and WC has a
%   size of LxNrxNr.
%
%   [WP,WC,P] = diagbfweights(HCHAN) returns the distributed power, P (in
%   linear scale), for each transmit element. If H is a matrix, P is a 1xNt
%   vector. If H is a 3-dimensional array, P is an LxNt vector.
%
%   [WP,WC,P,G] = diagbfweights(HCHAN) returns the subchannel gains in G.
%   If HCHAN is a matrix, G is a 1xR vector where R is the lessor of Nt and
%   Nr. If HCHAN is a 3-dimensional array, G has a size of LxR. G is
%   measured in linear units.
%
%   [WP,WC,P,G,C] = diagbfweights(HCHAN) returns the sum of capacity of the
%   channel in C (in bps/Hz). If HCHAN is a matrix, C is a scalar. If HCHAN
%   is a 3-dimensional array, C is an Lx1 vector.
%
%   [...] = diagbfweights(HCHAN,PT) specifies the total transmit power,
%   PT (in linear units), as a positive scalar or an L-element vector. PT
%   and P has the same units.
%
%   If PT is a scalar, then all subcarriers have the same transmit power.
%   If PT is a vector, its element specifies the transmit power for the
%   corresponding subcarrier. The total power is distributed evenly across
%   N transmit elements at each subcarrier and the result is included in
%   WP. The default value of PT is 1.
%
%   [...] = diagbfweights(HCHAN,PT,NPOW) specifies the noise power, NPOW,
%   in each receive antenna element as a scalar. All subcarriers are
%   assumed to have the same noise power. The default value of NPOW is 1.
%   NPOW shares the same unit as PT.
%
%   [...] = diagbfweights(HCHAN,PT,NPOW,POPTION) specifies how the total
%   transmit power, PT, is distributed among transmit array elements in
%   POPTION as one of 'Uniform' | 'Waterfill', where the default is
%   'Uniform'. If POPTION is 'Uniform', the transmit power is evenly
%   distributed across Nt channels. If POPTION is 'Waterfill', the transmit
%   power is distributed across the Nt channels using a waterfill 
%   algorithm.
%
%   % Examples:
%
%   % Example 1:
%   %   Given a channel matrix for a 4x4 MIMO channel, show that
%   %   diagonalization-based precoding and combining weights can achieve
%   %   spatial multiplexing, where the received signal at each element
%   %   matches the signal sent by the corresponding transmit element.
%   %   Assume the noise power of the receive channel is 0.01 watt.
%   
%   Hchan = rand(4,4);
%   x = 2*(randi(2,[20 4])-1.5);
%   [w_pre,w_comb] = diagbfweights(Hchan,4);
%   y = (x*w_pre*Hchan+0.1*randn(20,4))*w_comb;
%   for m = 1:4
%       subplot(4,1,m);
%       stem([x(:,m) y(:,m)]);
%       ylabel('Signal')
%   end
%   xlabel('Samples')
%   legend('Input','Recovered')
%
%   % Example 2:
%   %   Given a channel matrix for a 4x4 MIMO channel, show that
%   %   diagonalization-based precoding and combining weights can achieve
%   %   spatial multiplexing, where the received signal at each element
%   %   matches the signal sent by the corresponding transmit element.
%   %   Assume the noise power of receive channel is 0.1 watt. Also compare
%   %   the performance between different power allocation strategies.
%   
%   Hchan = rand(4,4);
%   x = 2*(randi(2,[20 4])-1.5);
%   [w_pre1,w_comb1,P1] = diagbfweights(Hchan,4);
%   npow = 0.1;
%   noise = sqrt(npow)*randn(20,4);
%   y1 = ((x.*P1)*w_pre1*Hchan+noise)*w_comb1;
%   [w_pre2,w_comb2,P2] = diagbfweights(Hchan,4,npow,'waterfill');
%   y2 = ((x.*P2)*w_pre2*Hchan+noise)*w_comb2;
%   for m = 1:4
%       subplot(4,1,m);
%       stem([x(:,m) y1(:,m) y2(:,m)]);
%       ylabel('Signal');
%   end
%   xlabel('Samples');
%   legend('Input','Uniform','Waterfill');
%   subplot(4,1,1);
%   title(['P_{uniform} = [',num2str(P1),...
%       '], P_{waterfill} = [',num2str(P2),']'])
%
%   See also phased, waterfill.

%   Copyright 2016-2018 The MathWorks, Inc.

%   Reference
%   [1] David Tse and Pramod Viswanath, Fundamentals of Wireless
%   Communication, Cambridge, 2005
%   [2] Arogyaswami Paulraj, et al. Introduction to Space-Time Wireless
%   Communications, Cambridge, 2003

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

narginchk(1,4);
validateattributes(Hchann_in,{'double'},{'nonnan','nonempty','finite','3d'},...
    'diagbfweights','HCHAN');

isHmatrix = ismatrix(Hchann_in);

if isHmatrix
    Hchann = Hchann_in.';
    chansize = size(Hchann);
    Nt = chansize(2);
    Nr = chansize(1);
else
    Hchann = permute(Hchann_in,[3 2 1]);  % NrxNtxL
    chansize = size(Hchann);
    L = chansize(3);
    Nt = chansize(2);
    Nr = chansize(1);
end
    
if nargin < 2
    if isHmatrix
        P = 1;
    else
        P = ones(L,1);
    end
else
    if isHmatrix
        sigdatatypes.validatePower(P_in,'diagbfweights','PT',{'scalar'});
        P = P_in;
    else
        if isscalar(P_in)
            sigdatatypes.validatePower(P_in,'diagbfweights','PT',{'scalar'});
            P = P_in*ones(L,1);
        else
            sigdatatypes.validatePower(P_in,'diagbfweights','PT',{'vector','numel',L});
            P = P_in(:);
        end
    end
end

if isHmatrix
    [U,S,V] = svd(Hchann);
    w_pre = V.';
    w_comb = conj(U);
    G = abs(diag(S).').^2;
    Nrk = sum(G>0);
else % 3D
    w_pre = zeros(L,Nt,Nt,'like',1+1i);
    w_comb = zeros(L,Nr,Nr,'like',1+1i);
    G = zeros(L,min(Nt,Nr));
    Nrk = zeros(L,1);
    for m = 1:L
        [U,S,V] = svd(Hchann(:,:,m));
        w_pre(m,:,:) = V.';
        w_comb(m,:,:) = conj(U);
        G(m,:) = abs(diag(S).').^2;
        Nrk(m) = sum(G(m,:)>0);
    end
end

if nargin < 3
    if isHmatrix
        NpowSubchan = 1;
    else
        NpowSubchan = ones(L,1);
    end
else
    sigdatatypes.validatePower(NpowSubchan_in,'diagbfweights','NPOW',{'scalar'});
    if isHmatrix
        NpowSubchan = NpowSubchan_in;
    else
        NpowSubchan = NpowSubchan_in*ones(L,1);
    end
end

if nargin < 4
    Popt = 'Uniform';
else
    Popt = validatestring(Popt_in,{'Uniform','Waterfill'},'diagbfweights','POPTION');
    Popt = convertStringsToChars(Popt);
end
    
if strcmp(Popt,'Waterfill')
    if isHmatrix
        % waterfill
        Nvec = NpowSubchan./G;
        Nvec(isinf(Nvec)) = realmax;
        Pi = waterfill(P,Nvec);
    else
        Nvec = bsxfun(@rdivide,NpowSubchan,G);
        Nvec(isinf(Nvec)) = realmax;
        Pi = waterfill(P,Nvec);
    end
else
    if isHmatrix
        Pi = repmat(P/Nt,1,Nt);
    else
        Pi = repmat(P/Nt,1,Nt);
    end
    Nvec = bsxfun(@rdivide,NpowSubchan,G);
end

if isHmatrix
    Nsub = min(Nrk,sum(Pi>0));
    C = sum(log2(1+Pi(1:Nsub)./Nvec(1:Nsub)));
else
    C = zeros(L,1);
    for m = 1:L
        Nsub = min(Nrk(m),sum(Pi(m,:)>0));
        C(m) = sum(log2(1+Pi(m,1:Nsub)./Nvec(m,1:Nsub)));
    end
end
