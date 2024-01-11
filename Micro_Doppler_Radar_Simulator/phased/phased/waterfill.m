function varargout = waterfill(P_in,Nvec_in)
%waterfill   Power distribution using water-filling algorithm
%   P = waterfill(PT,NPOW) distributes the total power, PT (in linear
%   units), into multiple channels to maximize the channel capacity.
%
%   PT can be either a positive scalar or an L-element vector where L is
%   number of subcarriers. If PT is a scalar, then all subcarriers have the
%   same total power. If PT is a vector, the vector elements represent the
%   total power across channels in the corresponding subcarrier.
%
%   NPOW is an N-element vector or an LxN matrix where N is the number of
%   channels. If NPOW is a vector, each element in NPOW represents the
%   noise power (same unit as P) in the corresponding channel and that
%   noise power is the same across all subcarriers. If NPOW is a matrix,
%   each element in NPOW represents the noise power in the corresponding
%   channel at the corresponding subcarrier.
%
%   P is the power (same unit as PT) allocated in each channel and has a
%   dimension of LxN.
%
%   waterfill(PT,NPOW) plots the water-filling diagram.
%
%   % Example
%   %   Using a waterfill algorithm, distribute 10 watts of power into four
%   %   channels with noise power of 3, 7, 2, and 8 watts, respectively.
%   
%   Pt = 10;
%   Npow = [3 7 2 8];
%   P = waterfill(Pt, Npow)
%
%   % visualize
%   waterfill(Pt, Npow)
%   
%   See also phased, diagbfweights.

%   Copyright 2016-2017 The MathWorks, Inc.

%   Reference
%   [1] Thoma Cover and Joy Thomas, Elements of Information Theory, Wiley,
%   1991
%   [2] Arogyaswami Paulraj, et al. Introduction to Space-Time Wireless
%   Communications, Cambridge, 2003

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

narginchk(2,2);

sigdatatypes.validatePower(P_in,'waterfill','PT',{'vector','positive'});
sigdatatypes.validatePower(Nvec_in,'waterfill','NPOW',{'2d'});

if isscalar(P_in)
    if size(Nvec_in) == 1
        P = P_in;
        Nvec = Nvec_in;
    else
        P = repmat(P_in,size(Nvec_in,1),1);
        Nvec = Nvec_in;
    end
else
    P = P_in(:);
    if size(Nvec_in,1)==1
        Nvec = repmat(Nvec_in,numel(P_in),1);
    else
        Nvec = Nvec_in;
    end
end

cond = (numel(P)~=size(Nvec,1));
if cond
    coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers','NPOW',numel(P));
end

L = numel(P);
[NiSorted,NiIdx] = sort(Nvec,2);
NiSteps = [zeros(L,1) diff(NiSorted,1,2)];
PtotRequired = cumsum(bsxfun(@times,NiSteps,(0:size(Nvec,2)-1)),2); 
Pw = zeros(size(Nvec));
for m = 1:L
    rk_in = find(PtotRequired(m,:)<P(m),1,'last');
    if ~isempty(rk_in)
        rk = rk_in(1);
        Pw(m,1:rk) = (P(m)-PtotRequired(m,rk))/rk+(NiSorted(m,rk)-NiSorted(m,1:rk));
        [~,NiInvIdx] = sort(NiIdx(m,:));
        Pw(m,:) = Pw(m,NiInvIdx);
    end
end
if nargout > 0
    varargout{1} = Pw;
else
    cond = isempty(coder.target);
    if ~cond
        coder.internal.assert(cond,'phased:ambgfun:invalidCodegenOutput');
    end
    Nvec_plot = Nvec.';
    Pw_plot = Pw.';
    b = bar([[Nvec_plot(:);nan], [Pw_plot(:);nan]],1,'Stacked');
    b(1).FaceColor = 'y';
    b(2).FaceColor = 'b';
    xlim([0.5,numel(Nvec_plot(:))+0.5]);
    xlabel('Channel');
    ylabel('Power');
    legend('Channel noise','Allocated power',...
        'Location','NorthOutside','Orientation','Horizontal');
end



