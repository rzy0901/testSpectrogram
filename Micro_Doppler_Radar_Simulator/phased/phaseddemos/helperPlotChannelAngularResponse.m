function helperPlotChannelAngularResponse(txarraypos,rxarraypos,Hmat,channelstr)
% This function is only in support of ArrayProcessingForMIMOExample. It may
% be removed in a future release.

%   Copyright 2016 The MathWorks, Inc.

ang = -90:90;
txstv = steervec(txarraypos,ang);
rxstv = steervec(rxarraypos,ang);
resp = abs(txstv.'*Hmat*conj(rxstv));  % matrix NtxNr
% respspan = max(resp(:))-min(resp(:));
% C = (resp-min(resp(:)))/respspan;
h = surf(ang,ang,resp);
h.EdgeColor = 'none';
view(0,90);
xlabel('Arrival Angle (deg)');
ylabel('Departure Angle (deg)');
title(sprintf('Angular Response of %s Channel',channelstr));
axis tight

