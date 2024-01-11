function varargout = plotRadiationPattern(UseArray,H,varargin)
%This function is for internal use only. It may be removed in the future.

%PLOTRADIATIONPATTERN plot radiation pattern
%   PLOTRADIATIONPATTERN is used by element and array to plot response
%   patterns.

%   Copyright 2011-2017 The MathWorks, Inc.

freq = varargin{1};

% check if the right arguments are passed with Array Object and Element
% system object
if UseArray
    v = varargin{2};
    pvpair = varargin(3:end);
else
    pvpair = varargin(2:end);
end


format = '';
respcut = '';
cutangle = [];
unit = '';
normalizeresponse = [];
overlayfreq = [];
steerangle = [];
polarization = '';
azimuthangles = [];
elevationangles = [];
ugrid = [];
vgrid = [];
weights = [];
funcname = '';
twosidedelevation = false;
elementweights = [];

sigutils.pvparse(pvpair{:});

directivityflag = strcmp(unit,'dbi');
if directivityflag
    unit = 'db';
end

%check polarization capability
if ~directivityflag
    polflag = ~strcmp(polarization,'None');
else
    polflag = false;
end

if UseArray
    if isa(H,'phased.ReplicatedSubarray')
        element = getElementHandle(H.Subarray);
    elseif isa(H,'phased.PartitionedArray')
        element = getElementHandle(H.Array);
    else
        element = getElementHandle(H);
    end
else
    element = H;
end

if iscell(element)  % heterogeneous array case
    for idx = 1:numel(element)
        if isElementFromAntenna(element{idx})
            isValidFreq = true;
        else
            [isValidFreq,frange] = searchValidFrequency(element{idx},freq);
        end
        if any(~isValidFreq)
            warning(message('phased:system:array:NoResponseFrequency', ...
                num2str(idx),num2str(frange)));
        end
    end
else
    if isElementFromAntenna(element)
        isValidFreq = true;
    else
        [isValidFreq,frange] = searchValidFrequency(element,freq);
    end
    if any(~isValidFreq)
        warning(message('phased:system:element:NoResponseFrequency',num2str(frange)));
    end
end

if UseArray && isa(H,'phased.internal.AbstractSubarray') && ...
        ~strncmp(H.SubarraySteering,'None',1)
    SteerSubarray = true;
    if strncmp(H.SubarraySteering,'Custom',1)
        steerangle = elementweights; % reuse steerangle variable
    end
else
    SteerSubarray = false;
end

if UseArray
    N = getDOF(H);
    % Validate that weights is in a proper vector or matrix format.
    if strncmpi(respcut,'3d',2)
        validateattributes(weights,{'numeric'},{'vector','nonnan','finite',...
            'nonempty','size',[N 1]},funcname,'Weights');
    end
    if iscolumn(weights)
        % If weights is entered as a column vector, ensure that
        % the number of rows is equal to the number of elements.
        validateattributes(weights,{'numeric'},{'nonnan','finite',...
            'nonempty','size',[N 1]},funcname,'Weights');
    else
        % If weights is entered as a matrix, ensure that it contains
        % exactly one column per frequency point and one row per
        % element.
        if(numel(freq)==1)
            rfreq = repmat(freq,1,size(weights,2));
            validateattributes(weights,{'numeric'},{'nonnan','finite',...
                'nonempty','size',[N length(rfreq)]},funcname,'Weights');
        else
        validateattributes(weights,{'numeric'},{'nonnan','finite',...
            'nonempty','size',[N length(freq)]},funcname,'Weights');
        end
    end
    temp = cloneSensor(H);
    if ~directivityflag
        har = phased.ArrayResponse('SensorArray',temp,...
            'PropagationSpeed',v,'WeightsInputPort',true,...
            'EnablePolarization',polflag);
    elseif ~strncmpi(format,'uv',2) && strncmpi(respcut,'3',1)
        % for complete 3D directivity plot in non-UV coordinates, get the
        % pattern then do the directivity computation directly to avoid
        % recompute the pattern
        complete3dflag = isComplete3DSpace(azimuthangles,elevationangles);
        if complete3dflag
            har = phased.ArrayResponse('SensorArray',temp,...
                'PropagationSpeed',v,'WeightsInputPort',true);
        else
            har = phased.internal.Directivity('Sensor',temp,...
                'PropagationSpeed',v,'WeightsInputPort',true);
        end
    else
        har = phased.internal.Directivity('Sensor',temp,...
            'PropagationSpeed',v,'WeightsInputPort',true);
    end        
    clear temp;
else
    if ~directivityflag
        har = clone(H);
    else
        har = phased.internal.Directivity('Sensor',H);
    end
    release(har);
end

if ~strncmpi(format,'uv',2)
    
      
    azangle = azimuthangles;
    elangle = elevationangles;
           
    if strncmpi(respcut,'a',1)
        
        multiangsliceflag = (numel(cutangle)>1);
        numcut = numel(cutangle);
        
        if ~multiangsliceflag
           
            elangle = cutangle*ones(size(azangle));

            if UseArray
                if iscolumn(weights)
                    % If weights is entered as a vector, apply the same
                    % set of weights to each frequency when step is called
                    if SteerSubarray
                        resp = step(har,freq,[azangle; elangle],weights,steerangle);
                    else
                        resp = step(har,freq,[azangle; elangle],weights);
                    end
                    if polflag
                        if strcmp(polarization,'H')
                            resp = resp.H;
                        elseif strcmp(polarization,'V')% V
                            resp = resp.V;
                        else
                            resp = hypot(resp.H,resp.V);
                        end
                    end
                else
                    % If weights is entered as a matrix, call step with a
                    % separate set of frequency points for each column of
                    % weights.
                    if numel(freq)==1
                        freq = rfreq;
                    end
                    if SteerSubarray
                        if polflag
                            for nfreq = length(freq):-1:1
                                tempresp = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq),steerangle);
                                resp.H(:,nfreq) = tempresp.H;
                                resp.V(:,nfreq) = tempresp.V;
                            end
                            if strcmp(polarization,'H')
                                resp = resp.H;
                            elseif strcmp(polarization,'V')% V
                                resp = resp.V;
                            else
                                resp = hypot(resp.H,resp.V);
                            end
                        else
                            for nfreq = length(freq):-1:1
                                resp(:,nfreq) = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq),steerangle);
                            end
                        end
                    else
                        if polflag
                            for nfreq = length(freq):-1:1
                                tempresp = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq));
                                resp.H(:,nfreq) = tempresp.H;
                                resp.V(:,nfreq) = tempresp.V;
                            end
                            if strcmp(polarization,'H')
                                resp = resp.H;
                            elseif strcmp(polarization,'V')% V
                                resp = resp.V;
                            else
                                resp = hypot(resp.H,resp.V);
                            end
                        else
                            for nfreq = length(freq):-1:1
                                resp(:,nfreq) = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq));
                            end
                        end
                    end
                end
            else
                if polflag
                    resp = step(har,freq,[azangle; elangle]);
                    if strcmp(polarization,'H')
                        resp = resp.H;
                    elseif strcmp(polarization,'V')% V
                        resp = resp.V;
                    else
                        resp = hypot(resp.H,resp.V);
                    end
                else
                    if ~directivityflag
                        resp = step(har,freq,[azangle; elangle]);
                        if isPolarizationEnabled(har)
                            resp = hypot(resp.H,resp.V);
                        end
                    else
                        resp = step(har,freq,[azangle; elangle]);
                    end
                end
            end
            if directivityflag
                hresp = phased.internal.DirectivityPattern2D(...
                    'Pattern',resp,'Angle',azangle.','SliceDir','Az',...
                    'Freq',freq,'CutAngle',cutangle);
                plotdata.resp = resp;
            else
                hresp = phased.internal.RespPattern2D(...
                    'Pattern',abs(resp),'Angle',azangle.','SliceDir','Az',...
                    'Freq',freq,'CutAngle',cutangle);
                plotdata.resp = abs(resp);
            end
            plotdata.az = azangle(:).';
            plotdata.el = cutangle(:).';
        
        else
            % cutangle is vector, freq is scalar, weights is also vector
            
            elangle = cutangle(:)*ones(size(azangle));
            
            if UseArray
                
                validateattributes(weights,{'double'},{'column'},funcname,'Weights');
                if SteerSubarray
                    for m = numcut:-1:1
                        resp(:,m) = step(har,freq,[azangle; elangle(m,:)],weights,steerangle);
                    end
                else
                    for m = numcut:-1:1
                        resp(:,m) = step(har,freq,[azangle; elangle(m,:)],weights);
                    end
                end
                if polflag
                    if strcmp(polarization,'H')
                        resp = [resp.H];
                    elseif strcmp(polarization,'V')% V
                        resp = [resp.V];
                    else
                        resp = hypot([resp.H],[resp.V]);
                    end
                end
                
            else
                if polflag
                    for m = numcut:-1:1
                        resp(:,m) = step(har,freq,[azangle; elangle(m,:)]);
                    end
                    if strcmp(polarization,'H')
                        resp = [resp.H];
                    elseif strcmp(polarization,'V')% V
                        resp = [resp.V];
                    else
                        resp = hypot([resp.H],[resp.V]);
                    end
                else
                    if ~directivityflag
                        for m = numcut:-1:1
                            resp(:,m) = step(har,freq,[azangle; elangle(m,:)]);
                        end
                        if isPolarizationEnabled(har)
                            resp = hypot(resp.H,resp.V);
                        end
                    else
                        for m = numcut:-1:1
                            resp(:,m) = step(har,freq,[azangle; elangle(m,:)]);
                        end
                    end
                end
            end
            if directivityflag
                hresp = phased.internal.DirectivityPattern2D(...
                    'Pattern',resp,'Angle',azangle.','SliceDir','Az',...
                    'Freq',freq,'CutAngle',cutangle);
                plotdata.resp = resp;
            else
                hresp = phased.internal.RespPattern2D(...
                    'Pattern',abs(resp),'Angle',azangle.','SliceDir','Az',...
                    'Freq',freq,'CutAngle',cutangle);
                plotdata.resp = abs(resp);
            end
            plotdata.az = azangle(:).';
            plotdata.el = cutangle(:).';
        end
    elseif strncmpi(respcut,'e',1)
      
        multiangsliceflag = (numel(cutangle)>1);
        numcut = numel(cutangle);
        
        if ~multiangsliceflag
            if strcmp(funcname,'plotResponse') || ...
                    (strcmp(funcname,'patternElevation') && ...
                    twosidedelevation)
                elangle1 = -180-fliplr(elangle(elangle>=-90&elangle<=0));
                elangle2 = elangle(elangle>-90&elangle<90);
                elangle3 = 180-fliplr(elangle(elangle>=0&elangle<=90));
                elangle_check = [elangle1 elangle2 elangle3];
                elangle = [wrapTo180(180-elangle1) elangle2 wrapTo180(180-elangle3)];

                azangle1 = wrapTo180(cutangle+180)*ones(size(elangle1));
                azangle2 = cutangle*ones(size(elangle2));
                azangle3 = wrapTo180(cutangle+180)*ones(size(elangle3));
                azangle = [azangle1 azangle2 azangle3];
            else
                azangle = cutangle*ones(size(elangle));
            end

            if UseArray
                if iscolumn(weights)
                    % If weights is entered as a vector, apply the same
                    % set of weights to each frequency when step is called
                    if SteerSubarray
                        resp = step(har,freq,[azangle; elangle],weights,steerangle);
                    else
                        resp = step(har,freq,[azangle; elangle],weights);
                    end
                    if polflag
                        if strcmp(polarization,'H')
                            resp = resp.H;
                        elseif strcmp(polarization,'V')% V
                            resp = resp.V;
                        else
                            resp = hypot(resp.H,resp.V);
                        end
                    end
                else
                    % If weights is entered as a matrix, call step with a
                    % separate set of frequency points for each column of
                    % weights.
                    if(numel(freq)==1)
                        freq = rfreq;
                    end
                    if SteerSubarray
                        if polflag
                            for nfreq = length(freq):-1:1
                                tempresp = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq),steerangle);
                                resp.H(:,nfreq) = tempresp.H;
                                resp.V(:,nfreq) = tempresp.V;
                            end
                            if strcmp(polarization,'H')
                                resp = resp.H;
                            elseif strcmp(polarization,'V')% V
                                resp = resp.V;
                            else
                                resp = hypot(resp.H,resp.V);
                            end
                        else
                            for nfreq = length(freq):-1:1
                                resp(:,nfreq) = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq),steerangle);
                            end
                        end
                    else
                        if polflag
                            for nfreq = length(freq):-1:1
                                tempresp = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq));
                                resp.H(:,nfreq) = tempresp.H;
                                resp.V(:,nfreq) = tempresp.V;
                            end
                            if strcmp(polarization,'H')
                                resp = resp.H;
                            elseif strcmp(polarization,'V')% V
                                resp = resp.V;
                            else
                                resp = hypot(resp.H,resp.V);
                            end
                        else
                            for nfreq = length(freq):-1:1
                                resp(:,nfreq) = ...
                                    step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq));
                            end
                        end
                    end
                end
            else
                if polflag
                    resp = step(har,freq,[azangle; elangle]);
                    if strcmp(polarization,'H')
                        resp = resp.H;
                    elseif strcmp(polarization,'V')% V
                        resp = resp.V;
                    else
                        resp = hypot(resp.H,resp.V);
                    end
                else
                    if ~directivityflag
                        resp = step(har,freq,[azangle; elangle]);
                        if isPolarizationEnabled(har)
                            resp = hypot(resp.H,resp.V);
                        end
                    else
                        resp = step(har,freq,[azangle; elangle]);
                    end
                end
            end
            if strcmp(funcname,'plotResponse') || ...
                    (strcmp(funcname,'patternElevation') && ...
                    twosidedelevation)
                min_check = any(elangle_check==-90);
                max_check = any(elangle_check==90);
                if(min_check&&max_check)
                    elangle = [elangle1.';elangle2.';elangle3.'];
                elseif (min_check&&~max_check)
                    elangle = [elangle1.';elangle2.';NaN;elangle3.'];
                    resp = [resp(1:numel(elangle1),:);resp(numel(elangle1)+1:numel(elangle1)+numel(elangle2),:);...
                                NaN*ones(1,size(resp,2));resp(end-numel(elangle3)+1:end,:)];
                elseif (max_check&&~min_check)
                    elangle = [elangle1.';NaN;elangle2.';elangle3.'];
                    resp = [resp(1:numel(elangle1),:);NaN*ones(1,size(resp,2));resp(numel(elangle1)+1:numel(elangle1)+numel(elangle2),:);...
                                resp(end-numel(elangle3)+1:end,:)];
                else
                    elangle = [elangle1.';NaN;elangle2.';NaN;elangle3.'];
                    resp = [resp(1:numel(elangle1),:);NaN*ones(1,size(resp,2));resp(numel(elangle1)+1:numel(elangle1)+numel(elangle2),:);...
                                NaN*ones(1,size(resp,2));resp(end-numel(elangle3)+1:end,:)];
                end
            else
                elangle = elangle.';
            end
            if directivityflag
                hresp = phased.internal.DirectivityPattern2D(...
                    'Pattern',resp,...
                    'Angle',elangle,...
                    'SliceDir','El',...
                    'Freq',freq,...
                    'CutAngle',cutangle);
                plotdata.resp = resp;
            else
                hresp = phased.internal.RespPattern2D(...
                    'Pattern',abs(resp),...
                    'Angle',elangle,...
                    'SliceDir','El',...
                    'Freq',freq,...
                    'CutAngle',cutangle);
                plotdata.resp = abs(resp);
            end
            plotdata.az = cutangle(:).';
            plotdata.el = elangle(:).';
            if strcmp(funcname,'patternElevation') && ~twosidedelevation
                valid_el_idx = ~isnan(plotdata.el) & plotdata.el>=-90 & ...
                    plotdata.el<=90;
                plotdata.el = plotdata.el(valid_el_idx);
                plotdata.resp = plotdata.resp(valid_el_idx,:);
            end
        
        else
            % cutangle is vector, freq is scalar, weights is also vector
            
            if strcmp(funcname,'plotResponse') || ...
                    (strcmp(funcname,'patternElevation') && ...
                    twosidedelevation)
                elangle1 = -180-fliplr(elangle(elangle>=-90&elangle<=0));
                elangle2 = elangle(elangle>-90&elangle<90);
                elangle3 = 180-fliplr(elangle(elangle>=0&elangle<=90));
                elangle_check = [elangle1 elangle2 elangle3];
                elangle = [wrapTo180(180-elangle1) elangle2 wrapTo180(180-elangle3)];

                azangle1 = wrapTo180(cutangle(:)+180)*ones(size(elangle1));
                azangle2 = cutangle(:)*ones(size(elangle2));
                azangle3 = wrapTo180(cutangle(:)+180)*ones(size(elangle3));
                azangle = [azangle1 azangle2 azangle3];
            else
                azangle = cutangle(:)*ones(size(elangle));
            end

            if UseArray
                
                validateattributes(weights,{'double'},{'column'},funcname,'Weights');
                if SteerSubarray
                    for m = numcut:-1:1
                        resp(:,m) = step(har,freq,[azangle(m,:); elangle],weights,steerangle);
                    end
                else
                    for m = numcut:-1:1
                        resp(:,m) = step(har,freq,[azangle(m,:); elangle],weights);
                    end
                end
                if polflag
                    if strcmp(polarization,'H')
                        resp = [resp.H];
                    elseif strcmp(polarization,'V')% V
                        resp = [resp.V];
                    else
                        resp = hypot([resp.H],[resp.V]);
                    end
                end
                
            else
                if polflag
                    for m = numcut:-1:1
                        resp(:,m) = step(har,freq,[azangle(m,:); elangle]);
                    end
                    if strcmp(polarization,'H')
                        resp = [resp.H];
                    elseif strcmp(polarization,'V')% V
                        resp = [resp.V];
                    else
                        resp = hypot([resp.H],[resp.V]);
                    end
                else
                    if ~directivityflag
                        for m = numcut:-1:1
                            resp(:,m) = step(har,freq,[azangle(m,:); elangle]);
                        end 
                        if isPolarizationEnabled(har)
                            resp = hypot(resp.H,resp.V);
                        end
                    else
                        for m = numcut:-1:1
                            resp(:,m) = step(har,freq,[azangle(m,:); elangle]);
                        end
                    end
                end
            end
            if strcmp(funcname,'plotResponse') || ...
                    (strcmp(funcname,'patternElevation') && ...
                    twosidedelevation)
                min_check = any(elangle_check==-90);
                max_check = any(elangle_check==90);
                if(min_check&&max_check)
                    elangle = [elangle1.';elangle2.';elangle3.'];
                elseif (min_check&&~max_check)
                    elangle = [elangle1.';elangle2.';NaN;elangle3.'];
                    resp = [resp(1:numel(elangle1),:);resp(numel(elangle1)+1:numel(elangle1)+numel(elangle2),:);...
                                NaN*ones(1,size(resp,2));resp(end-numel(elangle3)+1:end,:)];
                elseif (max_check&&~min_check)
                    elangle = [elangle1.';NaN;elangle2.';elangle3.'];
                    resp = [resp(1:numel(elangle1),:);NaN*ones(1,size(resp,2));resp(numel(elangle1)+1:numel(elangle1)+numel(elangle2),:);...
                                resp(end-numel(elangle3)+1:end,:)];
                else
                    elangle = [elangle1.';NaN;elangle2.';NaN;elangle3.'];
                    resp = [resp(1:numel(elangle1),:);NaN*ones(1,size(resp,2));resp(numel(elangle1)+1:numel(elangle1)+numel(elangle2),:);...
                                NaN*ones(1,size(resp,2));resp(end-numel(elangle3)+1:end,:)];
                end
            else
                elangle = elangle.';
            end
            if directivityflag
                hresp = phased.internal.DirectivityPattern2D(...
                    'Pattern',resp,...
                    'Angle',elangle,...
                    'SliceDir','El',...
                    'Freq',freq,...
                    'CutAngle',cutangle);
                plotdata.resp = resp;
            else
                hresp = phased.internal.RespPattern2D(...
                    'Pattern',abs(resp),...
                    'Angle',elangle,...
                    'SliceDir','El',...
                    'Freq',freq,...
                    'CutAngle',cutangle);
                plotdata.resp = abs(resp);
            end
            plotdata.az = cutangle(:).';
            plotdata.el = elangle(:).';
            if strcmp(funcname,'patternElevation') && ~twosidedelevation
                valid_el_idx = ~isnan(plotdata.el) & plotdata.el>=-90 & ...
                    plotdata.el<=90;
                plotdata.el = plotdata.el(valid_el_idx);
                plotdata.resp = plotdata.resp(valid_el_idx,:);
            end
        end

    else % 3D
        azangle = azangle.';
        elangle = elangle.';
        azangle_len = size(azangle,1); 
        elangle_len = size(elangle,1);
        resp = ones(elangle_len,azangle_len);
        
        temp = ones(azangle_len,1);
        temp_ang = [azangle temp];
        for m = 1:elangle_len
            temp_ang(:,2) = elangle(m)*temp;
            if UseArray
                if SteerSubarray
                    if polflag
                        tempresp = step(har,freq,temp_ang.',weights,steerangle);
                        if strcmp(polarization,'H')
                            resp(m,:) = tempresp.H;
                        elseif strcmp(polarization,'V')% V
                            resp(m,:) = tempresp.V;
                        else
                            resp(m,:) = hypot(tempresp.H,tempresp.V);
                        end
                    else
                        resp(m,:) = step(har,freq,temp_ang.',weights,steerangle);
                    end
                else
                    if polflag
                        tempresp = step(har,freq,temp_ang.',weights);
                        if strcmp(polarization,'H')
                            resp(m,:) = tempresp.H;
                        elseif strcmp(polarization,'V')% V
                            resp(m,:) = tempresp.V;
                        else
                            resp(m,:) = hypot(tempresp.H,tempresp.V);
                        end
                    else
                        resp(m,:) = step(har,freq,temp_ang.',weights);
                    end
                end
            else
                if polflag
                    tempresp = step(har,freq,temp_ang.');
                    if strcmp(polarization,'H')
                        resp(m,:) = tempresp.H;
                    elseif strcmp(polarization,'V')% V
                        resp(m,:) = tempresp.V;
                    else
                        resp(m,:) = hypot(tempresp.H,tempresp.V);
                    end
                else
                    if ~directivityflag
                        tempresp = step(har,freq,temp_ang.');
                        if isPolarizationEnabled(har)
                            resp(m,:) = hypot(tempresp.H,tempresp.V);
                        else
                            resp(m,:) = tempresp;
                        end
                    else
                        resp(m,:) = step(har,freq,temp_ang.');
                    end
                end
            end
        end
        
        if directivityflag
            if UseArray && complete3dflag
                resp_real = phased.internal.directivityFromPattern(...
                    resp,azangle,elangle);
                hresp = phased.internal.DirectivityPattern3D(...
                    'Pattern',resp_real,'AzAngle',azangle,'ElAngle',elangle);
                plotdata.resp = resp_real;
            else
                hresp = phased.internal.DirectivityPattern3D(...
                    'Pattern',resp,'AzAngle',azangle,'ElAngle',elangle);
                plotdata.resp = resp;
            end
        else
            hresp = phased.internal.RespPattern3D(...
                'Pattern',abs(resp),'AzAngle',azangle,'ElAngle',elangle);
            plotdata.resp = abs(resp);
        end
        plotdata.az = azangle(:).';
        plotdata.el = elangle(:).';
    end
else  %UV
    if strncmpi(respcut,'U',1)
        u = ugrid;
        azangle = asind(u);
        elangle = zeros(size(azangle));
        
        if UseArray
            if iscolumn(weights)
                % If weights is entered as a vector, apply the same
                % set of weights to each frequency when step is called
                if SteerSubarray
                    resp = step(har,freq,[azangle; elangle],weights,steerangle);
                else
                    resp = step(har,freq,[azangle; elangle],weights);
                end
                if polflag
                    if strcmp(polarization,'H')
                        resp = resp.H;
                    elseif strcmp(polarization,'V')% V
                        resp = resp.V;
                    else
                        resp = hypot(resp.H,resp.V);
                    end
                end
            else
                % If weights is entered as a matrix, call step with a
                % separate set of frequency points for each column of
                % weights.
                if(numel(freq)==1)
                    freq = rfreq;
                end
                if SteerSubarray
                    if polflag
                        for nfreq = length(freq):-1:1
                            tempresp = ...
                                step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq),steerangle);
                            resp.H(:,nfreq) = tempresp.H;
                            resp.V(:,nfreq) = tempresp.V;
                        end
                        if strcmp(polarization,'H')
                            resp = resp.H;
                        elseif strcmp(polarization,'V')% V
                            resp = resp.V;
                        else
                            resp = hypot(resp.H,resp.V);
                        end
                    else
                        for nfreq = length(freq):-1:1
                            resp(:,nfreq) = ...
                                step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq),steerangle);
                        end
                    end
                else
                    if polflag
                        for nfreq = length(freq):-1:1
                            tempresp = ...
                                step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq));
                            resp.H(:,nfreq) = tempresp.H;
                            resp.V(:,nfreq) = tempresp.V;
                        end
                        if strcmp(polarization,'H')
                            resp = resp.H;
                        elseif strcmp(polarization,'V')% V
                            resp = resp.V;
                        else
                            resp = hypot(resp.H,resp.V);
                        end
                    else
                        for nfreq = length(freq):-1:1
                            resp(:,nfreq) = ...
                                step(har,freq(nfreq),[azangle; elangle],weights(:,nfreq));
                        end
                    end
                end
            end
        else
            if polflag
                resp = step(har,freq,[azangle; elangle]);
                if strcmp(polarization,'H')
                    resp = resp.H;
                elseif strcmp(polarization,'V')% V
                    resp = resp.V;
                else
                    resp = hypot(resp.H,resp.V);
                end
            else
                if ~directivityflag
                    resp = step(har,freq,[azangle; elangle]);
                    if isPolarizationEnabled(har)
                        resp = hypot(resp.H,resp.V);
                    end
                else
                    resp = step(har,freq,[azangle; elangle]);
                end
            end
        end
        if directivityflag
            hresp = phased.internal.DirectivityPatternUV2D(...
                'Pattern',resp,'U',u.','Freq',freq);
            plotdata.resp = resp;
        else
            hresp = phased.internal.RespPatternUV2D(...
                'Pattern',abs(resp),'U',u.','Freq',freq);
            plotdata.resp = abs(resp);
        end
        plotdata.az = u(:).';
        plotdata.el = 0;
    else  %3D
        u = ugrid;
        v = vgrid;
        [u_plot,v_plot] = meshgrid(u',v');
        index = find(hypot(u_plot,v_plot)<=1);
        azel_plot = zeros(2,numel(u_plot));
        azel_plot(:,index) = uv2azel([u_plot(index).';v_plot(index).']);
        az_plot = reshape(azel_plot(1,:).',size(u_plot));
        el_plot = reshape(azel_plot(2,:).',size(v_plot));
        resp = NaN(size(az_plot));
        
         for m = 1:size(resp,1)
            temp_ang = [az_plot(m,:); el_plot(m,:)];
            if UseArray
                if SteerSubarray
                    if polflag
                        tempresp = step(har,freq,temp_ang,weights,steerangle);
                        if strcmp(polarization,'H')
                            resp(m,:) = tempresp.H;
                        elseif strcmp(polarization,'V') % V
                            resp(m,:) = tempresp.V;
                        else
                            resp(m,:) = hypot(tempresp.H,tempresp.V);
                        end
                    else
                        resp(m,:) = step(har,freq,temp_ang,weights,steerangle);
                    end
                else
                    if polflag
                        tempresp = step(har,freq,temp_ang,weights);
                        if strcmp(polarization,'H')
                            resp(m,:) = tempresp.H;
                        elseif strcmp(polarization,'V') % V
                            resp(m,:) = tempresp.V;
                        else
                            resp(m,:) = hypot(tempresp.H,tempresp.V);
                        end
                    else
                        resp(m,:) = step(har,freq,temp_ang,weights);
                    end
                end
            else
                if polflag
                    tempresp = step(har,freq,temp_ang);
                    if strcmp(polarization,'H')
                        resp(m,:) = tempresp.H;
                    elseif strcmp(polarization,'V')% V
                        resp(m,:) = tempresp.V;
                    else
                        resp(m,:) = hypot(tempresp.H,tempresp.V);
                    end
                else
                    if ~directivityflag
                        tempresp = step(har,freq,temp_ang);
                        if isPolarizationEnabled(har)
                            resp(m,:) = hypot(tempresp.H,tempresp.V);
                        else
                            resp(m,:) = tempresp;
                        end
                    else
                        resp(m,:) = step(har,freq,temp_ang);
                    end
                end
            end
        end
        resp(hypot(u_plot,v_plot)>1)=NaN;

        if(size(resp,2)<=3)
            if(any(u==-1)||any(u==1))
              error(message('phased:system:array:NeedPoints','UGrid'));  
            end
        elseif (size(resp,1)<=3)
            if(any(v==-1)||any(v==1))
              error(message('phased:system:array:NeedPoints','VGrid'));  
            end
        end
        
        if directivityflag
            hresp = phased.internal.DirectivityPatternUV3D(...
                'Pattern',resp,'u',u_plot,'v',v_plot);
            plotdata.resp = resp;
        else
            hresp = phased.internal.RespPatternUV3D(...
                'Pattern',abs(resp),'u',u_plot,'v',v_plot);
            plotdata.resp = abs(resp);
        end

        plotdata.az = u(:).';
        plotdata.el = v(:).';
    end
    
end

if ~directivityflag
    plotdata.resp = phased.internal.computePlotPattern(...
        plotdata.resp,normalizeresponse,unit);
end
            
if strncmpi(format,'p',1)
    if nargout
        if ~strncmp(funcname,'pattern',7)
            plotdata.fighandle = polar(hresp,'Units',unit,...
                'NormalizeResp',normalizeresponse);
        end
        varargout{1} = plotdata;
    else
        polar(hresp,'Units',unit,...
            'NormalizeResp',normalizeresponse);
    end
elseif overlayfreq %#ok<BDSCI,BDLGI>
    if nargout
        if ~strncmp(funcname,'pattern',7)
            plotdata.fighandle = plot(hresp,'Units',unit,...
                'NormalizeResp',normalizeresponse);
        end
        varargout{1} = plotdata;
    else
        plot(hresp,'Units',unit,...
                'NormalizeResp',normalizeresponse);
    end

else
    % If a line plot with without overlaid frequencies is
    % specified, invoke RespPattern2D's WATERFALL method.
    if nargout
        if ~strncmp(funcname,'pattern',7)
            plotdata.fighandle = waterfall(hresp,'Units',unit,...
                'NormalizeResp',normalizeresponse);
        end
        varargout{1} = plotdata;
    else
        waterfall(hresp,'Units',unit,...
            'NormalizeResp',normalizeresponse);
    end
end

end

function flag = isComplete3DSpace(azimuthangles,elevationangles)

stepthreshold = 5;

diff_az = diff(azimuthangles);
delta_az = mean(diff_az);
[min_az,max_az] = bounds(azimuthangles);
span_az = max_az-min_az;

diff_el = diff(elevationangles);
delta_el = mean(diff_el);
[min_el,max_el] = bounds(elevationangles);
span_el = max_el-min_el;

relerr = 0.1;
if max(abs(diff_az-delta_az))/delta_az < relerr && ... % nearly uniform sampling in az
        max(abs(diff_el-delta_el))/delta_el < relerr && ... % nearly uniform sampling in el
        360-span_az < min((1+relerr)*delta_az,stepthreshold) && ... % span 360 az
        180-span_el < min((1+relerr)*delta_el,stepthreshold)
    flag = true;
else
    flag = false;
end


end

% [EOF]
