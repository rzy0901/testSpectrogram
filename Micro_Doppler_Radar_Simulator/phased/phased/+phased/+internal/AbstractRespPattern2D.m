classdef (Hidden) AbstractRespPattern2D < phased.internal.AbstractRespPattern
%This class is for internal use only. It may be removed in the future.

%AbstractRespPattern2D   Define the AbstractRespPattern2D class.

%   Copyright 2011-2016 The MathWorks, Inc.

    properties

        %Pattern - Response pattern
        %   Pattern is a matrix containing the samples of the magnitude
        %   response pattern or patterns in each column.  The default
        %   pattern is a slice of isotropic pattern.
        Pattern
        %Freq - Row vector of frequencies corresponding to each pattern
        %   Freq holds the frequency at which each corresponding pattern was
        %   computed.  The default frequency is 1GHz.
        Freq
    end
    
    methods (Access = protected, Abstract)
        xlbl = getXLabel(obj,plotoptionobj);
        plotoptionobj = getPlotOption(obj,varargin);
        plotTitle = genPlotTitle(obj)    
        validatePattern(obj,pattern)
        angles = getPlotAngles(obj)
    end

    methods (Access = protected)

        function obj = AbstractRespPattern2D
            %AbstractRespPattern2D   Construct the AbstractRespPattern2D
            %class.

        end
    end
    
    methods (Access = protected)
        
        function g = getSliceGroup(obj) %#ok<MANU>
            g = 'Frequency';
        end
        
        function annotatePlot(obj,plotoptionobj)
            %ANNOTATEPLOT Annotate response pattern plots
            phased.internal.AbstractRespPattern.annotatePlotSetup;
            resplbl = getRespLabel(plotoptionobj);
            
            % Label the Y axis and Z axis for waterfall plots
            if strcmp(plotoptionobj.PlotType,'waterfall')
                zlabel(resplbl);
                ylabel(getGroupLabel(plotoptionobj,getSliceGroup(obj)))
            else
                ylabel(resplbl);
            end
            
            xlabel(getXLabel(obj,plotoptionobj));
            title(plotoptionobj.Title);
        end
        
        function addLegend(obj,hlines,varargin)
            %ADDLEGEND Add a legend to label each pattern's frequency
            %   ADDLEGEND(obj,hlines) adds a legend to the current axes
            %   with frequency labels corresponding to each line in handle
            %   array HLINES. If a single pattern is being plotted to an
            %   axes that does not contain a pattern previously labelled in
            %   a legend by AddLegend or saved to the AppData, the legend
            %   label and line handle are saved to the axes' AppData.  The
            %   AppData entries are used to populate a legend if ADDLEGEND
            %   is called again with the same current axes and the axes'
            %   hold state is set to 'on'.
            %
            %   ADDLEGEND(obj,hlines,hExistingLines,existingLabels) adds a
            %   legend with frequency labels corresponding to each line in
            %   HLINES.  Existing line handles in array HEXISTINGLINES and
            %   existing legend labels in EXISTINGLABELS are included in
            %   the legend.
            
            % Convert frequency values from Hz units to engineering units
            precision = 5; %number of digits to use after the decimal point.
            [convertedVals,unitPrefix] = convert2engstrs(obj.Freq,precision);
            
            % Ensure that CONVERTEDVALS is of class cell.  It will be a string
            % when convert2engstrs is called with a single frequency.
            if ~iscell(convertedVals)
                convertedVals = {convertedVals};
            end
            
            % Create a cell array to store legend labels
            labels = cell(size(convertedVals));
            % labels = cell(size(obj.Freq));

            % Determine if existing legend entries are to be appended
            appendFlag = nargin>2;
                
            % Determine if a single line plot is being created in an axes
            % without an existing legend.
            isNewSingleLinePlot =  ~appendFlag && isscalar(hlines);

            % Append each frequency value with its corresponding unit
            % prefix.  For multi-frequency plots, each line is given a
            % descriptive Tag.
            if isNewSingleLinePlot
                labels = {sprintf('%s %sHz',convertedVals{:},unitPrefix)};
            else
                for n=1:numel(labels)
                    labels{n} = ...
                        sprintf('%s %sHz',convertedVals{n},unitPrefix);
                    set(hlines(n),'tag',...
                        sprintf('Freq%s',num2str(obj.Freq(n))));
                end
            end
            
            %Label the plots according to Weights if there exists one
            %frequency corresponding to multiple patterns
            if isNewSingleLinePlot
                labels = {sprintf('%s %sHz',convertedVals{:},unitPrefix)};
            else
                if (~isscalar(obj.Freq)&& (numel(unique(obj.Freq))<numel(obj.Freq)))
                    for n= 1:numel(labels)
                        labels{n} = sprintf('Weights %s',num2str(n));
                        set(hlines(n),'tag',sprintf('Weights %s',num2str(n)));
                    end
                end
            end
                
            % This section is executed when the figure's previous legend
            % entries are going to be included in the new legend.  Previous
            % legend labels and line handles are appended to the newly
            % created ones in the LABELS and LINES cell arrays.
            if appendFlag
                % Get existing line handles and legend labels
                hExistingLines = varargin{1};
                existingLabels = varargin{2};

                % Append legend labels and line handles to any existing ones
                labels = [existingLabels(:); labels(:)];
                hlines = [hExistingLines ; hlines];
            end
            
            % If the current axes does not contain an existing legend and a
            % single frequency is being plotted, save the line handle and
            % label to AppData for possible later use.  Otherwise, create a
            % legend.
            if isNewSingleLinePlot
                 hLeg = legend(hlines,labels,'AutoUpdate','off');
                 set(hLeg,'Visible','off');
            else
                legend(hlines,labels,'location','BestOutside','AutoUpdate','off');
            end
        end
    end
    
    methods

        function varargout = plot(obj,varargin)
        %PLOT     Plot the response pattern
        %   plot(Hresp) plots the response pattern Hresp.
        %
        %   plot(Hresp, 'Units', UNIT) plots the response pattern using the
        %   unit specified in UNIT. UNIT can be any of the following:
        %   ['mag' | 'power' | {'db'}].
        %
        %   plot(..., 'NormalizeResp', NFLAG) plots the normalized response
        %   pattern if NFLAG is true. NFLAG can be any of the following:
        %   [true | {false}].
        %
        %   plot(..., 'Title', TITLE) uses TITLE as the title of the resulting
        %   plot.
        %
        %   plot returns an handle to the current plot object.
        %
        %   Example:
        %       % Construct a 2D isotropic array response pattern and then plot
        %       % it.
        %       hresp = phased.internal.RespPattern2D;
        %       plot(hresp)
            plotoption = getPlotOption(obj,varargin{:});
            plotoption.PlotType = 'plot';
            if isempty(plotoption.Title)
                plotoption.Title = genPlotTitle(obj);
            end
            
            [angles, response] = chkvectorinput(obj);
            
            response = privRespPattern(plotoption,response);
            
            % Get existing legend entries to be appended
            [hExistingLines,existingLabels] = obj.getLegendEntriesToAppend;
            
            if ishold && ~isempty(findall(gca,'type','surface'))
                error(message('phased:system:array:No2DPlotOver3D'));
            end
           
            h = plot(angles,response);
            annotatePlot(obj,plotoption);
            
            limitDynamicPlotRange(plotoption,h);
            
            % Create a legend and append existing legend entries to it if
            % needed.
            if isempty(hExistingLines)
                addLegend(obj,h);
            else
                % If appending previous legend entries, pass these
                % entries to ADDLEGEND.
                addLegend(obj,h,hExistingLines,existingLabels);
            end

            if nargout == 1,
                varargout{1} = h;
            end
        end
        
        function varargout = waterfall(obj,varargin)
        %WATERFALL  Plot response patterns in a mesh plot as functions of
        %           angle and frequency.
        %   mesh(Hresp) plots the response pattern Hresp.
        %
        %   waterfall(Hresp, 'Units', UNIT) plots the response pattern
        %   using the unit specified in UNIT. UNIT can be any of the
        %   following: ['mag' | 'power' | {'db'}].
        %
        %   waterfall(..., 'NormalizeResp', NFLAG) plots the normalized
        %   response pattern if NFLAG is true. NFLAG can be any of the
        %   following: [true | {false}].
        %
        %   waterfall(..., 'Title', TITLE) uses TITLE as the title of the
        %   resulting plot.
        %
        %   waterfall returns a handle to the current surfaceplot object if a
        %   vector of frequencies is being plotted.  If a single frequency
        %   is being plotted, waterfall return a handle to the current
        %   patch object.
        %
        %   Example:
        %       % Construct a 2D isotropic array response pattern and then plot
        %       % it.
        %       hresp = phased.internal.RespPattern2D;
        %       waterfall(hresp)
            plotoption = getPlotOption(obj,varargin{:});
            plotoption.PlotType = 'waterfall';
            if isempty(plotoption.Title)
                plotoption.Title = genPlotTitle(obj);
            end
            
            [angles, response] = chkvectorinput(obj);
            
            response = privRespPattern(plotoption,response);
            
            % Return a scaled frequency vector for plotting.  The frequency
            % vector is scaled to match an engineering unit (e.g. GHz or
            % MHz).
            freq = scalePlotFreqs(plotoption,obj.Freq);
            
            % If appending to an existing waterfall plot created with
            % plotResponse, append newly specified frequencies to the
            % existing frequencies.
            holdFlag = ishold;
            if holdFlag
                hSurf = findall(gca,'type',...
                    'surface','Tag','phasedWaterfall');
                if ~isempty(hSurf)
                    existingFreq = get(hSurf,'UserData');
                    freq = [existingFreq freq];
                elseif ~isempty(findall(gca,'type','line'))
                    error(message('phased:system:array:No3DPlotOver2D'));
                end
            end
            
            % Append to an existing plot if hold is on.  Otherwise Create a
            % mesh plot with Meshstyle set to 'Row'.  This generates a
            % waterfall style plot.
            if holdFlag
                set(hSurf,'XData',angles');
                set(hSurf,'YData',freq);
                existingZData = get(hSurf,'ZData');
                set(hSurf,'ZData',[existingZData ; response']);
            else
                hSurf = mesh(angles',freq,response');
                set(hSurf,'MeshStyle','Row')
            end
            
            % Tag the mesh surface and store frequencies so that it can be
            % appended to in future calls to plotResponse if hold is set to
            % on.
            set(hSurf,'Tag','phasedWaterfall');
            set(hSurf,'UserData',freq);
            
            % Set FaceAlpha below 1.0 to provide some transparency.  This
            % helps in visualizing partially obscured portions of a pattern
            % and avoids blocking  scene lighting from reaching curves
            % towards the back of the plot. Use less transparency (higher
            % FaceAlpha) for plots with more curves.  This helps reduce the
            % large amounts of visual clutter that would result in such
            % plots.
            numPatterns = length(freq);
            if numPatterns<=3
                set(hSurf,'FaceAlpha',0.2);
            elseif numPatterns<=5
                set(hSurf,'FaceAlpha',0.4);
            else
                set(hSurf,'FaceAlpha',0.9);
            end
            
            % Increase the line width to make the pattern more prominent.
            set(hSurf,'lineWidth',1.2)
            
            annotatePlot(obj,plotoption);
            
            limitDynamicPlotRange(plotoption,hSurf);
            
            % Remove ZData beyond Z-axis limits.  This ensures that nulls
            % are not drawn beyond axes boundaries.
            zLim = get(get(hSurf,'parent'),'ZLim');
            zData = get(hSurf,'Zdata');
            zData(zData<zLim(1))=NaN;
            set(hSurf,'Zdata',zData);

            if nargout == 1,
                varargout{1} = hSurf;
            end
        end
        
    end
    
    methods (Access = protected)
        function [angles, response] = chkvectorinput(obj)
            %return number of responses in each object
            Nobj = validateinputs(obj); 
            
            % Loop through each object and store angles and responses in
            % arrays. The third dimension corresponds to the object when
            % OBJ is an array of objects
            for cnt = 1:Nobj
                angles(:,1,cnt) = getPlotAngles(obj(cnt)); %#ok<AGROW>
                response(:,:,cnt) = obj(cnt).Pattern;  %#ok<AGROW>
            end
        end
        
    end
    
    methods
        function obj = set.Pattern(obj,value)
            validatePattern(obj,value);
            obj.Pattern = value;
        end
        function obj = set.Freq(obj,value)
            sigdatatypes.validateFrequency(value,...
                sprintf('%s.Freq',class(obj)),'Freq',{'vector'});
            obj.Freq = value;           
        end
    end
    
    methods (Static)
        function [hlines,labels] = getLegendEntriesToAppend
        %   GETLEGENDENTRIES returns the line handles and labels contained
        %   in the current axes' legend that need to be appended in the
        %   present plot.  T

        % Do not return legend entries if the current axes' hold state is
        % off.  In this case, existing legend entries will be erased.
            if ~ishold
                hlines = [];
                labels = [];
                return
            end
            
            % Get a handle to a legend associated with the current axes.
            drawnow expose;
            hLeg = get(gca,'Legend');
            
            if isempty(hLeg)
                hlines = [];
                labels = [];
            else
                
                hlines = get(hLeg,'PlotChildren');

                %Ensure that HLINES is a column vector.
                if size(hlines,2)>1
                    hlines = hlines';
                end

                labels = get(hLeg,'String');

            end
            
        end
    end
end

% Local utility functions

function Nobj = validateinputs(obj)
%VALIDATEINPUTS Validate inputs.
%   Validates inputs, and pre-allocates memory.

%calculate the number of patterns stored in obj and save in NRESP
Nobj = numel(obj);

Nresp = size(obj(1).Pattern,2);

dataLen = size(obj(1).Pattern,1); % Cache length of 1 data set to compare against.

% Validate inputs and error if necessary. If we support specifying multiple
% responses, we should verify that all the responses are the same length
% and each object contains the same number of responses. 
validateresponses(obj,Nobj,Nresp,dataLen);
end

function validateresponses(obj,Nobj,Nresp,dataLen) %#ok<INUSD>
%VALIDATERESPONSES Validate the array responses.
%   Verify that all the responses are the same length and each object
%   contains the same number of responses.

% ensure each object contains the same number of patterns and same pattern
% length.

% Currently we don't have the case for multiple responses.

% for k = 2:Nobj
%     if (size(obj(k).Pattern,1) ~= dataLen) || (size(obj(k).Pattern,2 ~= Nresp))
%         error(message('phased:phased:internal:RespPattern2D:InvalidResponseLength'));
%     end
% end

end

% [EOF]
