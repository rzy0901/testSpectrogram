classdef (Hidden) AbstractRespPatternPlotOption < sigutils.pvpairs
%This class is for internal use only. It may be removed in the future.

%AbstractRespPatternPlotOption   Hidden option class for plotting
%response patterns

%   Copyright 2009-2010 The MathWorks, Inc.

    properties 
        NormalizeResp = false;
        Units = 'db';
        Title;
    end
    
    methods
        function response = privRespPattern(obj,response)
            %PRIVRESPPATTERN Prepare response pattern
            %   PRIVRESPPATTERN is a protected method which prepares the
            %   response pattern object for plotting.
            
            response = phased.internal.computePlotPattern(...
                response,obj.NormalizeResp,obj.Units);
            
        end
        
        function resplbl = getRespLabel(obj)
            %getRespLabel Obtain response pattern label
            
            resplblPrefix = '';           
            
            if obj.NormalizeResp
                resplblPrefix = 'Normalized';
            end
            
            resplbl = sprintf('%s %s',resplblPrefix, getUnitLabel(obj));
            
        end

        function hax = limitDynamicPlotRange(this,hplotobj)
            hax = get(hplotobj(1),'Parent');
            if strcmpi(this.Units,'db')
                datalim = getDataLimit(this,hplotobj);
                if datalim(1) == datalim(2)
                    if any(isinf(datalim))
                        % set to -1 to 1 if all data is inf
                        rlim(1) = -1;
                        rlim(2) = 1;
                    else
                        rlim(1) = min(datalim(1)-1,0.9*datalim(1));
                        rlim(2) = max(datalim(2)+1,1.1*datalim(2));
                    end
                else
                    rlim(1) = max(datalim(1),datalim(2)-100); % min 
                    rlim(2) = datalim(2)+0.1*(datalim(2)-rlim(1)); %max
                end
                limlabel = getLimLabel(this);
                set(hax,limlabel,rlim);
            end
        end
        
    end
    
    methods (Access = protected)
        function unitlbl = getUnitLabel(obj)
            switch obj.Units
                case 'mag',
                    unitlbl = 'Magnitude';
                case 'power',
                    unitlbl = 'Power';
                case 'db'
                    unitlbl = 'Power (dB)';
            end
        end
    end
    
    methods (Access = protected, Abstract)
        limlabel = getLimLabel(this);
        datalim = getDataLimit(this,hplotobj);
    end
    
    methods
        
        function set.NormalizeResp(this,val)
            validateattributes(val,{'logical'},{'scalar'},...
                sprintf('%s.NormalizeResp',class(this)),'NormalizeResp');
            this.NormalizeResp = val;
        end
        
        function set.Units(this,val)
            val = validatestring(val,{'mag','power','db'},...
                sprintf('%s.Units',class(this)),'Units');
            this.Units = val;
        end
        
    end

end

% [EOF]


