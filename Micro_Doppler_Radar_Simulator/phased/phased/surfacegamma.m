function gamma = surfacegamma(TerrainType,freq)
%surfacegamma Gamma value for different terrains
%   G = surfacegamma(TerrainType) returns the gamma value G (in dB) for the
%   terrain specified in TerrainType at 10 GHz.
%
%   The supported TerrainType are:
%       'sea state 3', 'sea state 5', 'woods', 'metropolitan',
%       'rugged mountain', 'farmland','wooded hill','flatland'
%
%   G = surfacegamma(TerrainType,FREQ) specifies frequency (in Hz) in FREQ
%   as a vector. The default value of FREQ is 10e9.
%
%   surfacegamma displays a table of above terrain types with their
%   corresponding gamma values. These gamma values are for an operating
%   frequency of 10 GHz.
%
%   % Example:
%   %   Determine the gamma value for woods area and then simulate the
%   %   clutter return from the area. Assume the radar system uses a single
%   %   cosine pattern antenna element and an operating frequency of 300
%   %   MHz.
%     
%   fc = 300e6;
%   g = surfacegamma('woods',fc);
%   hclutter = phased.ConstantGammaClutter('Gamma',g,...
%           'Sensor',phased.CosineAntennaElement,'OperatingFrequency',fc);
%   x = step(hclutter);
%   r = (0:numel(x)-1)/(2*hclutter.SampleRate)*hclutter.PropagationSpeed;
%   plot(r,abs(x)); xlabel('Range (m)'); ylabel('Clutter Magnitude (V)');
%   title('Clutter Return vs. Range');
%   
%   See also horizonrange, grazingang, phased.ConstantGammaClutter.

%   Copyright 2010-2011 The MathWorks, Inc.

%   Reference
%   [1] Fred Nathanson, Radar Design Principle, 1999
%   [2] Long, Radar Reflectivity of Land and Sea, 2001

%#codegen
%#ok<*EMCA>

    phased.internal.narginchk(0,2,nargin);

    TerrainTypeList = {...
        'Sea State 3',... %[2]
        'Sea State 5',... %[2]
        'Woods',... %[2]
        'Metropolitan',... %[2]
        'Rugged Mountain',... %[2]
        'Farmland',... %[1]
        'Wooded Hill',... %[1]
        'Flatland',... %[1]
        };

    if nargin > 0
        if nargin < 2
            freq = 10e9;
        end
        sigdatatypes.validateFrequency(freq,'surfacegamma','FREQ',...
            {'vector'});
        TerrainType = validatestring(TerrainType,TerrainTypeList,...
                                     'surfacegamma','TerrainType');        
        gamma = getTerrainGamma(TerrainType);
        gamma = gamma + 5*log10(freq./10e9);
    else
        cond = nargout == 0;
        if ~cond
            coder.internal.assert(cond,'MATLAB:nargoutchk:tooManyOutputs');
        end
        numgamma = numel(TerrainTypeList);
        for m = 1:numgamma
            fprintf('%s  %sdB\n',TerrainTypeList{m},num2str(getTerrainGamma(TerrainTypeList{m})));
        end
    end


function gamma = getTerrainGamma(TerrainType)
    

   % TerrainGamma =
   %  'Sea State 3'    -40
   %  'Sea State 5'    -30
   %  'Woods',         -15
   %  'Metropolitan',   0
   %  'Rugged Mountain' 0
   %  'Farmland',      -15
   %  'Wooded Hill',   -10
   %  'Flatland',     -20

       switch TerrainType(1)
         case 'S'
           if TerrainType(11) == '3'    
               gamma = -40;    %'Sea State 3'
           else
               gamma = -30;    %'Sea State 5'
           end
         case 'W'
           if TerrainType(5) == 's'
               gamma = -15;    %'Woods'
           else
               gamma = -10;    %'Wooded Hill'
           end
         case 'F'
           if TerrainType(2) == 'a'
               gamma = -15;    %'Farmland'
           else
               gamma = -20;    %'Flatland'
           end
         otherwise 
           gamma = 0;  %'Metropolitan','Rugged Mountain'
       end

