function stapgpu(hcpugpuaction,hcurrenttime,hcluttersimtimecpu,...
    hcluttersimtimemex,hcluttersimtimegpu,hcluttersimtimemexspeedup,...
    hcluttersimtimegpuspeedup,receiversignal_line,receiverangdop_image,...
    hprocessedsignal_line,hweightsangdop_image,hwaitbar_text,...
    hwaitbar_fill,clutter_validate)
% This function stapgpu is only in support of STAPCPUGPUExample. It may be
% removed in a future release.

%   Copyright 2012-2016 The MathWorks, Inc.

load BasicMonostaticRadarExampleData.mat collector ...
    radiator receiver transmitter waveform;   % Load monostatic pulse radar

fc = radiator.OperatingFrequency; 
c = radiator.PropagationSpeed;
lambda = c/fc;
numpulse = 10; % Number of pulses

antenna = phased.IsotropicAntennaElement('FrequencyRange',[8e9 12e9]);  
ula = phased.ULA('Element',antenna,'NumElements',6,'ElementSpacing', lambda/2);

radiator.Sensor = ula;
collector.Sensor = ula;
sensormotion = phased.Platform('InitialPosition',[0; 0; 1000]);
arrayAxis = [0; 1; 0];
prf = waveform.PRF;
vr = ula.ElementSpacing*prf;  % in [m/s]
sensormotion.Velocity = vr/2*arrayAxis;

% Target
target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',1, ...
    'OperatingFrequency', fc);
tgtmotion = phased.Platform('InitialPosition',[1000; 1000; 0],...
    'Velocity',[30; 30; 0]);

% Jammer
jammer = phased.BarrageJammer('ERP',100);
Fs = waveform.SampleRate;
rngbin = c/2*(0:1/Fs:1/prf-1/Fs).';
jammer.SamplesPerFrame = numel(rngbin);
jammermotion = phased.Platform(...
    'InitialPosition',[1000; 1732; 1000]);

% Clutter
clutter = phased.gpu.ConstantGammaClutter('Sensor',ula,...
    'SampleRate',Fs,...
    'Gamma',-15,'PlatformHeight',1000,...
    'OperatingFrequency',fc,...
    'PropagationSpeed',c,...
    'PRF',prf,...
    'TransmitERP',transmitter.PeakPower*db2pow(transmitter.Gain),...
    'PlatformSpeed',norm(sensormotion.Velocity),...
    'PlatformDirection',[90;0],...
    'BroadsideDepressionAngle',0,...
    'MaximumRange',clutter_validate.MaximumRange,...
    'AzimuthCoverage',180,...
    'PatchAzimuthWidth',clutter_validate.PatchAzimuthWidth,...
    'OutputFormat','Pulses');


% Propagation path
tgtchannel = phased.FreeSpace('TwoWayPropagation',true,...
    'SampleRate',Fs,'OperatingFrequency', fc); 
jammerchannel = phased.FreeSpace('TwoWayPropagation',false,...
    'SampleRate',Fs,'OperatingFrequency', fc); 

% Angle-Doppler response
angdopresp_signal = phased.AngleDopplerResponse('SensorArray',ula,...
              'OperatingFrequency',fc, 'PropagationSpeed',c,...
              'PRF',prf, 'ElevationAngle',0);
angdopresp_weights = phased.AngleDopplerResponse('SensorArray',ula,...
              'OperatingFrequency',fc, 'PropagationSpeed',c,...
              'PRF',prf, 'ElevationAngle',0);

% SMI
stap = phased.STAPSMIBeamformer('SensorArray', ula, 'PRF', prf, ...
    'PropagationSpeed', c, 'OperatingFrequency', fc, ...
    'DirectionSource', 'Input port', 'DopplerSource', 'Input port', ...
    'WeightsOutputPort', true,...
    'NumGuardCells', 4, 'NumTrainingCells', 100);


tsig = zeros(size(rngbin,1), ula.NumElements, numpulse);
jsig = tsig; tjcsig = tsig; csig = tsig; 

total_iter = 5*numpulse;
        
for m = 1:total_iter
    if getappdata(hcpugpuaction,'cpugpusim_stop')
        set(hcpugpuaction,'Enable','off');
        updatewaitbar(hwaitbar_text,hwaitbar_fill,inf);
        break;
    end
    storeidx = rem(m-1,numpulse)+1;
    
    % Update sensor, target, and jammer positions
    [sensorpos,sensorvel] = sensormotion(1/prf);
    [tgtpos,tgtvel] = tgtmotion(1/prf);
    [jampos,jamvel] = jammermotion(1/prf);           

    % Calculate the target and jammer angles as seen by the radar
    [~,tgtang] = rangeangle(tgtpos,sensorpos);       
    [~,jamang] = rangeangle(jampos,sensorpos);    

    % Simulate propagation of pulse in direction of target
    pulse = waveform();
    [pulse,txstatus] = transmitter(pulse);
    pulse = radiator(pulse,tgtang);
    pulse = tgtchannel(pulse,sensorpos,tgtpos,sensorvel,tgtvel);
    
    % Collect target returns at sensor
    pulse = target(pulse);
    tsig(:,:,storeidx) = collector(pulse,tgtang);
    
    % Collect jammer signal at sensor
    jamsig = jammer();
    jamsig = jammerchannel(jamsig,jampos,sensorpos,jamvel,sensorvel);
    jsig(:,:,storeidx) = collector(jamsig,jamang);
    
    % Collect clutter signal at sensor
    tcsim = tic;
    csig(:,:,storeidx) = clutter();
    csimtime = toc(tcsim);
    set(hcurrenttime,'String',...
        num2str(round(csimtime*1e3),'%d')); drawnow;

    tcpu = str2double(get(hcluttersimtimecpu,'String'));
    tmex = str2double(get(hcluttersimtimemex,'String'));
    tgpu = str2double(get(hcluttersimtimegpu,'String'));
    if ~isnan(tcpu) && ~isnan(tgpu)
        set(hcluttersimtimegpuspeedup,'String',...
            num2str(round(tcpu/tgpu),'%d')); drawnow;
    end
    if ~isnan(tcpu) && ~isnan(tmex)
        set(hcluttersimtimemexspeedup,'String',...
            num2str(round(tcpu/tmex),'%d')); drawnow;
    end

    % Receive collected signals
    tjcsig(:,:,storeidx) = receiver(...
        tsig(:,:,storeidx)+jsig(:,:,storeidx)+csig(:,:,storeidx),...
        ~(txstatus>0)); % Target + jammer + clutter

    set(receiversignal_line,'YData',abs(tjcsig(:,1,storeidx))); drawnow;


    if ~rem(m,10)
        % *True Target Range, Angle and Doppler*
        % 
        % The target azimuth angle is 45 degrees, and the elevation angle
        % is about -35.27 degrees.
        tgtLocation = global2localcoord(tgtpos,'rs',sensorpos);
        tgtAzAngle = tgtLocation(1);
        tgtElAngle = tgtLocation(2);
        tgtRng = tgtLocation(3);

        % The target Doppler normalized frequency is about 0.21.
        sp = radialspeed(tgtpos, tgtmotion.Velocity, ...
                        sensorpos, sensormotion.Velocity);
        tgtDp = 2*speed2dop(sp,lambda);  % Round trip Doppler     
        tgtCellIdx = val2ind(tgtRng,c/(2*Fs));

        % Plot angle Doppler response of the target snapshot
        [sigad_resp,sigad_ang_grid,sigad_dop_grid] = ...
            angdopresp_signal(shiftdim(tjcsig(tgtCellIdx,:,:)));
        set(receiverangdop_image,'XData',sigad_ang_grid,...
            'YData',sigad_dop_grid/prf,...
            'CData',mag2db(abs(sigad_resp))); drawnow;

        tgtAngle = [tgtAzAngle; tgtElAngle];
        [y,w] = stap(tjcsig,tgtCellIdx,tgtAngle,tgtDp);

        set(hprocessedsignal_line,'YData',abs(y)); drawnow;

        % Plot angle Doppler response of the STAP weights
        [w_resp,w_ang_grid,w_dop_grid] = angdopresp_weights(w);
        set(hweightsangdop_image,'XData',w_ang_grid,...
            'YData',w_dop_grid/prf,...
            'CData',mag2db(abs(w_resp))); drawnow;

    end
    updatewaitbar(hwaitbar_text,hwaitbar_fill,m/total_iter);
end
end

function updatewaitbar(hwaitbar_text,hwaitbar_fill,percentdone)
    if isinf(percentdone)
        set(hwaitbar_text,'String','Stopped');
    else
        set(hwaitbar_fill,'Visible','on',...
            'Position',[0 0 percentdone 1]);
        set(hwaitbar_text,'String',...
            sprintf('%4.1f%%',percentdone*100));
    end
end
