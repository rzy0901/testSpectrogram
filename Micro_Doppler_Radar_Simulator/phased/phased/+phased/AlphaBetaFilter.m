classdef (Sealed) AlphaBetaFilter < matlabshared.tracking.internal.AbstractAlphaBetaFilter
%AlphaBetaFilter     Alpha-beta filter
%   H = phased.AlphaBetaFilter returns an alpha-beta filter object, H. This
%   object performs alpha-beta filer based tracking on measurements.
%
%   H = phased.AlphaBetaFilter(Name,Value) creates an alpha-beta filter
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   AlphaBetaFilter methods:
%
%   predict    - Perform state prediction using alpha-beta filter
%   correct    - Perform state correction using alpha-beta filter
%   distance   - Distances between measurements and filter prediction
%   residual   - Calculate the measurement residual and residual noise
%   clone      - Create alpha-beta filter object with same property
%                values
%
%   AlphaBetaFilter properties:
%
%   MotionModel      - Model of target motion
%   State            - Alpha-beta filter states
%   StateCovariance  - State estimation error covariance
%   ProcessNoise     - Process noise covariance 
%   MeasurementNoise - Measurement noise covariance 
%   Coefficients     - Alpha-beta filter coefficients
%
%   For more tracking filters and data association algorithms, please refer
%   to Sensor Fusion and Tracking Toolbox.
%
%   % Examples:
%
%   % Example 1:
%   %   Apply the alpha-beta filter to a moving target along x axis.
%
%   T = 0.1;  V0 = 100;  N = 100;
%   plat = phased.Platform('MotionModel','Velocity',...
%       'VelocitySource','Input port','InitialPosition',[100;0;0]);
%   abfilt = phased.AlphaBetaFilter('MotionModel','1D Constant Velocity');
%   Z = zeros(1,N);  Zp = zeros(1,N); Zc = zeros(1,N);
%   for m = 1:N
%       pos = plat(T,[100+20*randn;0;0]);
%       Z(m) = pos(1);
%       [~,~,Zp(m)] = predict(abfilt,T);
%       [~,~,Zc(m)] = correct(abfilt,Z(m));
%   end
%   t = (0:N-1)*T; plot(t,Z,t,Zp,t,Zc);
%   xlabel('Time (s)'); ylabel('Position (m)');
%   legend('True Track','Predicted Track','Corrected Track',...
%       'Location','Best');
%
%   % Example 2:
%   %   Apply the alpha-beta filter to an accelerating target along x axis.
%
%   T = 0.1;  a0 = 100;  N = 100;
%   plat = phased.Platform('MotionModel','Acceleration',...
%       'AccelerationSource','Input port','InitialPosition',[100;0;0]);
%   abfilt = phased.AlphaBetaFilter(...
%       'MotionModel','1D Constant Acceleration',...
%       'Coefficients',[0.5 0.5 0.1]);
%   Z = zeros(1,N);  Zp = zeros(1,N); Zc = zeros(1,N);
%   for m = 1:N
%       pos = plat(T,[100+20*randn;0;0]);
%       Z(m) = pos(1);
%       [~,~,Zp(m)] = predict(abfilt,T);
%       [~,~,Zc(m)] = correct(abfilt,Z(m));
%   end
%   t = (0:N-1)*T; plot(t,Z,t,Zp,t,Zc);
%   xlabel('Time (s)'); ylabel('Position (m)');
%   legend('True Track','Predicted Track','Corrected Track',...
%       'Location','Best');
%
%   % Example 3:
%   %   Apply the alpha-beta filter to a moving target in 3D space.
%
%   T = 0.1;  V0 = 100;  N = 100;
%   plat = phased.Platform('MotionModel','Velocity',...
%       'VelocitySource','Input port','InitialPosition',[100;0;0]);
%   abfilt = phased.AlphaBetaFilter('MotionModel',...
%       '3D Constant Velocity','State',zeros(6,1));
%   Z = zeros(3,N);  Zp = zeros(3,N); Zc = zeros(3,N);
%   for m = 1:N
%       Z(:,m) = plat(T,[V0+20*randn;0;0]);
%       [~,~,Zp(:,m)] = predict(abfilt,T);
%       [~,~,Zc(:,m)] = correct(abfilt,Z(:,m));
%   end
%   t = (0:N-1)*T; plot(t,Z(1,:),t,Zp(1,:),t,Zc(1,:));
%   xlabel('Time (s)'); ylabel('Position along X (m)');
%   legend('True Track','Predicted Track','Corrected Track',...
%       'Location','Best');
%
%   See also phased, phased.CFARDetector.

%   Copyright 2015-2018 The MathWorks, Inc.

%   Reference
%   [1] Blackman, Multiple-Target Tracking with Radar Applications, Artech
%   House, 1986
%   [2] Bar-Shalom et al., Estimation with Applications to Tracking and
%   Navigation, Theory, Algorithms, and Software. John Wiley & Sons, 2001


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods
        function obj = AlphaBetaFilter(varargin)
            % Parse the inputs.
            obj@matlabshared.tracking.internal.AbstractAlphaBetaFilter(...
                '1D Constant Velocity',varargin{:});
        end

    end
    
    methods

        function [x_pred,P_pred,z_pred] = predict(obj,T)
        %predict    Predicts the state and measurement
        %   X_pred = predict(H,T) returns the predicted state, X_pred, in a
        %   column vector using alpha-beta filter. T is scalar specifying
        %   the elapsed time (in seconds) between the current prediction
        %   and last prediction/correction.
        %
        %   [X_pred,P_pred] = predict(H,T) also returns the state estimate
        %   covariance in a matrix, P_pred.
        %
        %   [X_pred,P_pred,Z_pred] = predict(H,T) also returns the
        %   predicted measurement in a column vector, Z_pred.
        %
        %   % Example:
        %   %   Apply the alpha-beta filter to a moving target along x 
        %   %   axis. Note that in the second half of the simulation,
        %   %   measurements are not available so only predict() is used to
        %   %   track the target.
        %
        %   T = 0.1;  V0 = 100;  N = 100;
        %   plat = phased.Platform('MotionModel','Velocity',...
        %       'VelocitySource','Input port','InitialPosition',[100;0;0]);
        %   abfilt = phased.AlphaBetaFilter(...
        %       'MotionModel','1D Constant Velocity');
        %   Z = zeros(1,N);  Zp = zeros(1,N); Zc = zeros(1,N);
        %   for m = 1:N
        %       pos = step(plat,T,[100+20*randn;0;0]);
        %       Z(m) = pos(1);
        %       [~,~,Zp(m)] = predict(abfilt,T);
        %       if m <= N/2
        %           [~,~,Zc(m)] = correct(abfilt,Z(m));
        %       else
        %           Zc(m) = Zp(m);
        %       end
        %   end
        %   t = (0:N-1)*T; plot(t,Z,t,Zp,t,Zc);
        %   xlabel('Time (s)'); ylabel('Position (m)');
        %   legend('True Track','Predicted Track','Corrected Track',...
        %       'Location','Best');
        %
        %   See also phased, phased.AlphaBetaFilter/correct
        
            [x_pred,P_pred,z_pred] = predict@matlabshared.tracking.internal.AbstractAlphaBetaFilter(obj,T);
        end
        
        function [x_corr,P_corr,z_corr] = correct(obj,z,varargin)
        % correct    Corrects the state and measurement
        %   X_corr = correct(H,Z) returns the corrected state, X_corr, in a
        %   column vector using alpha-beta filter based on the measurement
        %   vector, Z. Z is also a column vector.
        %
        %   [X_corr,P_corr] = correct(H,Z) also returns the corrected state
        %   covariance matrix, P_corr.
        %
        %   [X_corr,P_corr,Z_corr] = correct(H,Z) also returns the
        %   corrected measurement in a column vector, Z_corr.
        %
        %   % Example:
        %   %   Apply the alpha-beta filter to a moving target along x 
        %   %   axis.
        %
        %   T = 0.1;  V0 = 100;  N = 100;
        %   plat = phased.Platform('MotionModel','Velocity',...
        %       'VelocitySource','Input port','InitialPosition',[100;0;0]);
        %   abfilt = phased.AlphaBetaFilter(...
        %       'MotionModel','1D Constant Velocity');
        %   Z = zeros(1,N);  Zp = zeros(1,N); Zc = zeros(1,N);
        %   for m = 1:N
        %       pos = step(plat,T,[100+20*randn;0;0]);
        %       Z(m) = pos(1);
        %       [~,~,Zp(m)] = predict(abfilt,T);
        %       [~,~,Zc(m)] = correct(abfilt,Z(m));
        %   end
        %   t = (0:N-1)*T; plot(t,Z,t,Zp,t,Zc);
        %   xlabel('Time (s)'); ylabel('Position (m)');
        %   legend('True Track','Predicted Track','Corrected Track',...
        %       'Location','Best');
        %
        %   See also phased, phased.AlphaBetaFilter/predict
        
            [x_corr,P_corr,z_corr] = correct@matlabshared.tracking.internal.AbstractAlphaBetaFilter(obj,z,varargin{:});
            
        end
    end
end