function [rng, ang] = rangeangle(pos, varargin)
%rangeangle Range and angle calculation
%   [RNG,ANG] = rangeangle(POS) returns the range RNG (in meters) and angle
%   ANG (in degrees) of the input position POS with respect to the origin.
%   POS must be a 3xN matrix with each column specifying a position in the
%   form of [x; y; z] (in meters) coordinates. 
%
%   [RNG,ANG] = rangeangle(POS,REFPOS) returns the ranges and angles of POS
%   with respect to the reference position specified in REFPOS. REFPOS is a
%   3x1 vector specifying the reference position in the form of [x; y; z]
%   (in meters) coordinates.
%
%   [RNG,ANG] = rangeangle(POS,REFPOS,REFAXES) returns the ranges and
%   angles of POS in the local coordinate system whose origin is REFPOS and
%   whose axes are defined in REFAXES. REFAXES must be a 3x3 matrix with
%   each column specifying the direction of an axis for the local
%   coordinate system in the form of [x; y; z] coordinates.
%
%   [...] = rangeangle(...,MODEL) specifies the model used when computing
%   the range and the angle between two positions as one of 'freespace' |
%   'two-ray', where the default is 'freespace'. 
%
%   When you specify MODEL as 'freespace', the computation only involves
%   direct path between the two positions. Under this setting, RNG is a 1xN
%   vector whose entries are the ranges for the corresponding positions
%   specified in POS. ANG is a 2xN matrix whose columns are the angles, in
%   the form of [azimuth; elevation], for the corresponding positions
%   specified in POS.
%
%   When you specify MODEL as 'two-ray', the computation involves both the
%   direct path and the ground reflection path. Under this setting, RNG is
%   a 1x2N vector whose adjacent entries are the ranges of the direct paths
%   and the ground reflection paths for the corresponding positions
%   specified in POS. ANG is a 2x2N matrix whose adjacent columns are the
%   angles of the direct paths and the ground reflection paths, in the form
%   of [azimuth; elevation], for the corresponding positions specified in
%   POS. Therefore, for the Nth position in POS, its ranges of the direct
%   path and the ground reflection path are in (2N-1)th and (2N)th entries
%   in RNG and its angles of the direct path and the ground reflection path
%   are in (2N-1)th and (2N)th columns.
%
%   Note that when the two-ray model is used, the coordinate system in
%   which the positions are specified is fixed. The xy-plane of the
%   coordinate system is assumed to be the flat earth and the z coordinate
%   represents the height above the ground.
%
%   % Examples:
%
%   % Example 1:
%   %   A target is located at [500; 0; 0] meters and a transmitter is 
%   %   located at [100; 100; 100] meters. The transmitter's local 
%   %   coordinates are given by [0 1 0;0 0 1;1 0 0]. Determine the range 
%   %   and angle of the target with respect to the transmitter expressed 
%   %   in the transmitter's local coordinate system.
%
%   posTgt = [500;0;0]; posTx = [100;100;100]; axTx = [0 1 0;0 0 1;1 0 0];
%   [tgt_rng,tgt_ang] = rangeangle(posTgt,posTx,axTx)
%
%   % Example 2:
%   %   The setting is the same as what described in Example 1 but now
%   %   compute the range and angle of the transmitter with respect to the
%   %   target. This information is necessary to compute the reflection 
%   %   when the target RCS is angle dependent. Assume the target's local
%   %   coordinate system is given by [1 0 0;0 1 0;0 0 1].
%
%   posTgt = [500;0;0]; posTx = [100;100;100]; axTgt = [1 0 0;0 0 1;1 0 0];
%   [tx_rng,tx_ang] = rangeangle(posTx,posTgt,axTgt)
%
%   % Example 3:
%   %   A target is located at [500; 0; 20] meters and a receiver is 
%   %   located at [100; 100; 100] meters. The receiver's local coordinates
%   %   are given by [0 1 0;0 0 1;1 0 0]. Determine the ranges and angles  
%   %   of the direct path and the ground reflection path of the target  
%   %   with respect to the receiver expressed in the receiver's local
%   %   coordinate system.
%
%   [tgt_rng,tgt_ang] = rangeangle([500; 0; 20],[100; 100; 100],...
%                           [0 1 0;0 0 1;1 0 0],'two-ray')
%
%   See also phased, global2localcoord, local2globalcoord.

%   Copyright 2010-2018 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>
    
phased.internal.narginchk(1,4,nargin);

if nargin > 1 && (ischar(varargin{end}) || isStringScalar(varargin{end}))
    model = varargin{end};
    Narg = nargin-2;
else
    model = 'freespace';
    Narg = nargin-1;
end

if Narg < 2
    refaxes = eye(3);
else
    refaxes = varargin{2};
end

if Narg < 1
    refpos = [0; 0; 0];
else
    refpos = varargin{1};
end

eml_assert_no_varsize(2:3, pos, refpos, refaxes);
sigdatatypes.validate3DCartCoord(pos,'rangeangle','Pos');
sigdatatypes.validate3DCartCoord(refpos,'rangeangle','RefPos',{'column'});
sigdatatypes.validate3DCartCoord(refaxes,'rangeangle','RefAxes',{'size',[3 3]});

model = validatestring(model,{'freespace','two-ray','tworay'},...
    'rangeangle','Model');

if model(1) == 'f'  % free space
    pos_in = pos;
else                % two-ray
    if any(pos(3,:)*refpos(3) < 0)
        error(message('phased:rangeangle:invalidHeight','Pos','RefPos'));
    end
    pos_in = reshape([pos;pos(1:2,:);-pos(3,:)],3,[]);
end
lclcoord = phased.internal.global2localcoord(pos_in,'rs',refpos,refaxes);
rng = lclcoord(3,:);
ang = lclcoord(1:2,:);


% [EOF]
