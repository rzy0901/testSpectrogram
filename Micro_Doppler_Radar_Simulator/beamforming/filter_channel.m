% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     filter_channel.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       December 2015 created
%
%  *************************************************************************************
%    Description:
%
%    Performs rotation of angles to coordination system associated with
%    Antenna Pattern
%
%    [ ch ] = filter_channel(ant_type, ch)
%
%    Inputs:
%
%       1. ant_type - type of antenna chose by user
%       2. ch       - channel structure
%       3. varargin - angles of Tx, Rx rotation in azimuth plane along LOS
%       direction for PAA only
%
%    Outputs:
%
%       1. ch - modified channel structure
%
%  *************************************************************************************/
function [ ch ] = filter_channel(ant_type, ch, varargin)

switch(ant_type)
    case 1, % phase antenna array with single polarization
        % Rotation of Tx, Rx for NLOS case
        phiTx = varargin{1};
        phiRx = varargin{2};
        
        [ch.tx_az, ch.tx_el] = basic2rot(ch.tx_az, ch.tx_el, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az, ch.rx_el] = basic2rot(ch.rx_az, ch.rx_el, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az > 90 & ch.tx_az < 270);
        ind2 = find(ch.rx_az > 90 & ch.rx_az < 270);
        ind = union(ind1, ind2);
        
        ch.am(ind) = 0;
    case 2, % phased antenna array with dual polarization
        % Rotation of Tx, Rx for NLOS case
        phiTx = varargin{1};
        phiRx = varargin{2};
        
        [ch.tx_az, ch.tx_el] = basic2rot(ch.tx_az, ch.tx_el, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az, ch.rx_el] = basic2rot(ch.rx_az, ch.rx_el, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az > 90 & ch.tx_az < 270);
        ind2 = find(ch.rx_az > 90 & ch.rx_az < 270);
        ind = union(ind1, ind2);
        
        ch.am_h11(ind) = 0;
        ch.am_h12(ind) = 0;
        ch.am_h21(ind) = 0;
        ch.am_h22(ind) = 0;
    case 3, % double phased antenna array
        % Rotation of Tx, Rx for NLOS case
        phiTx = varargin{1};
        phiRx = varargin{2};
        
        shift_tx11 = varargin{3};
        shift_rx11 = varargin{4};
        shift_tx22 = varargin{5};
        shift_rx22 = varargin{6};
        
        [ch.tx_az_11, ch.tx_el_11] = basic2rot(ch.tx_az_11, ch.tx_el_11, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_11, ch.rx_el_11] = basic2rot(ch.rx_az_11, ch.rx_el_11, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_11 > 90 + shift_tx11 & ch.tx_az_11 < 270 + shift_tx11);
        ind2 = find(ch.rx_az_11 > 90 + shift_rx11 & ch.rx_az_11 < 270 + shift_rx11);
        ind = union(ind1, ind2);
        
        ch.am_11(ind) = 0;

        [ch.tx_az_12, ch.tx_el_12] = basic2rot(ch.tx_az_12, ch.tx_el_12, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_12, ch.rx_el_12] = basic2rot(ch.rx_az_12, ch.rx_el_12, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_12 > 90 + shift_tx11 & ch.tx_az_12 < 270 + shift_tx11);
        ind2 = find(ch.rx_az_12 > 90 + shift_rx22 & ch.rx_az_12 < 270 + shift_rx22);
        ind = union(ind1, ind2);
        
        ch.am_12(ind) = 0;
        
        [ch.tx_az_21, ch.tx_el_21] = basic2rot(ch.tx_az_21, ch.tx_el_21, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_21, ch.rx_el_21] = basic2rot(ch.rx_az_21, ch.rx_el_21, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_21 > 90 + shift_tx22 & ch.tx_az_21 < 270 + shift_tx22);
        ind2 = find(ch.rx_az_21 > 90 + shift_rx11 & ch.rx_az_21 < 270 + shift_rx11);
        ind = union(ind1, ind2);
        
        ch.am_21(ind) = 0;
        
        [ch.tx_az_22, ch.tx_el_22] = basic2rot(ch.tx_az_22, ch.tx_el_22, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_22, ch.rx_el_22] = basic2rot(ch.rx_az_22, ch.rx_el_22, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_22 > 90 + shift_tx22 & ch.tx_az_22 < 270 + shift_tx22);
        ind2 = find(ch.rx_az_22 > 90 + shift_rx22 & ch.rx_az_22 < 270 + shift_rx22);
        ind = union(ind1, ind2);
        
        ch.am_22(ind) = 0;
    case 4, % single phased antenna array on Tx, Rx; Rx receives w/ dual polarization
        % Rotation of Tx, Rx for NLOS case
        phiTx = varargin{1};
        phiRx = varargin{2};
        shift_tx11 = varargin{3};
        shift_rx11 = varargin{4};
        shift_tx22 = varargin{5};
        shift_rx22 = varargin{6};
        
        [ch.tx_az_11, ch.tx_el_11] = basic2rot(ch.tx_az_11, ch.tx_el_11, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_11, ch.rx_el_11] = basic2rot(ch.rx_az_11, ch.rx_el_11, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_11 > 90 + shift_tx11 & ch.tx_az_11 < 270 + shift_tx11);
        ind2 = find(ch.rx_az_11 > 90 + shift_rx11 & ch.rx_az_11 < 270 + shift_rx11);
        ind = union(ind1, ind2);
        
        ch.am_11vv(ind) = 0;
        ch.am_11vh(ind) = 0;
        ch.am_11hv(ind) = 0;
        ch.am_11hh(ind) = 0;

        [ch.tx_az_12, ch.tx_el_12] = basic2rot(ch.tx_az_12, ch.tx_el_12, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_12, ch.rx_el_12] = basic2rot(ch.rx_az_12, ch.rx_el_12, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_12 > 90 + shift_tx11 & ch.tx_az_12 < 270 + shift_tx11);
        ind2 = find(ch.rx_az_12 > 90 + shift_rx22 & ch.rx_az_12 < 270 + shift_rx22);
        ind = union(ind1, ind2);
        
        ch.am_12vv(ind) = 0;
        ch.am_12vh(ind) = 0;
        ch.am_12hv(ind) = 0;
        ch.am_12hh(ind) = 0;
        
        [ch.tx_az_21, ch.tx_el_21] = basic2rot(ch.tx_az_21, ch.tx_el_21, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_21, ch.rx_el_21] = basic2rot(ch.rx_az_21, ch.rx_el_21, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_21 > 90 + shift_tx22 & ch.tx_az_21 < 270 + shift_tx22);
        ind2 = find(ch.rx_az_21 > 90 + shift_rx11 & ch.rx_az_21 < 270 + shift_rx11);
        ind = union(ind1, ind2);
        
        ch.am_21vv(ind) = 0;
        ch.am_21vh(ind) = 0;
        ch.am_21hv(ind) = 0;
        ch.am_21hh(ind) = 0;
        
        [ch.tx_az_22, ch.tx_el_22] = basic2rot(ch.tx_az_22, ch.tx_el_22, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az_22, ch.rx_el_22] = basic2rot(ch.rx_az_22, ch.rx_el_22, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az_22 > 90 + shift_tx22 & ch.tx_az_22 < 270 + shift_tx22);
        ind2 = find(ch.rx_az_22 > 90 + shift_rx22 & ch.rx_az_22 < 270 + shift_rx22);
        ind = union(ind1, ind2);
        
        ch.am_22vv(ind) = 0;
        ch.am_22vh(ind) = 0;
        ch.am_22hv(ind) = 0;
        ch.am_22hh(ind) = 0;
        
    case 5,
        % Rotation of Tx, Rx for NLOS case
        phiTx = varargin{1};
        phiRx = varargin{2};
        
        [ch.tx_az, ch.tx_el] = basic2rot(ch.tx_az, ch.tx_el, mod(phiTx + 90, 360), 90, 0);
        [ch.rx_az, ch.rx_el] = basic2rot(ch.rx_az, ch.rx_el, mod(phiRx + 90, 360), 90, 0);
        
        % Check boundaries for AoA, AoD and delete not existing rays
        ind1 = find(ch.tx_az > 90 & ch.tx_az < 270);
        ind2 = find(ch.rx_az > 90 & ch.rx_az < 270);
        ind = union(ind1, ind2);
        
        ch.am_h11(ind) = 0;
        ch.am_h12(ind) = 0;
end
