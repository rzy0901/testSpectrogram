% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_gen_intra_cls.m
%    Authors:       A. Lomayev, R. Maslennikov
%    Version:       1.0
%    History:       May 2010 created
%
%  *************************************************************************************
%    Description:
% 
%    function returns intra-cluster space-temporal parameters: rays amplitudes,
%    times of arrival in [ns] relative to cluster time, azimuth/elevation angles
%    in [deg] relative to cluster direction for TX/RX in CR environment
% 
%    [incls] = cr_gen_intra_cls(toa_cls)
%
%    Inputs:
%
%       1. toa_cls - cluster time of arrival 
%
%    Outputs:
%
%       1. incls.am    - amplitudes array
%       2. incls.toa   - times of arrival array
%       3. incls.tx_az - TX azimuths array
%       4. incls.tx_el - TX elevations array
%       5. incls.rx_az - RX azimuths array
%       6. incls.rx_el - RX elevations array
%
%  *************************************************************************************/
function [incls] = cr_gen_intra_cls(toa_cls)

% intra-clusters parameters
K_f_dB = 10;
K_b_dB = 14.2;

N_f = 6;
N_b = 8;

lambda_f = 0.37;
lambda_b = 0.31;

gamma_f  = 3.7;
gamma_b  = 4.5;

sigma_tx_az = 5;
sigma_tx_el = 5;
sigma_rx_az = 5;
sigma_rx_el = 5;

% generate time of arrivals
toa_f = [];
if (N_f > 0)
    toa_f = sort(-poisson(lambda_f, N_f));
end

% check if the absolute ray time is negative 
% which is very rare but probable event

ix = find((toa_f + toa_cls) < 0);
if (~isempty(ix))
    toa_f(ix) = [];
end

toa_b = [];
if (N_b > 0)
    toa_b = poisson(lambda_b, N_b);
end

% generate amplitudes
K_f = 10^(K_f_dB/10); % Convert from [dB]
K_b = 10^(K_b_dB/10); % Convert from [dB]


ave_pow_f = exp(toa_f/gamma_f)/K_f;
ave_pow_b = exp(-toa_b/gamma_b)/K_b;

toa = [toa_f; 0; toa_b];

amp_f = rayleigh(ave_pow_f);
amp_b = rayleigh(ave_pow_b);

amp = [amp_f; 1; amp_b];

phase = rand(size(amp)); % Random phase
imp_resp = amp .* exp(2*pi*1i*phase);

% normalize power
imp_resp = imp_resp./norm(imp_resp);

% generate angles
tx_az_f = randn(N_f, 1) * sigma_tx_az;
tx_az_f(ix) = [];
tx_az_b = randn(N_b, 1) * sigma_tx_az;
tx_az   = [tx_az_f; 0; tx_az_b];

tx_el_f = randn(N_f, 1) * sigma_tx_el;
tx_el_f(ix) = [];
tx_el_b = randn(N_b, 1) * sigma_tx_el;
tx_el   = [tx_el_f; 0; tx_el_b];

rx_az_f = randn(N_f, 1) * sigma_rx_az;
rx_az_f(ix) = [];
rx_az_b = randn(N_b, 1) * sigma_rx_az;
rx_az   = [rx_az_f; 0; rx_az_b];

rx_el_f = randn(N_f, 1) * sigma_rx_el;
rx_el_f(ix) = [];
rx_el_b = randn(N_b, 1) * sigma_rx_el;
rx_el   = [rx_el_f; 0; rx_el_b];

% fill intra-cluster output structure
incls.am    = imp_resp;
incls.toa   = toa;
incls.tx_az = tx_az;
incls.tx_el = tx_el;
incls.rx_az = rx_az;
incls.rx_el = rx_el;

% ------------------------------------------------------------------------ % 

function y = poisson(lambda, numSamples)
% generates poisson process
% lambda - rate
% numSmaples - number of elements
sigma = sqrt(1/lambda/2);
x = (sigma*randn(numSamples,1)).^2 + (sigma*randn(numSamples,1)).^2;
y = cumsum(x, 1);

% ------------------------------------------------------------------------ % 

function y = rayleigh(p)
% generates rayleigh amplitudes
% p - average power
y = sqrt(p).*sqrt(abs(randn(size(p))+1i*randn(size(p))).^2/2);