% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     normalize.m
%    Authors:       A. Lomayev, R. Maslennikov, Y. Gagiev
%    Version:       1.0
%    History:       May 2010 created
%                   April 2016 updated
%
%  *************************************************************************************
%    Description:
% 
%    function returns normalized channel responses
%
%    [imp_res] = normalize(ant_type, Pnorm, imp_res)
%
%    Inputs:
%
%    1. ant_type    - antenna type that defines internal structure of imp_res
%    2. Pnorm       - parameter defining apply normalization
%    3. imp_res     - structure with impulse responses after beamforming
% 
%    Outputs:
%
%    1. imp_res - structure with digitized impulse responses
%
%    Update: Normalization made in accordance with Evaluation methodology
%  *************************************************************************************/
function [ imp_res ] = normalize(ant_type, Pnorm, imp_res)

% Normalization is done using Frobenius norm taking into account number of
% Rx chains as defined in Evaluation methodology for 11ay standard
if ( Pnorm )
switch(ant_type)
    case num2cell(1:3), % conf#1, conf#2, conf#3
         Nrx = 2; % number of Rx chains
         H11 = sum(abs(imp_res.h11).^2);
         H12 = sum(abs(imp_res.h12).^2);
         H21 = sum(abs(imp_res.h21).^2);
         H22 = sum(abs(imp_res.h22).^2);
         normH = sqrt((H11 + H12 + H21 + H22)/Nrx);
         imp_res.h11 = imp_res.h11 ./ normH;
         imp_res.h12 = imp_res.h12 ./ normH;
         imp_res.h21 = imp_res.h21 ./ normH;
         imp_res.h22 = imp_res.h22 ./ normH;
         
    case 4, % conf#4
         Nrx = 4; % number of Rx chains
         H11 = sum(abs(imp_res.h11).^2);
         H12 = sum(abs(imp_res.h12).^2);
         H21 = sum(abs(imp_res.h21).^2);
         H22 = sum(abs(imp_res.h22).^2);

         H33 = sum(abs(imp_res.h33).^2);
         H34 = sum(abs(imp_res.h34).^2);
         H43 = sum(abs(imp_res.h43).^2);
         H44 = sum(abs(imp_res.h44).^2);
         
         H31 = sum(abs(imp_res.h31).^2);
         H32 = sum(abs(imp_res.h32).^2);
         H41 = sum(abs(imp_res.h41).^2);
         H42 = sum(abs(imp_res.h42).^2);
         
         H13 = sum(abs(imp_res.h13).^2);
         H14 = sum(abs(imp_res.h14).^2);
         H23 = sum(abs(imp_res.h23).^2);
         H24 = sum(abs(imp_res.h24).^2);
         
         normH = sqrt((H11 + H12 + H21 + H22 + ...
                      H33 + H34 + H43 + H44 + ...
                      H31 + H32 + H41 + H42 + ...
                      H13 + H14 + H23 + H24)/Nrx);
         
         imp_res.h11 = imp_res.h11 ./ normH;
         imp_res.h12 = imp_res.h12 ./ normH;
         imp_res.h21 = imp_res.h21 ./ normH;
         imp_res.h22 = imp_res.h22 ./ normH;
         
         imp_res.h13 = imp_res.h13 ./ normH;
         imp_res.h14 = imp_res.h14 ./ normH;
         imp_res.h23 = imp_res.h23 ./ normH;
         imp_res.h24 = imp_res.h24 ./ normH;
         
         imp_res.h31 = imp_res.h31 ./ normH;
         imp_res.h32 = imp_res.h32 ./ normH;
         imp_res.h41 = imp_res.h41 ./ normH;
         imp_res.h42 = imp_res.h42 ./ normH;
         
         imp_res.h33 = imp_res.h33 ./ normH;
         imp_res.h34 = imp_res.h34 ./ normH;
         imp_res.h43 = imp_res.h43 ./ normH;
         imp_res.h44 = imp_res.h44 ./ normH;
         
    case 5,% conf#5
         Nrx = 2; % number of Rx chains
         H11 = sum(abs(imp_res.h11).^2);
         H12 = sum(abs(imp_res.h12).^2);
         normH = sqrt((H11 + H12)/Nrx);
         imp_res.h11 = imp_res.h11 ./ normH;
         imp_res.h12 = imp_res.h12 ./ normH;
         
    otherwise,
       error('Prohibited value of "ant_type" parameter');
end
   
end
