function privrocplot(funcname,xval,yval,param,plotfun,xlabeltext,paramtext,paramunit,typeintitle)
%PRIVROCPLOT plot roc curve
%   PRIVROCPLOT(xval,yval,param,plotfun,xlabeltext,typeintitle) plots the
%   roc curve as yval vs. xval with changing parameter param.  The plot
%   uses the function handle plotfun as the plotting method. xlabeltext,
%   paramtext, paramunit and typeintitle are used to form x axis label and
%   figure title.

%   Copyright 2008-2010 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, pg 329

%#codegen

if ~isempty(coder.target)
    coder.internal.assert(false,'phased:rocsnr:invalidCodegenOutput',funcname);
end
xval = xval(:);
paramlen = length(param);
legendstr{paramlen}=0;  % preallocate
h = plotfun(xval,yval);
ax = get(h,'parent');
if iscell(ax)
    ax = ax{1};
end
ylim = get(ax,'YLim');
textyinc = (ylim(2)-ylim(1))*0.5/paramlen;
textystart = (ylim(2)+ylim(1))*0.5;
for k = 1:paramlen,
    legendstr{k} = [paramtext,'=',num2str(param(k)),paramunit];
    texty = textystart + textyinc*(k-1);
    [~,textxidx] = min(abs(yval(:,k)-texty));
    textx = xval(textxidx);
    text(textx,texty,legendstr{k});
end
grid on;

typeintitle = regexprep(typeintitle,'Nonfluctuating','Nonfluctuating ');
title({typeintitle, ' Receiver Operating Characteristic (ROC) Curves'});
xlabel(xlabeltext);
ylabel('P_d');


% [EOF]
