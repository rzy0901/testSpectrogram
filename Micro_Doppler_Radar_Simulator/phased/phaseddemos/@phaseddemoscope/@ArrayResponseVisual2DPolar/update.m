function update(this)
%UPDATE   

%   Copyright 2010-2011 The MathWorks, Inc.

% Cache the new frame data
source = this.Application.DataSource;

% get data, we only show the last frame

newData = source.getRawData;
if isempty(newData)
    return;
end

%% line display, add data to scope

newData = newData{1};
theta1 = -180:180;
theta = theta1.*pi./180;
if all(newData==0)
    [xx, yy] = pol2cart(theta',1);
else
    [xx, yy] = pol2cart (theta',newData);
end

if isempty(this.Lines)
    this.Lines = line('Parent',this.Axes,...
        'XData',xx ,'YData',yy,'Color','b');
else
    set(this.Lines,'XData',xx,'YData',yy,'Color','b');
end
end


% [EOF]
