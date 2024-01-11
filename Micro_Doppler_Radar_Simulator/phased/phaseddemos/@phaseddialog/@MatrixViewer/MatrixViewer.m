function this = MatrixViewer(block,varargin)
%MatrixViewer constructs phaseddialog.MatrixViewer object.
%   MatrixViewer(block) creates a phaseddialog.MatrixViewer dialog object associated with
%   a MatrixViewer block.  block may be an MatrixViewer block path (gcb), 
%   block handle (gcbh) or associated UDD object.

% Copyright 2014 The MathWorks, Inc.

this = phaseddialog.MatrixViewer(block);
this.init(block);
this.loadFromBlock; 

%[EOF]

