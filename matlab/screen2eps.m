function screen2eps(filename,format)
%SCREEN2JPEG Generate a JPEG file of the current figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2JPEG('filename') saves the current figure to the
% JPEG file "filename".
%
% Sean P. McCarthy
% Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved

if nargin < 1
error('Not enough input arguments!')
end

set(gcf,'PaperPositionMode','auto');
set(gca, 'xcolor', [0 0 0],...
         'ycolor', [0 0 0],...
         'color', 'none');
% set(gcf, 'color', 'none',...
%          'inverthardcopy', 'off');
     
if(strcmp(format,'eps'))
    print('-depsc2', filename, '-r100');
elseif(strcmp(format,'png'))
    print('-dpng', filename, '-r100');
end

end