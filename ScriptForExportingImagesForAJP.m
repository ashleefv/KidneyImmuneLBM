%insert this script filename after whichever figure you want to export from
%MATLAB to tiff for AJP Renal 
% Change the filename accordingly, then use Image Magick for the last step
% from the command prompt
% Use the pdfs in the manuscript and the tiff files for the journal
fig= gcf;

% Set the desired output size in inches (e.g., 3.4 inches wide for single column)
% For Figs 4-6, widthInches = 20; heightInches = 10; for 2x4 grid 
% For Figs 7-12, widthInches = 15; heightInches = 10; for 2x3 grid or 1x3
% as height scales by export graphics
% For Figs 13-15, widthInches = 10; heightInches = 10; for 2x2 grid 
%widthInches = 20;
heightInches = 10; % Adjust based on your figure's aspect ratio


% Get current size in inches
set(fig, 'Units', 'Inches');
figPos = get(fig, 'Position');


% Set figure size
set(fig, 'Position', [1, 1, widthInches, heightInches]);
set(fig, 'PaperUnits', 'Inches');
set(fig, 'PaperSize', [widthInches, heightInches]);
figPos = get(fig, 'Position');
set(fig, 'PaperPositionMode', 'manual');
fig = gcf;


% Apply to all subplots
ax = findall(gcf, 'Type', 'axes');
for k = 1:length(ax)
    box(ax(k), 'on'); % or 'off' depending on your preference
end

exportgraphics(fig,[filename, 'big.tiff'],'Resolution',1200)
exportgraphics(fig,[filename, '.pdf'],'Resolution',1200)

% do the following in the command prompt with image magick installed
%magick Fig4big.tiff -resize 36.90% -density 1200 Fig4.tiff

% Figs 4-6 big are 16.26x9.23, scale 6/16.26 = 36.90%
% Fig 7 big is 12.36x9.23, scale by same
% Fig 8 big is 12.44x9.23
% Figs 8-9 big are 12.44x4.49, scale by same
% Figs 10-12 big are 12.48x9.23
% Figs 13-14 big are 8.50x9.23 and 8.54x9.23
% Fig 15 big is 8.53x9.16