%insert this script figname after whichever figure you want to export from
%MATLAB to tiff for AJP Renal 
% Change the figname accordingly, then use Image Magick for the last step
% from the command prompt
% Use the pdfs in the manuscript and the tiff files for the journal
fig= gcf;

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

exportgraphics(fig,[figname, '.tiff'],'Resolution',600)
exportgraphics(fig,[figname, '.pdf'],'Resolution',600)
