function export_graph(figHandle, fileName, folderName, pngEn)

% Create folder if it doesn't exist
if ~exist(folderName,'dir')
    mkdir(folderName);
end

% Change background color
oldColor                 = figHandle.Color;
figHandle.Color          = 'none';
figHandle.InvertHardcopy = 'off';

% Paper size
set(figHandle,'PaperPositionMode','auto');

% Export to EPS
print(figHandle,'-painters','-depsc2',[folderName,'/',fileName,'.eps']);

% Export to PNG
if pngEn == true
    figHandle.InvertHardcopy = 'on';
    print(figHandle,'-dpng','-r150',[folderName,'/',fileName,'.png']);
end

% Change back background color
figHandle.Color = oldColor;
figHandle.InvertHardcopy = 'on';

end
