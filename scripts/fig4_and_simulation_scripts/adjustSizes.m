function adjustSizes(ax, lineWidth, fntSize)
if nargin<2
    lineWidth = 1;
end
if nargin<3
    fntSize = 14;
end

set(ax,'FontSize', fntSize);
set(ax,'LineWidth', lineWidth);
set(get(ax,'Title'),'FontSize', fntSize);
set(get(ax,'XLabel'),'FontSize', fntSize);
set(get(ax,'YLabel'),'FontSize', fntSize);
end