function h = plotPatch(mn, er, col)
    xs = 1:length(mn);
    h = plot(xs,mn,'-','Color', col, 'LineWidth',2);
    patch([xs, fliplr(xs)],[mn+er, fliplr(mn-er)],col,'FaceAlpha',0.5,'EdgeColor','none');
end
