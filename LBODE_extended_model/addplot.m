function addplot(h,whichCurve,test_i,whichsubplot,ax1)
    Curve = h(whichCurve);
    x = get(Curve,'Xdata');
    y = get(Curve,'Ydata');
    % Copy all relevant properties
    props = {'Color', 'LineStyle', 'LineWidth', 'Marker', 'MarkerSize', ...
             'MarkerEdgeColor', 'MarkerFaceColor', 'DisplayName'};
    f2=figure(6);
    subplot(whichsubplot)
    hold on
    newLine=plot(x,y); hold on
    for p = 1:length(props)
        try
            set(newLine, props{p}, get(Curve, props{p}));
        catch
            % Some properties might not be set, so we skip errors
        end
    end
end
