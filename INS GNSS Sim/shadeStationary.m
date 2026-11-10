%% Local function (must sit at the end of the script)
function shadeStationary(ax, t, iStart, iEnd)
    yl = ylim(ax);
    for j = 1:numel(iStart)
        x0 = t(iStart(j));  x1 = t(iEnd(j));
        p = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                  'HandleVisibility', 'off');
        uistack(p, 'bottom');        % keep bands behind the data
    end
    ylim(ax, yl);                    % stop the patch changing the limits
end