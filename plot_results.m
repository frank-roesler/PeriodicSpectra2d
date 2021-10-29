function plot_results(L, Spectrum, Spectrum_temp, window, Determinant)
%     plots the current spectrum and search window
    cast(Spectrum, 'like', 1+1i);
    cast(Spectrum_temp, 'like', 1+1i);
    plot(Spectrum+0.000000001i,'.','MarkerEdgeColor',[0.3 0.3 1], 'MarkerSize',10);
    xlim([min(real(L(:))) max(real(L(:)))]);
    ylim([min(imag(L(:))) max(imag(L(:)))]);
    hold on
    scatter(real(window(:)), imag(window(:)),'filled','SizeData',5);
    plot(Spectrum_temp+0.000000001i,'.','MarkerEdgeColor',[0 0 0], 'MarkerSize',10);
    alpha(0.4)
%     contour(real(window), imag(window),log(abs(Determinant)),10);
    hold off
    drawnow;
    pause(0.001)
end
