function  flag = plot_semicircular_arch(coords, disp, figname)
    hf = figure('visible','off');
    plot(coords(:,1), coords(:,2), 'k-', "linewidth", 2);
    hold on
    plot(coords(:,1)+disp(1:3:end), coords(:,2)+disp(2:3:end), 'b-', "linewidth", 2)
    axis([-150 150 -150 150])
%    labels=legend("Original","Deformed");
    legend("Original","Deformed", "location", "northwest");
%    legend("fontsize", 18);
    print(hf, figname, "-dpdf");
    flag = true;
endfunction