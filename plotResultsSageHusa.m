function plotResultsSageHusa(ds, sh, outputDir)
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    % 1. Harita Üzerinde İzleme
    f1 = figure('Name', 'Sage-Husa Map Result', 'Position', [100 100 1200 800]);
    geoplot(ds.gps.lat, ds.gps.lon, 'k.', 'DisplayName', 'Ham GPS', 'MarkerSize', 8); hold on;
    geoplot(sh.lat, sh.lon, 'b-', 'LineWidth', 2, 'DisplayName', 'Sage-Husa AKF');
    geobasemap('streets'); 
    title('Sage-Husa Adaptif Kalman Filtresi - Tren Yörünge Analizi');
    legend('Location', 'best');
    saveas(f1, fullfile(outputDir, 'sage_husa_map.png'));

    % 2. Pozisyon Karşılaştırma
    f2 = figure('Name', 'Sage-Husa Position Comparison', 'Position', [150 150 1000 600]);
    
    subplot(2,1,1);
    plot(sh.time, sh.n, 'b', 'LineWidth', 1.5); hold on;
    title('Kuzey (North) Pozisyonu - Sage-Husa'); 
    grid on; legend('SH-AKF');

    subplot(2,1,2);
    plot(sh.time, sh.e, 'b', 'LineWidth', 1.5); hold on;
    title('Doğu (East) Pozisyonu - Sage-Husa'); 
    grid on; legend('SH-AKF');
    
    saveas(f2, fullfile(outputDir, 'sage_husa_pos.png'));
end
