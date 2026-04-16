function plotResultsSTF(ds, stf, outputDir)
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    % 1. Harita Üzerinde İzleme
    f1 = figure('Name', 'STF Map Result', 'Position', [100 100 1200 800]);
    geoplot(ds.gps.lat, ds.gps.lon, 'k.', 'DisplayName', 'Ham GPS (RTK)', 'MarkerSize', 8); hold on;
    geoplot(stf.lat, stf.lon, 'r-', 'LineWidth', 2, 'DisplayName', 'Strong Tracking Filter');
    geobasemap('streets'); 
    title('Strong Tracking Filter (STF) - Tren Yörünge Analizi');
    legend('Location', 'best');
    saveas(f1, fullfile(outputDir, 'stf_map.png'));

    % 2. Eksen Bazlı Karşılaştırma
    f2 = figure('Name', 'STF Axis Comparison', 'Position', [150 150 1000 600]);
    
    % Metre dönüşümü için referanslar
    Re = 6378137; deg2rad = pi/180;
    lat0 = ds.gps.lat(1); lon0 = ds.gps.lon(1);
    gps_n = (ds.gps.lat - lat0) * deg2rad * Re;
    gps_e = (ds.gps.lon - lon0) * deg2rad * Re .* cos(lat0 * deg2rad);

    subplot(2,1,1);
    plot(stf.time, stf.n, 'r', 'LineWidth', 1.5); hold on;
    plot(ds.gps.time, gps_n, 'k:', 'LineWidth', 1);
    title('Kuzey (North) Pozisyonu - STF'); 
    grid on; legend('STF', 'GPS Raw');

    subplot(2,1,2);
    plot(stf.time, stf.e, 'r', 'LineWidth', 1.5); hold on;
    plot(ds.gps.time, gps_e, 'k:', 'LineWidth', 1);
    title('Doğu (East) Pozisyonu - STF'); 
    grid on; legend('STF', 'GPS Raw');
    
    saveas(f2, fullfile(outputDir, 'stf_pos_comparison.png'));
end
