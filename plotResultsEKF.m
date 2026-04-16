function plotResultsEKF(ds, ekf,outputDir)

    f = figure('Visible','on', 'Position', [100 100 1200 800]);
    

    geoplot(ds.gps.lat, ds.gps.lon, 'ko', 'DisplayName', 'GPS Raw', 'MarkerSize', 3); hold on;
    

    geoplot(ekf.lat, ekf.lon, 'g.-', 'LineWidth', 2.5, 'MarkerSize', 10, 'DisplayName', 'EKF');

    
    geobasemap('streets'); 
    title('Extended Kalman Filter');
    legend('Location', 'best');
    
    saveas(f, fullfile(outputDir, 'ekf_map.png'));
    
    f2 = figure('Visible','off');
    plot(ekf.time, ekf.n, 'g', 'LineWidth', 1.5); hold on;

    
    plot(ds.gps.time, ds.gps.lat, 'k:', 'LineWidth', 2);
    
    title('Kuzey (North) Pozisyonu Karşılaştırması');
    xlabel('Zaman (s)'); ylabel('Mesafe (m)');
    legend('Classic EKF', 'GPS Raw');
    grid on;
    saveas(f2, fullfile(outputDir, 'north_pos_comparison.png'));
    

   