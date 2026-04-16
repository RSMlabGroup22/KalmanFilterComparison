function [ds, T0] = loadAndSyncSensors(dataDir)    
    imu_file = fullfile(dataDir, 'clean_imu.csv');
    gnss_file = fullfile(dataDir, 'clean_gnss.csv');
    if ~exist(imu_file, 'file')
        error('Error: Couldnt find %s!', imu_file);
    end
    if ~exist(gnss_file, 'file')
        error('Error: Couldnt find %s!', gnss_file);
    end
    opts_imu = detectImportOptions(imu_file, 'VariableNamingRule', 'preserve');
    raw_imu = readtable(imu_file, opts_imu);
    opts_gnss = detectImportOptions(gnss_file, 'VariableNamingRule', 'preserve');
    raw_gnss = readtable(gnss_file, opts_gnss);
    timeCol_imu = findTimeCol(raw_imu);
    t_imu = raw_imu.(timeCol_imu);    
    t_imu = double(t_imu); 
    timeCol_gnss = findTimeCol(raw_gnss);
    t_gnss = double(raw_gnss.(timeCol_gnss));
    
    T0 = min(t_imu(1), t_gnss(1));
    
    ds.time = t_imu - T0;
    
    if length(ds.time) > 1
        dt_vector = builtin('diff', ds.time); 
        ds.dt = mean(dt_vector);
    else
        ds.dt = 0.01;
    end
    ds.datetime = datetime('now') + seconds(ds.time); 
    ds.ax = findColumnData(raw_imu, {'ax', 'acc_x', 'accX', 'accel_x'});
    ds.ay = findColumnData(raw_imu, {'ay', 'acc_y', 'accY', 'accel_y'});
    ds.az = findColumnData(raw_imu, {'az', 'acc_z', 'accZ', 'accel_z'});
    
    ds.gx = findColumnData(raw_imu, {'gx', 'gyr_x', 'gyrX', 'gyroX'});
    ds.gy = findColumnData(raw_imu, {'gy', 'gyr_y', 'gyrY', 'gyroY'});
    ds.gz = findColumnData(raw_imu, {'gz', 'gyr_z', 'gyrZ', 'gyroZ'});
    ds.gps.hacc = findColumnData(raw_gnss, {'hAcc'});
    ds.gps.heading = findColumnData(raw_gnss, {'heading'});
    ds.gps.lat = findColumnData(raw_gnss, {'lat'});
    ds.gps.lon = findColumnData(raw_gnss, {'lon'});
    ds.gyr_z = findColumnData(raw_imu, {'gyr_z'});


    ds.acc = [ds.ax, ds.ay, ds.az];
    
    ds.gps.time = t_gnss - T0;
    ds.gps.lat = findColumnData(raw_gnss, {'lat', 'latitude'});
    ds.gps.lon = findColumnData(raw_gnss, {'lon', 'lng', 'longitude'});
    ds.gps.h   = findColumnData(raw_gnss, {'alt', 'altitude', 'h', 'height'});

    ds.gps.hacc = findColumnData(raw_gnss, {'acc', 'accuracy', 'hacc'});
    if all(ds.gps.hacc == 0)
        ds.gps.hacc = ones(length(t_gnss), 1) * 10;
    end
    
    ds.gps.datetime = datetime('now') + seconds(ds.gps.time);
end


function colName = findTimeCol(tbl)
    cols = tbl.Properties.VariableNames;
    idx = contains(cols, 'time', 'IgnoreCase', true) | contains(cols, 't', 'IgnoreCase', true);
    if any(idx)
        colName = cols{find(idx, 1)}; 
    else
        error('Veri setinde "time" veya "t" isimli bir zaman sütunu bulunamadı!');
    end
end

function colData = findColumnData(tbl, keywords)
    cols = tbl.Properties.VariableNames;
    colData = zeros(height(tbl), 1); % Bulunamazsa 0 ile doldur
    
    for i = 1:length(keywords)
        idx = strcmpi(cols, keywords{i}) | contains(cols, keywords{i}, 'IgnoreCase', true);
        if any(idx)
            colName = cols{find(idx, 1)};
            colData = tbl.(colName);
            return;
        end
    end
end