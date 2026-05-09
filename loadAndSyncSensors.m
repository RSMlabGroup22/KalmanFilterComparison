function [ds, T0] = loadAndSyncSensors(dataDir)
% Yeni format kontrolü
isNewFormat = exist(fullfile(dataDir, 'Accelerometer.csv'), 'file');

if isNewFormat
    % --- YENİ FORMAT OKUMA ---
    acc_tbl = readtable(fullfile(dataDir, 'Accelerometer.csv'));
    gyr_tbl = readtable(fullfile(dataDir, 'Gyroscope.csv'));
    loc_tbl = readtable(fullfile(dataDir, 'Location.csv'));

    % Zaman (seconds_elapsed en güveniliridir)
    ds.time = acc_tbl.seconds_elapsed;
    ds.dt = mean(diff(ds.time));

    % İvme (x,y,z sütunlarını eşle)
    % Not: Yeni veriler linear acceleration (yerçekimsiz) ise
    % mevcut EKF kodunla uyum için Z eksenine 9.81 ekliyoruz.
    ds.ax = acc_tbl.x;
    ds.ay = acc_tbl.y;
    ds.az = acc_tbl.z + 9.81; % Uyumluluk için yerçekimi eklendi

    % Jiroskop
    ds.gx = gyr_tbl.x;
    ds.gy = gyr_tbl.y;
    ds.gz = gyr_tbl.z;
    ds.gyr_z = gyr_tbl.z; % EKF'nin beklediği ek alan

    % GNSS Verileri
    ds.gps.time = loc_tbl.seconds_elapsed;
    ds.gps.lat = loc_tbl.latitude;
    ds.gps.lon = loc_tbl.longitude;
    ds.gps.h   = loc_tbl.altitude;
    ds.gps.hacc = loc_tbl.horizontalAccuracy;

    % Heading (bearing sütunu -1 ise 0 kabul et veya Orientation'dan çek)
    ds.gps.heading = loc_tbl.bearing;
    ds.gps.heading(ds.gps.heading == -1) = 0;

    T0 = 0; % seconds_elapsed zaten 0'dan başlar
else
    % --- ESKİ FORMAT OKUMA (Mevcut mantığın korunmuş hali) ---
    imu_file = fullfile(dataDir, 'clean_imu.csv');
    gnss_file = fullfile(dataDir, 'clean_gnss.csv');

    opts_imu = detectImportOptions(imu_file, 'VariableNamingRule', 'preserve');
    raw_imu = readtable(imu_file, opts_imu);
    opts_gnss = detectImportOptions(gnss_file, 'VariableNamingRule', 'preserve');
    raw_gnss = readtable(gnss_file, opts_gnss);

    timeCol_imu = findTimeCol(raw_imu);
    t_imu = double(raw_imu.(timeCol_imu));
    timeCol_gnss = findTimeCol(raw_gnss);
    t_gnss = double(raw_gnss.(timeCol_gnss));

    T0 = min(t_imu(1), t_gnss(1));
    ds.time = t_imu - T0;
    ds.dt = mean(diff(ds.time));

    ds.ax = findColumnData(raw_imu, {'ax', 'acc_x', 'x'});
    ds.ay = findColumnData(raw_imu, {'ay', 'acc_y', 'y'});
    ds.az = findColumnData(raw_imu, {'az', 'acc_z', 'z'});
    ds.gx = findColumnData(raw_imu, {'gx', 'gyr_x', 'x'});
    ds.gy = findColumnData(raw_imu, {'gy', 'gyr_y', 'y'});
    ds.gz = findColumnData(raw_imu, {'gz', 'gyr_z', 'z'});
    ds.gyr_z = ds.gz;

    ds.gps.time = t_gnss - T0;
    ds.gps.lat = findColumnData(raw_gnss, {'lat', 'latitude'});
    ds.gps.lon = findColumnData(raw_gnss, {'lon', 'longitude'});
    ds.gps.h   = findColumnData(raw_gnss, {'h', 'altitude'});
    ds.gps.hacc = findColumnData(raw_gnss, {'hacc', 'horizontalAccuracy'});
    ds.gps.heading = findColumnData(raw_gnss, {'heading', 'bearing'});
end

ds.acc = [ds.ax, ds.ay, ds.az];
ds.datetime = datetime('now') + seconds(ds.time);
ds.gps.datetime = datetime('now') + seconds(ds.gps.time);
end
