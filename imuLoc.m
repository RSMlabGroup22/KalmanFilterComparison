function imu_out = imuLoc(ds)
    % Veri uzunluğu ve ön bellek tahsisi
    N = length(ds.time);
    imu_out.n = zeros(N, 1);
    imu_out.e = zeros(N, 1);
    imu_out.vn = zeros(N, 1);
    imu_out.ve = zeros(N, 1);
    imu_out.psi = zeros(N, 1);
    imu_out.lat = zeros(N, 1);
    imu_out.lon = zeros(N, 1);
    
    dt = ds.dt;
    v_max = 38.5; % Maksimum hız sınırı (m/s)

    % --- MANUEL 30Hz LOW-PASS FILTER (Bütünleşik) ---
    fs = 1 / dt;
    fc = 5;
    gamma = tan(pi * fc / fs);
    D = gamma^2 + sqrt(2) * gamma + 1;
    b = [gamma^2/D, 2*gamma^2/D, gamma^2/D];
    a = [1, 2*(gamma^2-1)/D, (gamma^2 - sqrt(2)*gamma + 1)/D];

    % Filtrelenmiş veri dizileri
    ax_f = zeros(N, 1);
    ay_f = zeros(N, 1);
    gz_f = zeros(N, 1);

    % İvme Y (ay) kontrolü
    if isfield(ds, 'ay'), ay_raw = ds.ay; else, ay_raw = zeros(N, 1); end

    % Manuel Filtreleme Döngüsü (Hata payını sıfıra indirmek için inline)
    ax_f(1:2) = ds.ax(1:2);
    ay_f(1:2) = ay_raw(1:2);
    gz_f(1:2) = ds.gz(1:2);
    
    for i = 3:N
        ax_f(i) = b(1)*ds.ax(i) + b(2)*ds.ax(i-1) + b(3)*ds.ax(i-2) - a(2)*ax_f(i-1) - a(3)*ax_f(i-2);
        ay_f(i) = b(1)*ay_raw(i) + b(2)*ay_raw(i-1) + b(3)*ay_raw(i-2) - a(2)*ay_f(i-1) - a(3)*ay_f(i-2);
        gz_f(i) = b(1)*ds.gz(i) + b(2)*ds.gz(i-1) + b(3)*ds.gz(i-2) - a(2)*gz_f(i-1) - a(3)*gz_f(i-2);
    end
    % ------------------------------------------------

    % 1. Başlangıç Konumu
    lat0 = ds.gps.lat(1);
    lon0 = ds.gps.lon(1);
    imu_out.lat(1) = lat0;
    imu_out.lon(1) = lon0;
    imu_out.n(1) = 0;
    imu_out.e(1) = 0;
    
    % 2. Başlangıç Yönelimi ve Hızı
    [n2, e2] = lla2ned_local_internal(ds.gps.lat(5), ds.gps.lon(5), lat0, lon0);
    imu_out.psi(1) = atan2(e2, n2);
    
    v_n_init = n2 / (dt * 4);
    v_e_init = e2 / (dt * 4);
    
    v_total_init = sqrt(v_n_init^2 + v_e_init^2);
    if v_total_init > v_max
        ratio = v_max / v_total_init;
        v_n_init = v_n_init * ratio;
        v_e_init = v_e_init * ratio;
    end
    imu_out.vn(1) = v_n_init;
    imu_out.ve(1) = v_e_init;
    
    % 3. IMU Dead Reckoning Döngüsü
    for k = 2:N
        % a. Yönelim Güncellemesi (Filtrelenmiş)
        imu_out.psi(k) = imu_out.psi(k-1) + gz_f(k) * dt;
        
        % b. İvme Dönüşümü (Filtrelenmiş)
        a_n = ax_f(k) * cos(imu_out.psi(k)) - ay_f(k) * sin(imu_out.psi(k));
        a_e = ax_f(k) * sin(imu_out.psi(k)) + ay_f(k) * cos(imu_out.psi(k));
        
        % c. Hız Güncellemesi
        v_n_new = imu_out.vn(k-1) + a_n * dt;
        v_e_new = imu_out.ve(k-1) + a_e * dt;
        
        % Hız Sınırlama
        v_mag = sqrt(v_n_new^2 + v_e_new^2);
        if v_mag > v_max
            scale = v_max / v_mag;
            v_n_new = v_n_new * scale;
            v_e_new = v_e_new * scale;
        end
        
        imu_out.vn(k) = v_n_new;
        imu_out.ve(k) = v_e_new;
        
        % d. Konum Güncellemesi
        imu_out.n(k) = imu_out.n(k-1) + imu_out.vn(k) * dt;
        imu_out.e(k) = imu_out.e(k-1) + imu_out.ve(k) * dt;
        
        % e. LLA Dönüşümü
        [imu_out.lat(k), imu_out.lon(k)] = ned2lla_local_internal(imu_out.n(k), imu_out.e(k), lat0, lon0);
    end
end

% --- YARDIMCI FONKSİYONLAR (Burayı silmeyin) ---
function [n, e] = lla2ned_local_internal(lat, lon, lat0, lon0)
    Re = 6378137; deg2rad = pi/180;
    n = (lat - lat0) * deg2rad * Re;
    e = (lon - lon0) * deg2rad * Re * cos(lat0 * deg2rad);
end

function [lat, lon] = ned2lla_local_internal(n, e, lat0, lon0)
    Re = 6378137; rad2deg = 180/pi;
    lat = lat0 + (n / Re) * rad2deg;
    lon = lon0 + (e / (Re * cos(lat0 * pi/180))) * rad2deg;
end