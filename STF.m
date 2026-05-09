function imu_out = imuLoc(ds)
    % Veri uzunluğu ve ön bellek tahsisi
    N = length(ds.time);
    imu_out.n = zeros(N, 1);
    imu_out.e = zeros(N, 1);
    imu_out.vn = zeros(N, 1);
    imu_out.ve = zeros(N, 1);
    imu_out.lat = zeros(N, 1);
    imu_out.lon = zeros(N, 1);
    imu_out.time = ds.time;
    
    % Referans noktası (İlk konum verisinden çekiyoruz)
    ref_lat = ds.lat(1);
    ref_lon = ds.lon(1);
    
    % Durum Vektörü: x = [Kuzey; Doğu; Hız_N; Hız_E]
    x = [0; 0; 0; 0]; 
    
    dt = ds.dt;
    
    % Başlangıç Yönelimi ve Hızı (İlk iki noktadan kestirim)
    [n2, e2] = lla2ned_local_internal(ds.lat(2), ds.lon(2), ref_lat, ref_lon);
    psi = atan2(e2, n2); 
    
    x(3) = n2 / dt;
    x(4) = e2 / dt;
    
    % Birim dönüşümleri (STF kodundaki mantık)
    scale_acc = 1.0; scale_gyro = 1.0;
    if mean(abs(ds.az)) < 2.0, scale_acc = 9.80665; end
    if mean(abs(ds.gz)) > 10, scale_gyro = pi/180; end
    
    % --- DÖNGÜ ---
    for k = 1:N
        % 1. Yönelim (Yaw) Tahmini
        gz = ds.gz(k) * scale_gyro;
        if k > 1 % Başlangıç noktasındaki açıyı bozmamak için
            psi = psi + gz * dt;
        end
        
        % 2. İvme Okumaları ve NED Eksenine Dönüşüm
        ax = ds.ax(k) * scale_acc; 
        
        ay = 0; % Eğer yanal ivme verisi dataset'te varsa eklenebilir
        if isfield(ds, 'ay')
            ay = ds.ay(k) * scale_acc;
        end
        
        aN = ax * cos(psi) - ay * sin(psi);
        aE = ax * sin(psi) + ay * cos(psi);
        
        % 3. Durum Güncellemesi (Euler İntegrasyonu)
        if k > 1
            % Konum: x = x_prev + v_prev * dt + 0.5 * a * dt^2
            x(1) = x(1) + x(3) * dt + 0.5 * aN * dt^2;
            x(2) = x(2) + x(4) * dt + 0.5 * aE * dt^2;
            
            % Hız: v = v_prev + a * dt
            x(3) = x(3) + aN * dt;
            x(4) = x(4) + aE * dt;
        end
        
        % Çıktıları kaydet
        imu_out.n(k) = x(1); 
        imu_out.e(k) = x(2);
        imu_out.vn(k) = x(3); 
        imu_out.ve(k) = x(4);
        [imu_out.lat(k), imu_out.lon(k)] = ned2lla_local_internal(x(1), x(2), ref_lat, ref_lon);
    end
end

% --- YARDIMCI FONKSİYONLAR (Dosya içine gömülü) ---
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