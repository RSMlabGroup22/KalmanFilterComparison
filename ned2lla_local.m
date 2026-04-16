function [lat, lon] = ned2lla_local(n, e, lat0, lon0, h0)
    Re = 6378137; 
    f = 1/298.257223563;
    deg2rad_const = pi/180;
    
    rn = Re * (1 - f*(2-f)) / (1 - (f*(2-f)*sin(lat0*deg2rad_const)^2))^1.5;
    rm = Re / sqrt(1 - (f*(2-f)*sin(lat0*deg2rad_const)^2));
    
    dLat = n / (rn + h0);
    dLon = e / ((rm + h0) * cos(lat0*deg2rad_const));
    
    lat = lat0 + dLat / deg2rad_const;
    lon = lon0 + dLon / deg2rad_const;
end
