function d=lldistkm(latlon1,latlon2)
% format: [d1km d2km]=lldistkm(latlon1,latlon2)
% Distance:
% d1km: distance in km based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)

%--------------------------------------------------------------------------

    % Distance in km based on Haversine formula
    R = 6371; % Earth radius in km
    lat1 = deg2rad(latlon1(1));
    lat2 = deg2rad(latlon2(1));
    lon1 = deg2rad(latlon1(2));
    lon2 = deg2rad(latlon2(2));
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2;
    c = 2*atan2(sqrt(a), sqrt(1-a));
    d = R * c;
end
