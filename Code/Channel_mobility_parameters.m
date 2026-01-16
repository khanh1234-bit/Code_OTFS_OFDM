% maximum user equipment speed
max_UE_speed_kmh=100;
max_UE_speed = max_UE_speed_kmh*(1000/3600);
% maximum Doppler spread (one-sided)
nu_max = (max_UE_speed*fc)/(c);
% maximum normalized Doppler spread (one-sided)
k_max = nu_max/Doppler_resolution;