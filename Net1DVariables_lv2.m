% Network 1DVAR+RTTOV retrieval: prepares variables' attributes as in the HDPC2 level2 format
%
%
% Net1DVariables_lv2 prepares attributes for variables in the HDPC2 level2 format
%
% Variables = Net1DVariables_lv2
%
% 

function Variables = Net1DVariables_lv2


i = 0;

i = i + 1;
Variables(i).Name = 'time';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'seconds since 1970-01-01 00:00:00 UTC';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'time';
Variables(i).Attributes(3).Name = 'bounds';
Variables(i).Attributes(3).Value = 'time_bnds';

i = i + 1;
Variables(i).Name = 'time_bnds';

i = i + 1;
Variables(i).Name = 'lat';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'degree_north';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'latitude';
          
i = i + 1;
Variables(i).Name = 'lon';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'degree_east';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'longitude';
           
i = i + 1;
Variables(i).Name = 'zsl';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'm';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'altitude';
Variables(i).Attributes(3).Name = 'long_name';
Variables(i).Attributes(3).Value = 'altitude above mean sea level';

i = i + 1;
Variables(i).Name = 'azi';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'degree';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'sensor_azimuth_angle';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = '0=North, 90=East, 180=South, 270=West';

i = i + 1;
Variables(i).Name = 'ele';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'degree';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'retrieval elevation angle';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This variable specifies the elevation angle at which retrievals have been derived.';

i = i + 1;
Variables(i).Name = 'height';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'm';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'height';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';

i = i + 1;
Variables(i).Name = 'ta';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'K';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'air_temperature';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'ta profiles are vertical profiles over the measurement site';

i = i + 1;
Variables(i).Name = 'ta_offset';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'K';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'air_temperature offset correction based on brightness temperature offset';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '0';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'In order to obtain the un-corrected ta profile, add this offset to ta. This variable is intended for expert use only.';

i = i + 1;
Variables(i).Name = 'ta_err';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'K';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'standard error of air_temperature';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';

% for hua

i = i + 1;
Variables(i).Name = 'ele_ret';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'degree';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'retrieval elevation angle';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This variable specifies the elevation angle at which retrievals have been derived.';

i = i + 1;
Variables(i).Name = 'hua';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-3';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'absolute humidity';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'hua profiles are vertical profiles over the measurement site';

i = i + 1;
Variables(i).Name = 'hua_offset';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-3';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'absolute humidity offset correction based on brightness temperature offset';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '0';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'In order to obtain the un-corrected hua profile, add this offset to hua. This variable is intended for expert use only.';

i = i + 1;
Variables(i).Name = 'hua_err';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-3';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'standard error of absolute humidity';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This variable specifies the uncertainty of hua as a function of height above ground';

return