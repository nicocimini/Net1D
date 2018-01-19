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

i = i + 1;
Variables(i).Name = 'ta_sys';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'K';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'systematic error of air_temperature';
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

i = i + 1;
Variables(i).Name = 'hua_sys';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-3';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'systematic error of absolute humidity';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This variable specifies the systematic uncertainty of hua as a function of height above ground';

i = i + 1;
Variables(i).Name = 'prw';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-2';
Variables(i).Attributes(2).Name = 'standard_name';
Variables(i).Attributes(2).Value = 'atmosphere_mass_content_of_water_vapor';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'prw is the precipitable water, i.e. the vertically integrated amount of water vapor from the surface to TOA, also called Total water vapor content (TWVC)';

i = i + 1;
Variables(i).Name = 'prw_err';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-2';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'standard error of pwr';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'pwr random uncertainty';

i = i + 1;
Variables(i).Name = 'prw_sys';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-2';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'standard error of pwr';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'pwr systematic uncertainty';

i = i + 1;
Variables(i).Name = 'flag';
Variables(i).Attributes(1).Name = 'long_name';
Variables(i).Attributes(1).Value = 'quality control flags';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'quality control flags';
Variables(i).Attributes(3).Name = 'flag_masks';
Variables(i).Attributes(3).Value = '[1     2     4     8    16    32    64   128   256   512  1024]';
Variables(i).Attributes(3).Name = 'flag_meanings';
Variables(i).Attributes(3).Value = 'manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3 iwv_lwp_threshold';
Variables(i).Attributes(4).Name = 'FillValue';
Variables(i).Attributes(4).Value = '0';
Variables(i).Attributes(5).Name = 'comment';
Variables(i).Attributes(5).Value = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. A Fillvalue of 0 means that data has not been flagged. Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; tb valid range: [  2.70, 330.00] in K; prw valid range: [   0.,  100.] in kgm-2; clwvi (zeroing not applied) valid range: [-0.2,  3.0] in kgm-2;';

i = i + 1;
Variables(i).Name = 'prw_offset';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-2';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'pwr offset based on Tb offset';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This field is not used in GAIA-CLIM';

i = i + 1;
Variables(i).Name = 'prw_off_zenith';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-2';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'off zenith path integrated water vapor';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This field is not used in GAIA-CLIM';

i = i + 1;
Variables(i).Name = 'prw_off_zenith_offset';
Variables(i).Attributes(1).Name = 'units';
Variables(i).Attributes(1).Value = 'kg m-2';
Variables(i).Attributes(2).Name = 'long_name';
Variables(i).Attributes(2).Value = 'prw_off_zenith offset based on Tb offset';
Variables(i).Attributes(3).Name = 'FillValue';
Variables(i).Attributes(3).Value = '-999';
Variables(i).Attributes(4).Name = 'comment';
Variables(i).Attributes(4).Value = 'This field is not used in GAIA-CLIM';

return