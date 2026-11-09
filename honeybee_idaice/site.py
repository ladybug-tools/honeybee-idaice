"""A module for functions related to IDM site and climate."""
from ladybug.location import Location
from ladybug.designday import ASHRAETau


def project_info_to_climate_header(project_info=None):
    """Get a string for the header with climate and project information.

    Args:
        project_info: An optional ladybug ProjectInfo object that describes criteria
            about the project site. When specified, this will result in the location
            and north_angle written into the IDM. If a weather file has been associated
            with the project_info, a .PRN file will be written from the EPW and
            will be referenced within the IDM. (Default: None).
    """
    # parse the project_info to get the location and north_angle
    location, climate_location, year = None, None, None
    summer_day, winter_day = None, None
    if project_info is not None:
        location = project_info.location
        if project_info.location:
            location = project_info.location.duplicate()
            year = location.source if location.source is not None else None
            for attr in ('city', 'country', 'station_id', 'source'):
                if not getattr(location, attr):
                    setattr(location, attr, 'Unknown')
        climate_location = project_info.default_epw_file.location \
            if project_info.default_epw_file is not None else None
        des_days = project_info.default_design_days
        if len(des_days) != 0:
            summer_days = [dd for dd in des_days if dd.day_type == 'SummerDesignDay']
            if len(summer_days) != 0:
                summer_day = max(summer_days, key=lambda dd: dd.dry_bulb_condition.dry_bulb_max)
            winter_days = [dd for dd in des_days if dd.day_type == 'WinterDesignDay']
            if len(winter_days) != 0:
                winter_day = min(winter_days, key=lambda dd: dd.dry_bulb_condition.dry_bulb_max)

    # use the default location of Sweden if none has been specified
    default_location = Location(
        city='Kalmar', country='Sweden',
        latitude=56.73, longitude=16.3, time_zone=1, elevation=16,
        source='ASHRAE Fundamentals 2013',
        station_id='SMHI weather station in Kalmar'
    )
    if location is None:
        location = default_location
    if climate_location is None:
        climate_location = default_location.duplicate()
        climate_location.latitude = 56.68
        climate_location.elevation = 0
    if year is None:
        year = '1968'

    # use the default summer design day of Sweden if none has been specified
    if summer_day is None:
        sum_db_min, sum_db_max, sum_wb_max = 12.7, 27.1, 18.7
        sum_wind_dir, sum_wind_vel = 180, 4.5
        sum_tau_b, sum_tau_d, sum_clear = 0.343, 2.423, 1
    else:
        sum_db, sum_humid = summer_day.dry_bulb_condition, summer_day.humidity_condition
        sum_wind, sum_sky = summer_day.wind_condition, summer_day.sky_condition
        sum_db_max = sum_db.dry_bulb_max
        sum_db_min = sum_db.dry_bulb_max - sum_db.dry_bulb_range
        if sum_humid.humidity_type == 'Wetbulb':
            sum_wb_max = sum_humid.humidity_value
        else:
            sum_wb_max = max(sum_humid.hourly_dew_point_values(sum_db))
        sum_wind_dir, sum_wind_vel = sum_wind.wind_direction, sum_wind.wind_speed
        if isinstance(sum_sky, ASHRAETau):
            sum_tau_b, sum_tau_d, sum_clear = sum_sky.tau_b, sum_sky.tau_d, 1
        else:
            sum_tau_b, sum_tau_d, sum_clear = 0.343, 2.423, sum_sky.clearness

    # use the default winter design day of Sweden if none has been specified
    if winter_day is None:
        win_db_min, win_db_max, win_wb_max = -14.3, -9.3, -10.2
        win_wind_dir, win_wind_vel = 340, 1.7
        win_tau_b, win_tau_d, win_clear = 0.394, 2.366, 1
    else:
        win_db, win_humid = winter_day.dry_bulb_condition, winter_day.humidity_condition
        win_wind, win_sky = winter_day.wind_condition, winter_day.sky_condition
        win_db_max = win_db.dry_bulb_max
        win_db_min = win_db.dry_bulb_max - win_db.dry_bulb_range
        if win_humid.humidity_type == 'Wetbulb':
            win_wb_max = win_humid.humidity_value
        else:
            win_wb_max = max(win_humid.hourly_dew_point_values(win_db))
        win_wind_dir, win_wind_vel = win_wind.wind_direction, win_wind.wind_speed
        if isinstance(win_sky, ASHRAETau):
            win_tau_b, win_tau_d, win_clear = win_sky.tau_b, win_sky.tau_d, 1
        else:
            win_tau_b, win_tau_d, win_clear = 0.394, 2.366, win_sky.clearness

    # assemble the header for the climate info
    climate_header = \
        '((LOCATION :N "{city}" :T LOCATION :D "Data from {source}")\n' \
        ' (:PAR :N COUNTRY :V "{country}")\n' \
        ' (:PAR :N CITY :V "{city}")\n' \
        ' (:PAR :N LATITUDE :V {lat})\n' \
        ' (:PAR :N LONGITUDE :V {long})\n' \
        ' (:PAR :N ELEVATION :V {elev})\n' \
        ' (:PAR :N TIME_ZONE :V {tz})\n' \
        ' ((DESIGN-DAY :N WINTER-DESIGN-DAY)\n' \
        '  (:PAR :N DRY-BULB-MIN :F 128 :V {win_db_min})\n' \
        '  (:PAR :N DRY-BULB-MAX :F 128 :V {win_db_max})\n' \
        '  (:PAR :N WET-BULB-MAX :F 128 :V {win_wb_max})\n' \
        '  (:PAR :N WIND-DIR :F 128 :V {win_wind_dir})\n' \
        '  (:PAR :N WIND-VEL :F 128 :V {win_wind_vel})\n' \
        '  (:PAR :N TAU_B :F 128 :V {win_tau_b})\n' \
        '  (:PAR :N TAU_D :F 128 :V {win_tau_d})\n' \
        '  (:PAR :N RED-FAC :F 128 :V {win_clear}))\n' \
        ' ((DESIGN-DAY :N SUMMER-DESIGN-DAY)\n' \
        '  (:PAR :N DRY-BULB-MIN :F 128 :V {sum_db_min})\n' \
        '  (:PAR :N DRY-BULB-MAX :F 128 :V {sum_db_max})\n' \
        '  (:PAR :N WET-BULB-MAX :F 128 :V {sum_wb_max})\n' \
        '  (:PAR :N WIND-DIR :F 128 :V {sum_wind_dir})\n' \
        '  (:PAR :N WIND-VEL :F 128 :V {sum_wind_vel})\n' \
        '  (:PAR :N TAU_B :F 128 :V {sum_tau_b})\n' \
        '  (:PAR :N TAU_D :F 128 :V {sum_tau_d})\n' \
        '  (:PAR :N RED-FAC :F 128 :V {sum_clear}))\n' \
        ' (:ORES :N CLIMATE-DEF :V "{c_city}-{year}")\n' \
        ' (:PAR :N DESIGN-WEATHER-FILE :V "ASHRAE 2013\\026720.tbl"))\n' \
        '((CLIMATE-DEF :N "{c_city}-{year}" :T CLIMATE-DEF :D "{st_id} weather station in {c_city}")\n' \
        ' (:PAR :N FILENAME :V "{c_city}.prn")\n' \
        ' (:PAR :N STATION :V "{c_city}")\n' \
        ' (:PAR :N COUNTRY :V "{c_country}")\n' \
        ' (:PAR :N LATITUDE :V {c_lat})\n' \
        ' (:PAR :N LONGITUDE :V {c_long})\n' \
        ' (:PAR :N ELEVATION :V {c_elev})\n' \
        ' (:PAR :N TIME_ZONE :V {c_tz})\n' \
        ' (:PAR :N WIND-HEIGHT :V 10))\n'.format(
            city=location.city, country=location.country, year=year,
            lat=location.latitude, long=-location.longitude,
            elev=location.elevation, tz=-location.time_zone,
            st_id=climate_location.station_id, source=climate_location.source,
            c_city=climate_location.city, c_country=climate_location.country,
            c_lat=climate_location.latitude, c_long=-climate_location.longitude,
            c_elev=climate_location.elevation, c_tz=-climate_location.time_zone,
            win_db_min=win_db_min, win_db_max=win_db_max, win_wb_max=win_wb_max,
            win_wind_dir=win_wind_dir, win_wind_vel=win_wind_vel,
            win_tau_b=win_tau_b, win_tau_d=win_tau_d, win_clear=win_clear,
            sum_db_min=sum_db_min, sum_db_max=sum_db_max, sum_wb_max=sum_wb_max,
            sum_wind_dir=sum_wind_dir, sum_wind_vel=sum_wind_vel,
            sum_tau_b=sum_tau_b, sum_tau_d=sum_tau_d, sum_clear=sum_clear
        )
    return climate_header
