import os
import glob
import pandas as pd
import xarray as xr
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
from scipy.ndimage import label, generate_binary_structure
from skimage.measure import perimeter
from shapely.geometry import Polygon, Point
from scipy.stats import linregress
from numpy.linalg import eig
import warnings
import gc
import matplotlib as mpl

warnings.filterwarnings("ignore")

class Track_reflectivity_regions:
    
       """
    Identify, track, and analyze radar reflectivity regions (PS) and their interactions with fire polygons.

    Parameters
    ----------
    events_data : pd.DataFrame
        Table with event information (storm start, PFDF coordinates, etc.).
    directory : str
        Path to directory with radar reflectivity NetCDF files (radar mosaics).
    fire_shapefile : str
        Path to shapefile of fire polygons.
    folder_regionplots : str
        Output directory for reflectivity region plots.
    folder_fireplots : str
        Output directory for plots over fire polygons.
    folder_csv : str
        Output directory for property CSV files.
    """

    # ------------------------------------------------------
    # Initialization
    # ------------------------------------------------------
    
    def __init__(self, events_data, directory, fire_shapefile, folder_regionplots, folder_fireplots, folder_csv):
        self.events_data = events_data
        self.directory = directory
        self.fire_shapefile = fire_shapefile
        self.folder_regionplots = folder_regionplots
        self.folder_fireplots = folder_fireplots
        self.folder_csv = folder_csv
        
    # ------------------------------------------------------
    # Utility functions: datetime UTC and storm timing
    # ------------------------------------------------------

    def hour_utc(self,date):
        
    """
    Convert a local (US/Pacific) datetime to a UTC timestamp string.

    Handles string inputs and ambiguous DST times by attempting tz_localize
    with ambiguous=False and falling back to ambiguous=True if needed.

    Returns
    -------
    str
        UTC datetime formatted as '%Y-%m-%d %H:%M'.
    """
        
        if isinstance(date, str):
            date = pd.to_datetime(date)
            
        datetime_local = pd.to_datetime(date)
        
        try:
            
            datetime_local = pd.to_datetime(date).tz_localize("US/Pacific", ambiguous = False)
        
        except pd._libs.tslibs.tz.AmbiguousTimeError:
            
            datetime_local = pd.to_datetime(date).tz_localize("US/Pacific", ambiguous = True)
        
        hour_UTC = datetime_local.tz_convert('UTC')

        hour_UTC_str = hour_UTC.strftime('%Y-%m-%d %H:%M')
    
        return hour_UTC_str
    
    def start_end_storm_df(self,):
        
    """
    Extract the start and end of a storm event.

    """
    
        df = pd.DataFrame()
        event_data=self.events_data
        
        # compute end time from start + duration (hours)
        
        end_storm = pd.to_datetime(event_data.Storm_StartDate) + pd.to_timedelta(event_data.Storm_Duration, unit='H')
        
        #start tracking 10 min before
        
        start_storm_utc = [pd.to_datetime(self.hour_utc(event_data.Storm_StartDate)) - dt.timedelta(minutes=10)]
        end_storm_utc = [pd.to_datetime(self.hour_utc(end_storm))]
        
        print(start_storm_utc, end_storm_utc)
        
        return start_storm_utc, end_storm_utc
    
    def get_reflectivity_value(self, netcdf, latitude, longitude, time):
        
    """
    Return the reflectivity value at the nearest grid point to (lat, lon, time).

    This uses xarray .sel with method='nearest' to find the closest grid cell.
    """
    
        reflectivity = netcdf.sel(time= time, lon = longitude, lat = latitude, method = 'nearest')
        
        return reflectivity.values
    
    def select_files(self, start_timestamp, end_timestamp):
        
     """
    Select NetCDF files whose filename timestamps fall between start and end.

    Expects filenames with pattern containing '_<YYYYmmddHHMM>.nc' (your original parser).
    Attempts to open each file (xr.open_dataset) to skip files that cannot be read.

    Returns
    -------
    (selected_files_sorted, times_sorted)
    """
        start = pd.to_datetime(start_timestamp)
        end = pd.to_datetime(end_timestamp)
        all_files = glob.glob(os.path.join(self.directory, '**', '*.nc'), recursive=True)
        selected_files = []
        times = []
        for file in all_files:
            try:
                filename = os.path.basename(file)
                file_timestamp_str = filename.split('_')[1].split('.')[0]
                file_timestamp = pd.to_datetime(file_timestamp_str, format='%Y%m%d%H%M')
                
                 # try opening to ensure file is readable; skip if open fails
                if start <= file_timestamp <= end:
                    try:
                        with xr.open_dataset(file) as ds:
                            pass  # If opening succeeds, do nothing and move to the next step
                    except Exception as open_error:
                        print(f"Skipping file {file} due to read error: {open_error}")
                        continue  # Skip this file and move to the next one
                
                    selected_files.append(file)
                    times.append(file_timestamp)

            except Exception as e:
                print(f"Error processing file {file}: {e}")
                 
        return sorted(selected_files), sorted(times)

    def get_fire_polygons(self, map_ids):
    """
    Read fire polygons shapefile (relative to home) and return polygons with Map_ID in map_ids.
    """
        
        fire_polygons = gpd.read_file(os.path.join(os.path.expanduser('~'),self.fire_shapefile))
        fire_polygons = fire_polygons.to_crs(epsg=4326)
        polygons = fire_polygons[fire_polygons.Map_ID.isin(map_ids)]
        return polygons
    

    def case_study_domain(self,centroid, casestudy_data):
    """
    Extract a square domain around the fire centroid from the case-study dataset.

    The domain is +/-5 degrees in lat & lon around the centroid.
    """
        centroid_x = centroid.x.item()
        centroid_y = centroid.y.item()
        casestudy_domain = casestudy_data.sel(lat=slice(centroid_y-5, centroid_y+5), lon=slice(centroid_x-5, centroid_x+5))
        return casestudy_domain

    # ------------------------------------------------------
    # Region segmentation and filtering
    # ------------------------------------------------------
    
    def segment_regions(self,data):
        
    """
    Label contiguous regions in a 2D reflectivity field.

    Parameters
    ----------
    data : np.ndarray or xarray.DataArray
        Input reflectivity field, where valid pixels are non-NaN.

    Returns
    -------
    labeled_array : np.ndarray
        2D array with unique integer IDs assigned to each region.
    num_features : int
        Number of detected regions.
    """
        
        binary_data = ~np.isnan(data) # Mask out NaN values
        s = generate_binary_structure(2, 2) # 8-connectivity
        labeled_array, num_features = label(binary_data, structure=s)
        return labeled_array, num_features

    def filter_regions_by_size(self,labeled_array, min_size):
        
    """
    Remove regions smaller than a given size threshold.

    Parameters
    ----------
    labeled_array : np.ndarray
        2D labeled array from `segment_regions`.
    min_size : int
        Minimum number of pixels required for a valid region.

    Returns
    -------
    filtered_labeled_array : np.ndarray
        Labeled array with only valid (large enough) regions retained.
    """
        
        unique, counts = np.unique(labeled_array, return_counts=True)
        filtered_labels = unique[counts >= min_size]
        filtered_labeled_array = np.isin(labeled_array, filtered_labels) * labeled_array
        return filtered_labeled_array

    def region_within_polygon(self,labeled_array, polygon):
        
    """
    Check which region IDs intersect with fire polygon 

    Parameters
    ----------
    labeled_array : xarray.DataArray
        Labeled array of reflectivity regions 
    polygon : geopandas.GeoDataFrame
        Fire perimeter or polygon geometry.

    Returns
    -------
    np.ndarray
        Unique region IDs that overlap with the polygon.
    """
        
        labeled = labeled_array.rio.write_crs("EPSG:4326")
        labeled.rio.set_spatial_dims(x_dim='lon', y_dim='lat', inplace=True)
        polygon_mask = labeled.rio.clip(polygon.geometry)
        polygon_mask = polygon_mask.where(~polygon_mask.isin([-2147483648, 0]))
        return np.unique(polygon_mask)
    
    # -------------------------------
    # Property computation utilities
    # -------------------------------

    def is_split(self,matches):
        # Check if it's a split
        if len(matches) > 1:
            return True
        return False

    def is_merge(self,area_t_i, area_t_i1, corr, prev_corr, region_id_t_i, region_1, prev_areas):
        
    """
    Check whether two regions have merged into one.

    Criteria
    --------
    - Area increases between timesteps.
    - Correlation drops below 95% of previous correlation.
    - Region pair already tracked in `prev_areas`.

    Returns
    -------
    bool
        True if merge occurred, False otherwise.
    """
    
        if area_t_i1 > area_t_i and corr < 0.95 * prev_corr and (region_id_t_i, region_1) in prev_areas:
            return True
        return False
    
    def compute_cg_idl(self, reflectivity_netcdf, all_segmented_arrays, region_id, t):
        
    """
    Compute reflectivity-weighted center of gravity for a region.

    Returns
    -------
    tuple : (lon_cold, lat_cold)
        Weighted center of gravity coordinates.
    """
        
        # Extract the relevant reflectivity data
        dbz_values = reflectivity_netcdf[t, :, :].where(all_segmented_arrays[t] == region_id)

        # Extract latitude and longitude from the dataset
        lat_values = dbz_values.lat.values
        lon_values = dbz_values.lon.values

        # Create 2D arrays for latitude and longitude using meshgrid
        lat_2d, lon_2d = np.meshgrid(lat_values, lon_values, indexing='ij')

        # Filter out NaNs
        valid_mask = ~np.isnan(dbz_values.values)

        # Compute the number of valid pixels
        num_valid_pixels = valid_mask.sum()

        # Compute the unweighted center of gravity
        xcg = lon_2d[valid_mask].sum() / num_valid_pixels
        ycg = lat_2d[valid_mask].sum() / num_valid_pixels

        # Weighted by reflectivity (inverse weighting: stronger dbz = lower weight)
        temp_weights = 90 - dbz_values.values
        valid_mask = ~np.isnan(temp_weights)

        lat_weights = lat_2d[valid_mask] * temp_weights[valid_mask]
        lon_weights = lon_2d[valid_mask] * temp_weights[valid_mask]
        weight_sum = temp_weights[valid_mask].sum()

        lat_cold = lat_weights.sum() / weight_sum
        lon_cold = lon_weights.sum() / weight_sum

        return lon_cold, lat_cold
    
    def compute_cg(self, region_id, t):
        
        #this one can be used too, gave similar values
        
        dbz_values = reflectivity_netcdf[t, :, :].where(all_segmented_arrays[t] == region_id)
        xcg = ((dbz_values.lon * dbz_values).sum() / dbz_values.sum()).item()
        ycg = ((dbz_values.lat * dbz_values).sum() / dbz_values.sum()).item()
        
        return xcg, ycg
    
    def compute_speed_idl(self, cgi, cgi_1, ti, ti_1):
        
    """
    Compute speed and direction from centers of gravity between timesteps.

    Returns
    -------
    tuple : (speed_kmh, angle_deg, delx, dely)
    """
        
        lon_cgi = cgi[0]
        lat_cgi = cgi[1]
        lon_cgi_1 = cgi_1[0]
        lat_cgi_1 = cgi_1[1]

        time_diff = pd.to_datetime(ti_1.values) - pd.to_datetime(ti.values)
        time_diff_hours = time_diff / np.timedelta64(1, 'h')
        print(time_diff_hours)

        r = 6371 # Earth radius (km)

        delfi = (lon_cgi_1 - lon_cgi) * np.pi /180
        thetamed = ((lat_cgi_1 + lat_cgi)/2) * np.pi /180
        delx = r*np.cos(thetamed)*delfi

        dely = r*((lat_cgi_1 - lat_cgi) * np.pi /180)

        yvel = np.sqrt(delx**2 + dely**2)/time_diff_hours
        angle_deg1 = np.degrees(np.arctan2(dely, delx)) % 360
        

        return yvel, angle_deg1, delx, dely
    
    def max_fraction_fire(self, netcdf, polygon, thrhold, t):
        
    """
    Compute maximum reflectivity and fraction of pixels exceeding threshold within fire polygon.

    Returns
    -------
    tuple : (polygon_mask, max_value, fraction_exceeding)
    """
        
        netcdf.rio.set_spatial_dims(x_dim='lon', y_dim='lat', inplace=True)
        netcdf.rio.write_crs("EPSG:4326", inplace=True)
        polygon_mask = netcdf.rio.clip(polygon.geometry, all_touched=True)
        
        dbz_fire = polygon_mask[t,:,:]
        
        max_polygon = dbz_fire.max().values
        
        fraction = (dbz_fire >= thrhold).sum().values / dbz_fire.count().values
        #print(max_polygon, fraction)
        
        return polygon_mask, max_polygon, fraction
        
    def compute_mean_variance_min_max(self, reflectivity_netcdf, all_segmented_arrays, region_id, t):
        
    """
    Compute basic statistics and geometry for a reflectivity region.
    """

        pixel_perimeter = ((reflectivity_netcdf.lat[1] - reflectivity_netcdf.lat[0]) *111)
        pixel_area = (pixel_perimeter)**2

        region = all_segmented_arrays[t] == region_id
        dbz_values = reflectivity_netcdf[t,:,:].where(region)
        dbz_mean = dbz_values.mean()
        dbz_var = dbz_values.var()
        dbz_min = dbz_values.min()
        dbz_max = dbz_values.max()
        region_area = np.sum(region) * pixel_area
        region_perimeter = perimeter(region) * pixel_perimeter

        return dbz_mean.item(), dbz_var.item(), dbz_min.item(), dbz_max.item(), region_perimeter.item(), region_area.item()
    
    def orientation_least_squares(self, reflectivity_netcdf, all_segmented_arrays, region_id, t):
        
    """
    Estimate orientation angle using least-squares linear regression on lat/lon.
    """

        coords = np.column_stack(np.where(all_segmented_arrays[t] == region_id))
        lats = reflectivity_netcdf.lat.values[coords[:, 0]]
        lons = reflectivity_netcdf.lon.values[coords[:, 1]]

        slope, intercept, _, _, _ = linregress(lons, lats)

        angle_rad = np.arccos(-1/(np.sqrt(1+slope**2)))
        #angle_rad = np.arctan(slope)
        angle_deg = np.degrees(angle_rad)

        #if angle_deg < 0:
            #angle_deg += 180

        return angle_deg
    
    def orientation_eof(self, reflectivity_netcdf, all_segmented_arrays, region_id, t):
        
    """
    Estimate orientation and eccentricity using EOF (PCA).
    """
            
        coords = np.column_stack(np.where(all_segmented_arrays[t] == region_id))
        lats = reflectivity_netcdf.lat.values[coords[:, 0]]
        lons = reflectivity_netcdf.lon.values[coords[:, 1]]

        f_mn = np.column_stack((lats, lons))
        covar = np.cov(f_mn, rowvar= False)
        eig_val, eig_vec = eig(covar)
        
        # Sort eigenvalues/vectors
        sorted_indices = np.argsort(eig_val)[::-1]
        eig_val = eig_val[sorted_indices]
        eig_vec = eig_vec[:, sorted_indices]

        # orientation of 1st EOF
        first_eigenvector = eig_vec[:,0]
        angle = np.arccos(-1*first_eigenvector[0])  

        # Determine the quadrant
        #if first_eigenvector[1] < 0:
            #angle = -angle

        angle_deg = np.degrees(angle)

        #norm_eof1 = np.linalg.norm(eig_vec[:, 0])  
        #norm_eof2 = np.linalg.norm(eig_vec[:, 1])  

        eccentricity = np.sqrt(eig_val[1] / eig_val[0])

        #eccentricity = norm_eof2 / norm_eof1

        return angle_deg, eccentricity
    
    def track_and_compute_properties(self, reflectivity_netcdf, reflectivity_case, all_segmented_arrays, label_regions, time_netcdf, corr_threshold, 
                                     polygon, strong_rainfall_thrsh, events_data):

     
    """
    Tracks reflectivity regions over time and computes multiple physical properties.

    Parameters
    ----------
    reflectivity_netcdf : xarray.DataArray
        NetCDF array containing reflectivity data (subsetted if necessary, e.g., >20 dBZ).
    reflectivity_case : xarray.DataArray
        Full reflectivity data without thresholding, used for additional property extraction.
    all_segmented_arrays : list of np.ndarray
        List of 2D arrays containing labeled regions for each timestep.
    label_regions : list of np.ndarray
        List of labeled regions (integer IDs) for tracking at each timestep.
    time_netcdf : list of datetime-like
        Timestamps corresponding to each timestep in the NetCDF arrays.
    corr_threshold : float
        Minimum correlation required to consider a region matched between consecutive timesteps. e.g. 20 dBZ
    polygon : GeoDataFrame
        Polygon(s) defining a region of interest (e.g., burned area) for fraction/max calculations.
    strong_rainfall_thrsh : float
        Threshold of reflectivity (dBZ) to compute fraction of strong rainfall within polygon. e.g. 30 dBZ
    events_data : DataFrame
        Event data with latitude/longitude for extracting reflectivity at known points.

    Returns
    -------
    reflectivity_fire : xarray.DataArray
        Reflectivity values clipped to the polygon
    matched_regions : list of tuples
        Each tuple corresponds to a timestep and contains a list of matched region pairs:
        (region_id_t, region_id_t+1, correlation).
    properties_df : pandas.DataFrame
        Table of computed properties for each tracked region, including:
        - Center of gravity (unweighted & weighted by reflectivity)
        - Speed and direction
        - Orientation (least squares and EOF)
        - Reflectivity statistics (mean, variance, min, max)
        - Region perimeter and area
        - Fraction and max reflectivity within polygon
        - Event type (split, merge, fast, undetermined)
    """
    
        # Initialize lists and dictionaries to store matches, splits/merges, previous area/correlation
        matched_regions = [] # Stores matched regions per timestep
        #splits = []
        #merges = []
        prev_areas = {}  # Dictionary to store previous area values   
        prev_corrs = {}  # Dictionary to store previous correlation values

        properties = []  # List to store computed properties

        # Loop through each consecutive pair of timesteps
        for i in range(len(label_regions) - 1):
            regions_t_i = label_regions[i]
            regions_t_i1 = label_regions[i + 1]

            matched_regions_at_t = []

            # Loop through each region in current timestep
            for region_id_t_i in np.unique(regions_t_i):
                # Extract reflectivity values for current region
                region_t_i = reflectivity_netcdf[i, :, :].where(all_segmented_arrays[i] == region_id_t_i)

                matches = []

                # Compare to each region in the next timestep
                for region_id_t_i1 in np.unique(regions_t_i1):
                    region_t_i1 = reflectivity_netcdf[i + 1, :, :].where(all_segmented_arrays[i + 1] == region_id_t_i1)

                    # Compute correlation between regions
                    corr = xr.corr(region_t_i, region_t_i1).item()

                    # If correlation exceeds threshold, consider regions matched
                    if corr > corr_threshold:
                        matched_corr = corr
                        matches.append((region_id_t_i1, matched_corr))

                # If any matches found, compute properties for each match       
                if matches:
                    for match in matches:
                        print(i, region_id_t_i, match[0], matched_corr)
                        matched_regions_at_t.append((region_id_t_i, match[0], match[1]))

                        # Compute physical properties for current and next timestep
                        cg_t_i = self.compute_cg_idl(reflectivity_netcdf, all_segmented_arrays, region_id_t_i, i)
                        cg_t_i1 = self.compute_cg_idl(reflectivity_netcdf, all_segmented_arrays, match[0], i + 1)
                        speed_idl, direction, delx, dely = self.compute_speed_idl(cg_t_i, cg_t_i1, time_netcdf[i], time_netcdf[i + 1])
                        #speed = compute_speed(cg_t_i, cg_t_i1, time_netcdf[i], time_netcdf[i + 1])
                        orient_lsq = self.orientation_least_squares(reflectivity_netcdf, all_segmented_arrays, match[0], i + 1)
                        orient_eof_angle, eccentricities = self.orientation_eof(reflectivity_netcdf, all_segmented_arrays, match[0], i + 1)
                        reflectivity_df = self.get_reflectivity_value(reflectivity_case, events_data.LAT, events_data.LONG, 
                                                                 time_netcdf[i + 1].values)
                        mean, var, min_val, max_val, perim, area_t_i1 = self.compute_mean_variance_min_max(
                            reflectivity_netcdf, all_segmented_arrays, match[0], i + 1
                        )
                        reflectivity_fire, max_fire, fraction = self.max_fraction_fire(reflectivity_case, polygon, strong_rainfall_thrsh, i + 1)

                        print(speed_idl, direction, area_t_i1)

                        # Determine event type (split/merge/fast/undetermined)
                        if (i, region_id_t_i) in prev_areas:
                            area_t_i = prev_areas[(i, region_id_t_i)]
                            prev_corr = prev_corrs[(i, region_id_t_i)]

                            # Initial classification
                            split_flag = self.is_split(matches)
                            merge_flag = self.is_merge(area_t_i, area_t_i1, matched_corr, prev_corr, region_id_t_i, match[0], prev_areas)

                            # Determine event type
                            if speed_idl > 50:
                                if area_t_i > area_t_i1:
                                    event_type = 'fast_split'
                                elif area_t_i < area_t_i1:
                                    event_type = 'fast_merge'
                                elif area_t_i == area_t_i1:
                                    event_type = 'fast_noclear'
                            else:
                                if split_flag:
                                    event_type = 'split'
                                elif merge_flag:
                                    event_type = 'merge'
                                else:
                                    event_type = 'no split/merge'

                        else:
                            # If no previous area information, default to 'undetermined'
                            event_type = 'undetermined'

                        #print(event_type)
                        # Store current region area and correlation for next timestep
                        prev_areas[(i + 1, match[0])] = area_t_i1
                        prev_corrs[(i + 1, match[0])] = matched_corr

                        # Add computed properties to the list
                        properties.append({ 
                            'datetime_t': time_netcdf[i].values,
                            'time_t': i,
                            'region_id_t': region_id_t_i,
                            'datetime_t1': time_netcdf[i + 1].values,
                            'time_t1': i + 1,
                            'region_id_t1': match[0],
                            'corr': matched_corr,
                            'cg_t_i': cg_t_i,
                            'cg_t_i1': cg_t_i1,
                            'eccs': eccentricities,
                            'delx': delx,
                            'dely' : dely,
                            'speed_km_per_h':speed_idl,
                            'orient_lsq': orient_lsq,
                            'orient_eof': orient_eof_angle,
                            'direction':direction,
                            'mean': mean,
                            'variance': var,
                            'min_reg_val': min_val,
                            'max_reg_val': max_val,
                            'perimeter': perim,
                            'area': area_t_i1,
                            'reflectivity_df_t1': reflectivity_df,
                            'max_fire_val': max_fire,
                            'fraction_thrsh': fraction,
                            'event_type': event_type  
                        })

                  
            # Save matches for this timestep
            matched_regions.append((i + 1, matched_regions_at_t))
            
        # Convert list of dictionaries to DataFrame for easy analysis
        properties_df = pd.DataFrame(properties)

        return reflectivity_fire, matched_regions, properties_df
    
    def obtain_label_regions(self, reflectivity_netcdf, reflectivity_netcdf_time, reflectivity_netcdf_lat, reflectivity_netcdf_lon, polygon):
        
    """
    Identifies and tracks labeled reflectivity regions within a polygon (e.g., burned area) for each timestep.

    Keeps all segmented arrays in `all_segmented_arrays` for downstream computations.

    Parameters
    ----------
    reflectivity_netcdf : xarray.DataArray
        3D array (time, lat, lon) of reflectivity values.
    reflectivity_netcdf_time : array-like
        List/array of timestamps corresponding to the reflectivity data.
    reflectivity_netcdf_lat, reflectivity_netcdf_lon : array-like
        Latitude and longitude coordinates of the reflectivity grid.
    polygon : GeoDataFrame
        Polygon defining the area of interest.

    Returns
    -------
    label_regions : dict
        Dictionary with timestep as key and array of region IDs inside the polygon as values.
    all_segmented_arrays : list of np.ndarray
        List of 2D arrays for each timestep containing the labeled regions (including small ones).
    """
    
        label_regions = {}

        all_segmented_arrays = []

        for t in range(len(reflectivity_netcdf_time)):

            #Identify regions that have adjacent True values
            segmented, num_features = self.segment_regions(reflectivity_netcdf[t,:,:])
            #Save arrays with pixels corresponding to the label region
            all_segmented_arrays.append(segmented)
            #For each array, remove regions less than 180 pixels
            segmented = self.filter_regions_by_size(segmented, 180)
            # Convert arrays to xarrays for the identification of regions within fire polygon
            segmented = xr.DataArray(segmented,coords=[reflectivity_netcdf_lat, reflectivity_netcdf_lon], dims=['lat', 'lon'])
            # Extract regions within fire polygon
            regions_t_i = self.region_within_polygon(segmented, polygon)
            # Track only regions that are found in fire polygon 
            label_regions[t] = regions_t_i[~np.isnan(regions_t_i)]

            #print(t, regions_t_i, reflectivity_netcdf_time[t])

        return label_regions, all_segmented_arrays
    
    def format_dfpoints(self,lat_col='LAT', long_col='LONG', crs='EPSG:4326'):
    """
    Converts a DataFrame or Series of event points into a GeoDataFrame.

    Parameters
    ----------
    lat_col, long_col : str
        Names of the latitude and longitude columns in the DataFrame.
    crs : str
        Coordinate reference system for the output GeoDataFrame.

    Returns
    -------
    gdf : GeoDataFrame
        GeoDataFrame with event points as geometry.
    """
        
        #Convert a single row Series to a DataFrame if necessary.
        df_or_row = self.events_data
        if isinstance(df_or_row, pd.DataFrame):
            df = df_or_row
        elif isinstance(df_or_row, pd.Series):
            # Convert a single row to a DataFrame
            df = pd.DataFrame([df_or_row])
        else:
            raise TypeError("Input must be a pandas DataFrame or Series")
        
        # Check if latitude and longitude columns contain lists
        if df[lat_col].apply(lambda x: isinstance(x, list)).any():
            # Explode lists into separate rows
            df_exploded = df.explode([lat_col, long_col])
        else:
            df_exploded = df

        # Create geometry from latitude and longitude
        df_exploded['geometry'] = gpd.points_from_xy(df_exploded[long_col], df_exploded[lat_col])

        # Convert to GeoDataFrame
        gdf = gpd.GeoDataFrame(df_exploded, geometry='geometry', crs=crs)

        return gdf
    
    def process_reflectivity_for_plot(self, event_peak15time, reflectivity_data, time_offset='1H', resample_interval='10T'):
        """
        Extracts and resamples reflectivity data around a specific event time for plotting or analysis.

            Parameters
        ----------
        event_peak15time : str or pd.Timestamp
            Event peak time to center the window on.
        reflectivity_data : xarray.DataArray or xarray.Dataset
            Reflectivity dataset to extract from.
        time_offset : str
            Time window around the event (e.g., '1H' for ±1 hour).
        resample_interval : str
            Resampling interval for the data (e.g., '10T' for 10 minutes).

        Returns
        -------
        reflectivity_resampled : xarray.DataArray or xarray.Dataset
            Reflectivity data for the selected time window, resampled at the specified interval.
        """
        # Convert event time to datetime if not already
        event_time = pd.to_datetime(event_peak15time)

        # Define a time window around the event (default ±1 hour).
        start_time = event_time - pd.Timedelta(time_offset)
        end_time = event_time + pd.Timedelta(time_offset)

        # Select the time window from the reflectivity data
        reflectivity_window = reflectivity_data.sel(time=slice(start_time, end_time))

        # Resample the data at a defined interval (default 10 minutes).
        reflectivity_resampled = reflectivity_window.resample(time=resample_interval).nearest()

        return reflectivity_resampled
    
    def plot_regions_dbz(self, df_properties_timestep, reflectivity_array, reflectivity_fire, all_segmented_arrays, debrisflows, fire_polygon, highlight_time=None, n_cols=3, firename=None, iteration_idx=None):
        
        #Define NWS reflectivity colormap.
        nws_reflectivity_colors = [ "#646464",  "#04e9e7", "#019ff4",  "#0300f4", "#02fd02",  "#01c501",  "#008e00", 
            "#fdf802", "#e5bc00", "#fd9500", "#fd0000", "#d40000", "#bc0000",  "#f800fd",  "#9854c6", "#fdfdfd" # 75
            ]
        
        def radar_colormap():
            return mpl.colors.ListedColormap(nws_reflectivity_colors)
        """
        Plots reflectivity regions and properties in a grid of subplots.

        Parameters
        ----------
        df_properties_timestep : pd.DataFrame
            DataFrame containing properties of tracked regions for a timestep.
        reflectivity_array : xarray.DataArray
            Full reflectivity dataset.
        reflectivity_fire : xarray.DataArray
            Reflectivity clipped to fire polygon.
        all_segmented_arrays : list of np.ndarray
            Segmented region arrays per timestep.
        debrisflows : GeoDataFrame
            Locations of debris flows for overlay.
        fire_polygon : GeoDataFrame
            Fire perimeter polygon for overlay.
        highlight_time : datetime-like, optional
            Time to highlight on plots and in tracking CSV.
        n_cols : int
            Number of columns in the grid of subplots.
        firename : str, optional
            Name of the fire for tracking CSV logging.
        iteration_idx : int, optional
            Iteration index for tracking CSV logging.

        Returns
        -------
        fig1, fig2 : matplotlib.figure.Figure
            Figure objects for the segmented regions and fire reflectivity plots.
        """
        
        cmap = radar_colormap()
        norm = mpl.colors.Normalize(vmin=0, vmax=80)
    
        n_plots = len(df_properties_timestep)
        n_rows = (n_plots + n_cols - 1) // n_cols  # Calculate number of rows needed
        
        #Create two sets of subplots: Fig1: reflectivity segmented by tracked regions, Fig2: reflectivity within fire polygon
        
        fig1, axes1 = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 5), subplot_kw={'projection': ccrs.PlateCarree()})
        fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 5), subplot_kw={'projection': ccrs.PlateCarree()})

        # Flatten the axes array for easy iteration
        #axes = axes.flatten()
        axes1 = axes1.flatten()
        axes2 = axes2.flatten()

        df_properties_timestep['datetime_t1'] = pd.to_datetime(df_properties_timestep['datetime_t1'])
        
        # Optionally highlight the nearest timestep to a given time and store tracking info.
        tracking_df = pd.DataFrame(columns=['firename','iteration_idx','peak15min', 'nearest_time','corr','diff_exceeds_10min'])

        if highlight_time is not None:

            highlight_time = pd.to_datetime(highlight_time)
            times = df_properties_timestep['datetime_t1'].values
            
            highlight_idx = (abs(times - np.datetime64(highlight_time))).argmin()
            nearest_time = df_properties_timestep.iloc[highlight_idx]['datetime_t1']
            corr_atnearest_time = df_properties_timestep.iloc[highlight_idx]['corr']
            
            time_diff = abs((nearest_time - highlight_time).total_seconds()) / 60
            
            exceeds_10min = "Yes" if time_diff > 10 else "No"
            
            tracking_entry = pd.DataFrame({
            'firename': [firename],  # Passing the firename here
            'iteration_idx': [iteration_idx],
            'peak15min': [highlight_time],
            'nearest_time': [nearest_time],
            'corr': [corr_atnearest_time],
            'diff_exceeds_10min': [exceeds_10min]
            })
        
            tracking_df = pd.concat([tracking_df, tracking_entry], ignore_index=True)

        else:

            highlight_idx = None
            
            #print(highlight_idx)
            
        tracking_df.to_csv(os.path.join(os.path.expanduser('~'),'DATA/Project3/NEXRAD/peak15min_inpropertiesdf_2.csv'), mode='a', 
                           header=not os.path.exists(os.path.join(os.path.expanduser('~'),'DATA/Project3/NEXRAD/peak15min_inpropertiesdf_2.csv')), 
                           index=False)
        
        #Loop through each row in `df_properties_timestep`:
        for i, (idx, time_row) in enumerate(df_properties_timestep.iterrows()):
            print(i)
            
            time_datetime = time_row['datetime_t1']
            time_step = time_row['time_t1']
            region_id = time_row['region_id_t1']
            speed = time_row['speed_km_per_h']
            eccentricity = time_row['eccs']
            direction = time_row['direction']
            reflectivity = time_row['reflectivity_df_t1']
            
            ax1 = axes1[i]  # Axes for first figure (segmented regions)
            ax2 = axes2[i]  # Axes for second figure (only reflectivity)

            
            # Plot the reflectivity array for the region in Fig1.
            img1 = reflectivity_array[time_step,:,:].where(all_segmented_arrays[time_step]==region_id).plot(ax=ax1, cmap=cmap, norm=norm, add_colorbar=False)
            
            #Overlay debris flow points and fire polygon.
            debrisflows.plot(ax=ax1, markersize=5, color='black', label='Debris flow')
            fire_polygon.plot(ax=ax1, edgecolor='red', facecolor='none', label="Fire Polygon")
            cbar = plt.colorbar(img1, label = 'Reflectivity [dBZ]')
            cbar.ax.set_ylim(20,80)
            ax1.coastlines()
           
            
            # Add a text box with speed, eccentricity, and direction.
            textstr = (f'Speed: {speed:.1f} km/h\n'
                       f'Eccentricity: {eccentricity:.2f}\n'
                       f'Direction: {direction:.2f}\u00b0')
            props = dict(boxstyle='round', facecolor='white', alpha=0.7)
            ax1.text(0.05, 0.05, textstr, transform=ax1.transAxes, fontsize=10,
                    verticalalignment='bottom', horizontalalignment='left', bbox=props)
            
            # Optional: plot center of gravity if needed
            cg = time_row['cg_t_i1']
            ax1.plot(cg[0], cg[1], 'o', markersize=5, color='gray', label = 'Center of gravity')
            
            
            # Gridlines
            gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=.6, color='gray', alpha=0.5, linestyle='-.')
            gl.xlabel_style = {"size": 7}
            gl.ylabel_style = {"size": 7}
            gl.top_labels = False
            gl.right_labels = False
            
            #Repeat plotting for Fig2 using `reflectivity_fire`.
            #Check dBZ threshold here
            reflectivity_fire = reflectivity_fire.where(reflectivity_fire >=20)
            img2 = reflectivity_fire[time_step, :, :].plot(ax=ax2, cmap=cmap, norm=norm, add_colorbar = False)
            debrisflows.plot(ax=ax2, markersize=5, color='black', label='Debris flow')
            fire_polygon.plot(ax=ax2, edgecolor='red', facecolor='none')
            cbar = plt.colorbar(img2, label = 'Reflectivity [dBZ]', shrink=0.75)
            cbar.ax.set_ylim(20,80)
            ax2.coastlines()
            ax2.set_title(f"{time_datetime} UTC")
            
            
            # Gridlines
            gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=.6, color='gray', alpha=0.5, linestyle='-.')
            gl.xlabel_style = {"size": 7}
            gl.ylabel_style = {"size": 7}
            gl.top_labels = False
            gl.right_labels = False


            if highlight_time is not None and i == highlight_idx:
                ax1.spines['geo'].set_linewidth(2)
                ax2.spines['geo'].set_linewidth(2)

        # Hide unused subplots if there are fewer plots than grid cells
        for j in range(i + 1, len(axes1)):
            axes1[j].set_visible(False)
            axes2[j].set_visible(False)

        fig1.tight_layout()
        fig2.tight_layout()
        
        return fig1, fig2
        
        
def run(events_data, directory, fire_shapefile, reflectivity_threshold, corr_threshold, strong_rainfall_thrsh, folder_regionplots, folder_fireplots, folder_csv):
    
     """
    Main workflow to track reflectivity regions associated with a fire event and compute properties.

    Parameters
    ----------
    events_data : pandas Series
        Single row from the events DataFrame containing fire and debris flow info.
    directory : str
        Directory containing reflectivity NetCDF files (radar mosaics).
    fire_shapefile : str
        Path to shapefile with fire polygons.
    reflectivity_threshold : float
        Threshold for reflectivity to identify regions (dBZ).
    corr_threshold : float
        Correlation threshold for region tracking.
    strong_rainfall_thrsh : float
        Threshold to flag strong rainfall regions.
    folder_regionplots : str
        Directory to save region-based plots.
    folder_fireplots : str
        Directory to save fire polygon-based plots.
    folder_csv : str
        Directory to save properties CSV files.

    Returns
    -------
    None
        Saves CSV and plots to disk. No objects are returned.
    """
    #Initialize the Track_reflectivity_regions class with event data, directories, and shapefile.
    track = Track_reflectivity_regions(events_data, directory, fire_shapefile, folder_regionplots, folder_fireplots, folder_csv)
    
    #    2. Extract event information: FireMap_ID, latitude, longitude, peak intensity time (UTC), and fire name.
    
    firemap_id = events_data.FireMap_ID
    event_lat = events_data.LAT
    event_lon = events_data.LONG
    peak_intimeutc = track.hour_utc(events_data.Peak_I15_Date)
    fire_name = events_data.FireName
    
    #Format debris flow locations as a GeoDataFrame.
    debris_flows = track.format_dfpoints()
    
    print(peak_intimeutc, events_data.Peak_I15_Date, fire_name, debris_flows.geometry) 
    
    # Generate filenames for plots and CSV output based on fire name, time, and latitude.
    plot_filename = f"{fire_name}_{pd.to_datetime(peak_intimeutc).strftime('%Y%m%d_%H%M')}_{str(event_lat).replace('.', '_')}.png"
    csv_filename = f"{fire_name}_{pd.to_datetime(peak_intimeutc).strftime('%Y%m%d_%H%M')}_{str(event_lat).replace('.', '_')}.csv"
    #print(plot_filename, csv_filename)

    # Determine start and end timestamps of the storm.
    start_timestamp, end_timestamp = track.start_end_storm_df()
    
    #Select reflectivity files for the storm period and load them into an xarray dataset.
    selected_files, times = track.select_files(start_timestamp, end_timestamp)
    casestudy_data = xr.open_mfdataset(sorted(selected_files),combine='nested', concat_dim='time', coords='minimal').sel(lon=slice(-124, -114), lat=slice(32,44)) #close to California boundaries
    casestudy_data['time'] = times
    #print(casestudy_data)
    
    # Clip the reflectivity data to a study domain around the fire centroid.
    fire_polygons = track.get_fire_polygons([firemap_id])
    fire_centroid = fire_polygons.centroid
    #print(fire_centroid )
    reflectivity_casestudy = track.case_study_domain(fire_centroid, casestudy_data)
    
    # Threshold reflectivity data (e.g., > 20 dBZ) to isolate reflectivity regions.
    
    reflectivity_gt_th = reflectivity_casestudy.where(reflectivity_casestudy.composite_n0q >= reflectivity_threshold, drop=False)
    reflectivity_gt_th.load()
    
    # Identify labeled regions within the fire polygon for each timestep.
    label_regions, all_segmented_arrays = track.obtain_label_regions(reflectivity_gt_th.composite_n0q, reflectivity_gt_th.time, reflectivity_gt_th.lat, 
                                                               reflectivity_gt_th.lon, fire_polygons)
    print('label regions completed')
    
    #segmented_final = np.stack(all_segmented_arrays, dtype=np.int16)
     
    # Track the labeled regions over time and compute properties: Returns reflectivity at fire, correlation metrics, and a properties DataFrame.
    reflectivity_at_fire, regions_corr, properties_df = track.track_and_compute_properties(reflectivity_gt_th.composite_n0q, reflectivity_casestudy.composite_n0q, all_segmented_arrays, label_regions, reflectivity_gt_th.time, corr_threshold, fire_polygons, strong_rainfall_thrsh, events_data)
    print('tracking completed')
    
    #Save the properties DataFrame to CSV.
    if properties_df.empty:
        print("DataFrame is empty, skipping this iteration.")
        return 
    properties_df.to_csv(os.path.join(os.path.expanduser('~'),folder_csv, csv_filename))
    
    # Resample reflectivity data for plotting around the peak intensity time.
    resampled_reflectivity = track.process_reflectivity_for_plot(peak_intimeutc, reflectivity_gt_th)
    
    #Select properties corresponding to the resampled time steps.
    properties_timestep = properties_df[properties_df.datetime_t1.isin(resampled_reflectivity.time.values)]
    
    if properties_timestep.empty:
        print("DataFrame with peak intensity timestep is empty, no plot saved and skipping this iteration.")
        return 
    
    #Plot the reflectivity regions and properties in two figures: Segmented region plot and fire polygon plot.
    fig1, fig2 = track.plot_regions_dbz(properties_timestep, reflectivity_gt_th.composite_n0q, reflectivity_at_fire, all_segmented_arrays, debris_flows, fire_polygons, highlight_time=peak_intimeutc, n_cols=3, firename=fire_name, iteration_idx=index)
    
    #Save the figures to the respective folders.
    fig1.savefig(os.path.join(os.path.expanduser('~'),folder_regionplots, plot_filename))
    fig2.savefig(os.path.join(os.path.expanduser('~'),folder_fireplots, plot_filename))
    plt.clf()  # Clear the current figure
    plt.close()
    
    #Clean up memory by deleting large objects and calling garbage collection.
    del track, casestudy_data, fire_polygons, reflectivity_casestudy, reflectivity_gt_th, label_regions, all_segmented_arrays, regions_corr, properties_df, resampled_reflectivity, reflectivity_at_fire, fig1, fig2,
    properties_timestep
    gc.collect()
    
        
if __name__ == "__main__":
    
    # -------------------------------------------------------------------------
    # Load event data
    # -------------------------------------------------------------------------
    # Example: read a CSV with storm/debris-flow events. If want to run only select specific rows, use .iloc
    # dfs = pd.DataFrame(pd.read_csv("~/DATA/Project3/Non_PFDFStorms_95th_modified.csv")).iloc[374:]
    
    # Currently using a CSV with CAPE, IVT, IWV (atmospheric variables) and selecting row 122
    
    dfs = pd.DataFrame(pd.read_csv("~/DATA/cape_ivt_iwv.csv"))

    # -------------------------------------------------------------------------
    # Define directories and shapefile paths
    # -------------------------------------------------------------------------
    directory = '/home/silvana/DATA/Project3/NEXRAD/all_dates/'
    fire_shapefile = 'DATA/Project2/LandslidesJan2023/Fires_NewInventory.shp'
    folder_csv = 'DATA/Project3/NEXRAD/cases_propertiesdf_trial'
    folder_regionplots = 'DATA/Project3/NEXRAD/cases_plots_trial'
    folder_fireplots = 'DATA/Project3/NEXRAD/cases_fireplots_trial'
    
    # -------------------------------------------------------------------------
    # Loop over each event row in the DataFrame
    # -------------------------------------------------------------------------
    for index, row in dfs.iterrows():
        #print(index,row.LAT)
        
        ## Call the main run function for this event
        run(row, directory, fire_shapefile, 
            reflectivity_threshold=20, 
            corr_threshold=0.3, 
            strong_rainfall_thrsh=30, 
            folder_regionplots, folder_fireplots, folder_csv)
        
    # ---------------------------------------------------------------------
    # Notes / reminders for future runs:
    # ---------------------------------------------------------------------
    # 1. Check which column contains peak 15-min rainfall:
    #    - Could be `Peak_I15_Date`
    # 2. Ensure output directories exist and are writable:
    #    - folder_csv, folder_regionplots, folder_fireplots
    # 3. Verify reflectivity threshold (reflectivity_threshold) and correlation threshold (corr_threshold)
    # 4. Latitude/longitude bounds for California are hard-coded in `run()`:
    #    - lon=slice(-124, -114), lat=slice(32,44)
    #    - Adjust if using other regions
    # 5. The fire polygon box around the fire centroid is used in `case_study_domain()`
    # ---------------------------------------------------------------------
