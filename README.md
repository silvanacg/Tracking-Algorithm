# Post-Fire Debris Flow (PFDF) Dataset and Tracking Algorithm (modified MASCOTTE) Code

## Overview

This repository contains a dataset of post-fire debris flow events in California, along with a Python script for tracking reflectivity regions over time. 

## Contents

1. Peak15min_properties.csv

Contains structural properties and reflectivity metrics measured at the peak 15-minute rainfall intensity for each storm (e.g., area, speed, eccentricity, reflectivity metrics at debris flow locations and burned area.)

Includes relevant information about the post-fire debris flow inventory (e.g., event timing, location, rainfall intensity).

2. tracking_code.py

Implements a modified MASCOTTE algorithm to track radar reflectivity regions over time.

Regions are identified and matched across timesteps using a spatial correlation threshold.

Calculates both large-scale storm properties and localized metrics within burned areas.

## Post-Fire Debris Flow Event Dataset – Column Definitions

Repository – Source or citation for the dataset or observation, e.g., USGS, NASA, or a published paper.

EVENT_DATE – Date when the debris flow or landslide event occurred.

NEAREST_PL – Name of the nearest place, community, or landmark to the debris flow event.

INFO_SOURC – Organization or person reporting the event (e.g., USGS, California Geological Survey).

IMPACTS – Type of impacts reported (e.g., Multiple, Road-unpaved, Boulders, Building-residential).

DAMAGE_DES – Text description of observed damages or effects of the debris flow event.

LAT – Latitude of the debris flow event location.

LONG – Longitude of the debris flow event location.

GlobalID – Unique identifier for the event in the source database.

PHOTO_CRED – Link or credit for associated photographs documenting the event.

Acres – Fire perimeter area burned in acres (if associated with wildfire).

FireName – Name of the wildfire associated with the debris flow, if applicable.

StartYear / StartMonth / StartDay – Year, month, and day when the fire started (if applicable).

StartDate – Full date when the fire started.

RainGauge – Name or identifier of the rain gauge nearest to the debris flow location.

Station ID – Unique ID of the rain gauge or meteorological station.

LAT_Rg / LON_Rg – Latitude and longitude of the rain gauge or station.

Distance_km – Distance from the debris flow site to the rain gauge in kilometers.

FireMap_ID – Identifier for the wildfire perimeter polygon in GIS datasets.

Storm_StartDate – Date and time when the associated storm began.

Storm_Duration – Duration of the storm in hours.

Storm_Accum – Total accumulated rainfall during the storm (mm).

Storm_AvgIntensity – Average rainfall intensity over the storm duration (mm/h).

Peak_I15_mm/h / Peak_I30_mm/h / Peak_I60_mm/h – Peak 15-, 30-, and 60-minute rainfall intensities (mm/h).

Peak_I15_Date / Peak_I30_Date / Peak_I60_Date – Date and time when the peak 15-, 30-, and 60-minute rainfall intensities occurred.

CummRain_24h / CummRain_48h – Cumulative rainfall over the 24 and 48 hours preceding the debris flow (mm).

Certainty – Confidence level in the event occurrence (e.g., Certain, Probable).

reg_datetime_t1 – Peak 15-min intensity associated with the debris flow used in radar data.

eccentricity – Eccentricity of the region or PS at peak 15-min intensity(dimensionless, from radar analysis).

area_km – Area of the region or PS at peak 15-min intensity in square kilometers.

speed_km_h – Propagation speed of the region or PS at peak 15-min intensity in km/h.

PSmean_dbz / PSmin_dbz / PSmax_dbz – Mean, minimum, and maximum radar reflectivity (dBZ) of the region or PS peak 15-min intensity

reflect_pfdf – Radar reflectivity observed at the debris flow location at peak 15-min intensity

maxfiredbz_peak15 – Maximum radar reflectivity over the burned area at peak 15-min intensity

direction – Propagation direction of the region or PS at peak 15-min intensity in degrees

percentile_dbzval – Percentile rank of reflectivity at debris flow location at peak 15-min intensity (reflect_pfdf) relative to all values in storm.

percentile_firemax – Percentile rank of maximum reflectivity within burned area at peak 15-min intensity (maxfiredbz) relative to all values in storm.

Season – Season when the debris flow occurred (e.g., Fall/Winter, Spring/Summer).

fract10-15, fract15-20 … fract>70 – Fraction of burned area pixels within specified y radar reflectivity thresholds (dBZ). Each column corresponds to a dBZ range (e.g., 10–15 dBZ, 15–20 dBZ, …, >70 dBZ).
