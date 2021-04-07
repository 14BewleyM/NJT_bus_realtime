# Creating some procedures for bus location data collected from NJ Transit.

# %% notes
# using gtfs_functions module to work with gtfs: https://medium.com/swlh/python-for-transit-get-the-most-out-of-your-gtfs-with-this-python-package-44d0b732f657
# gtfs_functions documentation: https://pypi.org/project/gtfs-functions/#installation
# Partridge (by Remix) also may be helpful for GTFS processing: https://github.com/remix/partridge

# mapping linestrings with folium: https://anitagraser.com/2019/10/31/interactive-plots-for-geopandas-geodataframe-of-linestrings/

# some other useful tutorials
# https://geoscripting-wur.github.io/PythonVector/
# https://www.earthdatascience.org/tags/time-series/

# managing memory with chunking: https://towardsdatascience.com/loading-large-datasets-in-pandas-11bdddd36f7b

# %% TODO
# use resample to select measurements every minute or so: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.resample.html
# also here: https://towardsdatascience.com/using-the-pandas-resample-function-a231144194c4

# CTA uses the same vendor for their bus tracker! and they have an API even!!: https://www.transitchicago.com/developers/bustracker/
# found via github of Neil Freeman, who created an archiver for the API: https://github.com/fitnr/cta-bus-archive

# visualization example: https://twitter.com/fitnr/status/1091860127555219456

# nyc bus dashboard speed calculation: https://github.com/Bus-Data-NYC/nyc-bus-stats/blob/master/sql/speed.sql
# paper with some technical notes on calculating speeds for transit vehicles using AVL: https://sci-hub.do/10.1109/ictis.2017.8047744

# movingpandas for calculating speeds and things: https://anitagraser.com/movingpandas/

# %% modules
import pandas as pd
import geopandas as gpd
import numpy as np
import logging
from datetime import datetime
import gtfs_functions as gtfs
import webbrowser
import folium
import os

# %% some initial settings
generate_maps = 0

# %% directories and filepaths
gtfs_directory = "C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs"
data_path = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/testdata.csv"

# %% set up logging
logging.basicConfig(filename="log.log", 
                    level=logging.INFO) # set to INFO if you want info messages to log

logging.info(f"Logging started at {datetime.now()}")

# %% load data and setup dataframe
logging.info(f"Reading in data from {data_path}")
# buspositions = pd.read_csv(data_path)
buspositions = gpd.GeoDataFrame(pd.read_csv(data_path), 
                      geometry=gpd.points_from_xy(pd.read_csv(data_path).lat, pd.read_csv(data_path).lon))
buspositions.crs = {"init": "epsg:4326"}

# set timestamp column to be datetime       
buspositions.timestamp = pd.to_datetime(buspositions.timestamp)

# set other data types
buspositions.vehicle_id = buspositions.vehicle_id.astype(str).str.split(".").str[0]
buspositions.route_number = buspositions.route_number.astype(str).str.split(".").str[0] 

# drop observations for routes that returned no service
print(f"Dropping {buspositions.shape[0] - buspositions.has_service.sum()} observations from routes with no service")
buspositions = buspositions[buspositions.has_service==True]

# some headsigns had ampersand replaced with &amp;amp;
# put ampersand back
buspositions.headsign = buspositions.headsign.str.replace("amp;amp;", "")

# check for duplicate coordinates in reduced dataframe, to get a sense of the potential problem
# mark the first duplicate, assuming that every duplicate afterwards erroneously suggests vehicle isn't moving
duplicate_count = buspositions.duplicated(subset="geometry", keep="first").sum()
print(f"There are {duplicate_count} observations with duplicate geometries")
# wow, 942297 duplicates
# BUT at least some of these are legit (eg different vehicles on different runs somehow reporting exactly the same coordinates)
# so should subset by vehicle id as well as geometry
duplicate_count = buspositions.duplicated(subset=["vehicle_id", "geometry"], keep="first").sum()
print(f"There are {duplicate_count} observations with duplicate vehicle ids and geometries")

# dataframe for examining all duplicates next to each other
#buspositions[buspositions.duplicated(subset=["vehicle_id", "geometry"], keep=False)].sort_values(by=["vehicle_id", "lat", "timestamp"])

# drop all except first duplicate for each vehicle
print(f"Dropping {duplicate_count} duplicates")
buspositions = buspositions[~buspositions.sort_values(by="timestamp").duplicated(subset=["vehicle_id", "geometry"], keep="first")]


# %% load gtfs (including route shapes)
routes, stops, stop_times, trips, shapes = gtfs.import_gtfs(gtfs_directory)

# trying to determine whether there are different route shapes for each route
trips.groupby(trips.route_id).shape_id.unique()
# counting them up for each route
trips.groupby(trips.route_id).shape_id.unique().apply(len)

# the line frequency function produces a gdf with geometries for every stopping pattern, I think
line_freq = gtfs.lines_freq(stop_times=stop_times, trips=trips, shapes=shapes, routes=routes)
line_freq.crs = {"init": "epsg:4326"}

# map route shapes
if generate_maps==1:
    m = folium.Map(location=[40.733260, -74.164128], tiles="Stamen Toner")
    geojson = line_freq[line_freq.route_id=="1"].to_json()
    #folium.Choropleth( 
    #   geojson,
    #  line_weight=3,
    # line_color="blue").add_to(m)
    #folium.Choropleth(
    #   line_freq[line_freq.route_id=="1"],
    #  line_weight=3,
    # line_color="blue").add_to(m)
    # for adding interactive features: https://vverde.github.io/blob/interactivechoropleth.html
    # compare with code starting around line 1240 here, in gtfs_functions implementation: https://github.com/Bondify/gtfs_functions/blob/master/gtfs_functions/gtfs_funtions.py
    lines = folium.features.GeoJson(
        geojson,
        tooltip = folium.features.GeoJsonTooltip(
            fields = ['route_name']))
    m.add_child(lines)
    # consider folium.features.ColorLine() too

    # save and open html map in browser for viewing
    m.save("map.html")
    webbrowser.open_new("map.html")

# merge trips with shapes and then dissolve on route
shapes_routeids = trips.groupby(by="shape_id")[["route_id"]].max().reset_index()
route_shapes = pd.merge(left=shapes_routeids, right=shapes[["shape_id", "geometry"]], on="shape_id", how="left")  
route_shapes = pd.merge(left=route_shapes, right=routes[["route_id", "route_name"]], on="route_id", how="left")
route_shapes = gpd.GeoDataFrame(route_shapes, geometry="geometry")
lines_dissolved = route_shapes.dissolve(by="route_name").reset_index()
# make new map with pared down route shapes
if generate_maps==1:
    m = folium.Map(location=[40.733260, -74.164128], tiles="Stamen Toner")
    geojson = lines_dissolved.to_json()
    lines = folium.features.GeoJson(
        geojson,
        tooltip = folium.features.GeoJsonTooltip(
            fields = ['route_name']))
    m.add_child(lines)
    m.save("map.html")
    webbrowser.open_new("map.html")

# %% some basic summary

first_observation = buspositions.timestamp.min()
last_observation = buspositions.timestamp.max()
logging.info(f"Dataset spans {first_observation} to {last_observation}")

unique_routes = buspositions.route_number.unique()
logging.info(f"Number of unique routes: {len(unique_routes)}")


# %% calculate speeds
# associate each point with the nearest dissolved gtfs shape of the same route

# calculate distance covered between measurements
#distance_between_measurements = # need to attach to gtfs shapes first
## FIRST calculate length of each route (BUT this is complicated bc not every trip goes the same distance along a route...)
## THEN interpolate from points to points along each line (see here: https://gis.stackexchange.com/questions/306838/snap-points-shapefile-to-line-shapefile-using-shapely)
## (see also here: https://shapely.readthedocs.io/en/stable/manual.html)
## THEN can multiply fraction of total route length by total route length to get cumulative distance to each measured vehicle position
## use that for calculating distance traveled between each measurement

# MAY HAVE TO manually identify corridors where all service patterns are exactly the same
# OR identify routes with very few different service patterns
# could at least start by creating different route shapes for inbound and outbound, because that you at least have from the headsign message
# (though you'd have to know which direction corresponds to inbound and which to outbound in GTFS, which might not be so simple)
# actually it looks like the headsigns correspond with each service pattern, so may be able to associate shapes with particular headsigns
# this would be in the route_shapes gdf as you currently have it
# BUT not sure that the route shape linestrings correspond to different service patterns
# NO you can use the headsigns from GTFS just have to merge them like you did with the routeids
# can compare what's in the GTFS trips file with what's returned in buspositions.headsign
trips_full = pd.read_csv("C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs/trips.txt")
trips_full["service_id"] = trips_full["service_id"].astype(str)
service_headsigns = trips_full.groupby(by="trip_headsign")["service_id"].max().reset_index() # assumes each trip headsign is matched to only one service id
# WAIT WAIT probably you have to merge or join on combination of route id, service id, and direction id
# NO WAIT what are you trying to get
# trying to get geometry associated with headsign, so taht you can match the points you get from NJT with a specific route shape
# in shapes gdf, you alrfeady have shape_id and geometry
# in the trips GTFS file, you have shape_idf and trip_headsign
# so the two of those should be enough??
# number of unique headsigns per route given by
unique_headsigns = trips_full.groupby(by=["route_id"])["trip_headsign"].unique().apply(len).sort_values(ascending=False)
trips_full.groupby(by="route_id")["trip_headsign"].unique().apply(len).sum() # total
# number of unique combinations of service id, route id, and direction id given by
trips_full.groupby(by=["service_id", "route_id", "direction_id"]).size().reset_index().rename(columns={0:"count"})
trips_full.groupby(by=["route_id", "trip_headsign"]).size() # gives number of trips under each headsign for each route
trips_full.groupby(by=["route_id", "service_id", "trip_headsign"]).size()  
# service id doesn't correspond perfectly with headsign. The same headsign will appear under different service ideas for a single route
# DOESN'T HAVE ANYTING TO DO WITH SERVICE ID, which identifies dates when services are available, come on RTFM: https://developers.google.com/transit/gtfs/reference
# REMOVE THESE THINGS AND NOTE IN PUSH 

# calendar dates to sort out how service ids are distributed
pd.crosstab(trips_full.route_id, trips_full.service_id)

#trips = pd.merge(left=trips, right=trips_full[["service_id", "trip_headsign"]], on="service_id", how="left") # add headsign from original GTFS trips file to gdf created by gtfs_functions 
#trips.join(trips_full[["service_id", "trip_headsign"]], on="service_id", how="inner")

# group by date then by vehicle (and maybe then by run?) and sort ascending by time
# by date so that trips that measurements that straddle midnight don't confuse things
# this should give groups that can be iterated over to find time elapsed and distance covered btw measurements
daily_vehicle_groups = buspositions.sort_values(by="timestamp", ascending=True).groupby([buspositions.timestamp.dt.date, "vehicle_id"])

# calculate time elapsed between measurements
# ONLY FOR THOSE POINTS THAT WERE SUCCESSFULLY SNAPPED
#buspositions["time_between_measurements"] = daily_vehicle_groups.timestamp.diff() #.fillna(0)

# probably want to use or consult gtfs_functions' cut_gtfs at some point: https://github.com/Bondify/gtfs_functions/blob/master/gtfs_functions/gtfs_funtions.py

# %% do headsigns returned from web match up with GTFS headsigns?
# if they have headsigns in common, it'll be much easier to make the vehicle positions to GTFS service patterns
trips_full = pd.merge(trips_full.astype({"route_id": "string"}), routes[["route_id", "route_short_name"]], on="route_id")
trips_full[trips_full["route_id"]=="1"]["trip_headsign"].unique()
buspositions[buspositions["route_number"]==1]["headsign"].unique()
# seemingly they all match, except GTFS includes "-Exact Fare" on the ends
# could just split at "-" and have the same names, looks like
for headsign in trips_full[trips_full["route_id"]=="1"]["trip_headsign"].unique():
    print(headsign.rsplit("-")[0])
# ALMOST ALL the same, but some differences,
# for example, there is a pattern that's "1 NWRK 16TH ST" in buspositions, but "1 NEWARK 16TH ST" in GTFS
# SHOULD ALSO STRIP BEGINNING AND ENDING SPACES, bc the spacing around the "-" isn't completely consistent

# number of unique headsigns per route given by
trips_full.groupby(by="route_short_name")["trip_headsign"].unique().apply(len)
trips_full.groupby(by="route_id")["trip_headsign"].unique().apply(len).sum() # total
# number of unique combinations of service id, route id, and direction id given by
trips_full.groupby(by=["service_id", "route_id", "direction_id"]).size().reset_index().rename(columns={0:"count"})
# almost the same
