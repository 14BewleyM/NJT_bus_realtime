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
from shapely.geometry import Point, LineString
import logging
import datetime
import gtfs_functions as gtfs
import webbrowser
import folium
import os
import re

# %% some initial settings
generate_maps = 0

# setting some important directories
main_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper"
code_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/code"

# %% directories and filepaths
gtfs_directory = "C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs"
data_path = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/testdata.csv"

# %% set up logging
logging.basicConfig(filename="log.log", 
                    level=logging.INFO) # set to INFO if you want info messages to log

logging.info(f"Logging started at {datetime.datetime.now()}")

# %% load data and setup dataframe
logging.info(f"Reading in data from {data_path}")
# buspositions = pd.read_csv(data_path)
buspositions = pd.read_csv(data_path)
#TODO@14BewleyM make sure the bigger dataset reads in okay (modify path from testdata.csv to testdata_bigger.csv)
#buspositions = gpd.GeoDataFrame(pd.read_csv(data_path), 
 #                     geometry=gpd.points_from_xy(pd.read_csv(data_path).lon, pd.read_csv(data_path).lat))
buspositions = gpd.GeoDataFrame(buspositions, 
                      geometry=gpd.points_from_xy(buspositions.lon, buspositions.lat))
buspositions.crs = {"init": "epsg:4326"}
print(f"Initial dataset has {buspositions.shape[0]} records")

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
print(f"Final dataset has {buspositions.shape[0]} records") # when run all together, this prints a value different from the actual final row count...

# create time elapsed column
buspositions["month"] = buspositions.timestamp.dt.month
buspositions["day"] = buspositions.timestamp.dt.day
print(f"Creating elapsed time column")
buspositions["time_elapsed"] = buspositions.sort_values(by=["timestamp", "vehicle_id", "run_number"]).groupby(by=["month", "day", "vehicle_id", "run_number"])["timestamp"].diff().fillna(pd.Timedelta(seconds=0))
# testing for negative times
negative_times = (buspositions.time_elapsed < datetime.timedelta(0)).sum()
print(f"{negative_times} records calculated with negative elapsed time")
# create timedelta columns in seconds
buspositions["time_elapsed_seconds"] = buspositions.time_elapsed.apply(pd.Timedelta.total_seconds)
export_bus_positions = 1
if export_bus_positions == 0:
    # drop timedelta column (bc can only be exported to geopackage as seconds, not as timedelta type)
    os.chdir(main_directory)
    buspositions.drop(columns="time_elapsed").to_file("buspositions_cleaned.gpkg", driver="GPKG")
    os.chdir(code_directory)

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

# %% comparing returned headsigns with GTFS headsigns
# if they have headsigns in common, it'll be much easier to make the vehicle positions to GTFS service patterns
trips_full = pd.read_csv("C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs/trips.txt")
trips_full = trips_full.astype(str)
trips_full = pd.merge(trips_full, routes[["route_id", "route_short_name"]], on="route_id")

# number of unique headsigns per route from gtfs given by:
unique_headsigns_gtfs = trips_full.groupby(by=["route_short_name"])["trip_headsign"].unique().apply(len).reset_index().sort_values(by="route_short_name").rename(columns={"trip_headsign": "headsign_count"})
unique_headsigns_gtfs["headsigns"] = trips_full.groupby(by=["route_short_name"])["trip_headsign"].unique().reset_index()["trip_headsign"] # column to hold lists of headsigns
trips_full.groupby(by="route_id")["trip_headsign"].unique().apply(len).sum() # total
# by direction and headsign
trips_full.groupby(by=["direction_id", "trip_headsign"]).size().reset_index().rename(columns={0:"count"})

# headsigns from vehicle measurements
unique_headsigns_buspositions = buspositions.groupby(by="route_number")["headsign"].unique().apply(len).reset_index().sort_values(by="route_number").rename(columns={"headsign": "headsign_count"})
unique_headsigns_buspositions["headsigns"] = buspositions.groupby(by=["route_number"])["headsign"].unique().reset_index()["headsign"] 

# compare routes from gtfs with routes from vehicle measurements
# concerned there may be some route numbers returned in vehicle measurements that don't exist in gtfs as short names
route_diff_buspos_gtfs = set(unique_headsigns_buspositions.route_number.unique()) - set(unique_headsigns_gtfs.route_short_name.unique())
print(f"Following routes are present in bus position dataset but not in GTFS: {route_diff_buspos_gtfs}")
route_diff_gtfs_buspos = set(unique_headsigns_gtfs.route_short_name.unique()) - set(unique_headsigns_buspositions.route_number.unique())
# just drop those two routes, it's not many
print(f"Dropping those {len(route_diff_buspos_gtfs)} routes")
buspositions = buspositions[~buspositions.route_number.isin(route_diff_buspos_gtfs)]

# compare headsign counts for both sets of routes
# drop routes not present in buspositions
unique_headsigns_buspositions = unique_headsigns_buspositions[~unique_headsigns_buspositions.route_number.isin(route_diff_buspos_gtfs)]
# drop routes not present in GTFS, and strip variants of "-Exact fare" as well as trailing and leading spaces
unique_headsigns_gtfs = unique_headsigns_gtfs[unique_headsigns_gtfs.route_short_name.isin(unique_headsigns_buspositions.route_number.unique())]
def split_in_list(patt, repl, list):
    newlist = []
    for element in list: 
        newlist.append(re.sub(patt, repl, element))
        #list[element] = element
    return newlist
unique_headsigns_gtfs["headsigns"] = unique_headsigns_gtfs["headsigns"].apply(lambda x: split_in_list(r"[ ]*-[ ]*[Ee]x.*", "", x))
# new dataframe merging headsigns etc from both datasets
unique_headsigns_buspositions = unique_headsigns_buspositions.rename(columns={"route_number": "route_short_name"})
headsigns_merged = pd.merge(unique_headsigns_buspositions, unique_headsigns_gtfs, on="route_short_name", suffixes=("_buspositions", "_gtfs"))

# compare busposition and gtfs headsigns
# headsigns_merged.headsigns_buspositions.compare(headsigns_merged.headsigns_gtfs)
def intersect_pandas_series(left_series, right_series):
    """
    Takes in two pandas series (each containing lists and of the same length) and returns the intersections of all their elements as a pandas series.
    Elements you want to compare should have the same index values.
    """
    intersection_list = []
    left_list = list(left_series)
    right_list = list(right_series)

    for element in left_series.index:
        intersection_list.append(list(set(left_list[element]).intersection(set(right_list[element]))))
    new_series = pd.Series(intersection_list)
    return new_series

headsigns_merged["headsign_intersection"] = intersect_pandas_series(headsigns_merged.headsigns_buspositions, headsigns_merged.headsigns_gtfs)
headsigns_merged["intersection_count"] = headsigns_merged.headsign_intersection.apply(len)
# find routes where the headsigns match completely (ie lengths of both groups of headsigns and their intersections are all the same)
matching_headsigns = headsigns_merged[(headsigns_merged.headsign_count_buspositions == headsigns_merged.intersection_count) & (headsigns_merged.headsign_count_gtfs == headsigns_merged.intersection_count)]
print(f"Headsigns match for {matching_headsigns.shape[0]} of {headsigns_merged.shape[0]} routes")
matched_route_numbers = matching_headsigns.route_short_name.unique()
# these are the routes the analysis will be run with (before manual matching of headsigns can be done)
# because the matchup allows bus position measurements to be matched with gtfs route shapes

# %% match gtfs route shapes to headsigns and then to bus position measurements
buspositions = buspositions[buspositions.route_number.isin(matched_route_numbers)]
headsign_shapes = trips_full[trips_full.route_short_name.isin(matched_route_numbers)][["route_short_name", "trip_headsign", "direction_id", "shape_id"]]
headsign_shapes = headsign_shapes.drop_duplicates()
# edit headsigns in new dataframes to match buspositions format, using pattern from above
headsign_shapes["trip_headsign"] = headsign_shapes["trip_headsign"].apply(lambda x: re.sub(r"[ ]*-[ ]*[Ee]x.*", "", x))
# change dtypes to more appropriate types
headsign_shapes[["direction_id", "shape_id"]] = headsign_shapes[["direction_id", "shape_id"]].astype(str)
# there are still multiple shapes per headsign
# maybe (hopefully) these are just different parts of what is basically a single shape
# in that case, can just merge them to create one shape per headsign (assuming each headsign corresponds to only one direction)

# comparing shape ids to see if they all match up
len(set(headsign_shapes.shape_id.unique()) - set(shapes.shape_id.unique()))
len(set(headsign_shapes.shape_id.unique()) - set(trips.shape_id.unique()))
len(set(headsign_shapes.shape_id.unique()) - set(trips_full.shape_id.unique()))
# for some reason, there seem to be trip ids in the full trips csv that aren't present in the dfs produced by gtfs_functions

# just construct the shapes yourself, see here for basic procedure: https://www.stevencanplan.com/2016/02/converting-a-transit-agencys-gtfs-to-shapefile-and-geojson-with-qgis/
# should be easy enough with geopandas, group by route and create list from points in sequence order, and then LineString: https://stackoverflow.com/questions/51071365/convert-points-to-lines-geopandas
shapes_full_csv = pd.read_csv("C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs/shapes.txt")
shapes_full_csv.shape_id = shapes_full_csv.shape_id.astype(str)
geometry_gtfs = gpd.points_from_xy(x=shapes_full_csv.shape_pt_lon, y=shapes_full_csv.shape_pt_lat) # create points from latlon
shapes_full = gpd.GeoDataFrame(shapes_full_csv, geometry=geometry_gtfs) # add points to create geodataframe

shapes_full = shapes_full.sort_values(by="shape_pt_sequence").groupby(by=["shape_id"])["geometry"].apply(lambda x: LineString(x.tolist())).reset_index()
# check ids like above
len(set(headsign_shapes.shape_id.unique()) - set(shapes_full_csv.shape_id.unique()))
len(set(headsign_shapes.shape_id.unique()) - set(shapes_full.shape_id.unique()))
# good, complete overlap this time

# merge in geometry
headsign_shapes = pd.merge(headsign_shapes, shapes_full, on="shape_id")
headsign_shapes = gpd.GeoDataFrame(headsign_shapes, geometry="geometry")
headsign_shapes.crs = {"init": "epsg:4326"}
# dissolve on headsign
headsign_shapes_dissolved = headsign_shapes.dissolve(by="trip_headsign").reset_index()

# TODO@14BewleyM Need to make sure there are no headsigns associated with both inbound and outbound direction
# bc you're relying on that for matching the busposition data to the gtfs shapes
# are there any headsigns that have some value other than zero for both the 0 and 1 columns?
headsign_direction_crosstab = pd.crosstab(headsign_shapes.trip_headsign, headsign_shapes.direction_id).reset_index().rename(columns={"0": "direction0", "1": "direction1"})
number_headsigns_two_directions = headsign_direction_crosstab[(headsign_direction_crosstab.direction0!=0) & (headsign_direction_crosstab.direction1!=0)].shape[0]
print(f"There are {number_headsigns_two_directions} headsigns associated with more than one direction")
print(f"Headsigns with both inbound and outbound shapes are: \n {headsign_direction_crosstab[(headsign_direction_crosstab.direction0!=0) & (headsign_direction_crosstab.direction1!=0)].trip_headsign}")
# currently, there are 4 routes with headsigns that are associated with both inbound and outbound directions
# from examining the 29, can at least tell that 0 seems to be inbound and 1 outbound
# look into how geopandas or shapely or movingpandas would calculate the distance btw two points on a line to see if this would cause problems after you merge shapes

# should be safe to dissolve by headsign and direction id. Different routings in the same direction, if they exist, shouldn't pose as much of a risk of miscalculating distances.
headsign_shapes_dissolved_bydirection = headsign_shapes.dissolve(by=["trip_headsign", "direction_id"]).reset_index()
# define some additional dataframes so things don't get way too long
# first, take the headsigns from the crosstab defined above
headsigns_two_directions = headsign_direction_crosstab[(headsign_direction_crosstab.direction0!=0) & (headsign_direction_crosstab.direction1!=0)].trip_headsign.unique()
# then dissolve the shapes associated with those headsigns
headsigns_two_directions_dissolved = headsign_shapes_dissolved_bydirection[headsign_shapes_dissolved_bydirection.trip_headsign.isin(headsigns_two_directions)].dissolve(by="trip_headsign").reset_index()
# then get the indexes of those shapes
headsigns_two_directions_index = list(headsign_shapes_dissolved_bydirection[headsign_shapes_dissolved_bydirection.trip_headsign.isin(headsigns_two_directions)].index)
# use the indexes to delete rows from first headsigns dataframe and replace with dissolved shapes (ie rows where the two-direction headsigns are replaced with single-direction headsigns)
headsign_shapes_dissolved_bydirection = headsign_shapes_dissolved_bydirection.drop(headsigns_two_directions_index) # drop rows for headsigns that have two directions
headsign_shapes_dissolved_bydirection = pd.concat([headsign_shapes_dissolved_bydirection, headsigns_two_directions_dissolved]) # add newly dissolved rows for those headsigns back into the dataframe
# seems to produce the right number of rows

# export to geopackage for visual inspection
export_maps = 1
if export_maps == 0:
    os.chdir(main_directory)
    headsign_shapes.to_file("headsign_shapes.gpkg", driver="GPKG")
    headsign_shapes_dissolved.to_file("headsign_shapes_dissolved.gpkg", driver="GPKG")
    headsign_shapes_dissolved_bydirection.to_file("headsign_shapes_dissolved_bydirection.gpkg", driver="GPKG")
    os.chdir(code_directory)


# TODO@14BewleyM subset buspositions dataset so that you use only observations from routes in the final headsigns dataframe

# now can add shapeid (and also geometry) to buspositions dataset
# but is that what you need to do to interpolate?

# interpolate busposition points to lines based on headsign and shape ids
# https://stackoverflow.com/questions/33769860/pandas-apply-but-only-for-rows-where-a-condition-is-met
# using interpolate: https://gis.stackexchange.com/questions/306838/snap-points-shapefile-to-line-shapefile-using-shapely
# gotta iterate somehow, or use headsign_shapes dataframe as condition in lambda function idk
buspositions["geometry_interpolated"] = buspositions[["headsign", "geometry"]].apply(lambda row: headsign_shapes_dissolved_bydirection[headsign_shapes_dissolved_bydirection.trip_headsign == row.headsign].geometry.interpolate(headsign_shapes_dissolved_bydirection[headsign_shapes_dissolved_bydirection.trip_headsign == row.headsign].geometry.project(row.geometry)), axis=1)
testrow = buspositions.iloc[15]
testpoint = testrow.geometry
testline = headsign_shapes_dissolved_bydirection[headsign_shapes_dissolved_bydirection.trip_headsign==testrow.headsign].geometry
testline.interpolate(testline.project(testpoint))
# all of that works, so it should work when doing apply
# takes forever, though, so how to do this
# try to do for loop over the headsigns
testgdf = gpd.GeoDataFrame(columns=["geometry"])
for headsign in headsign_shapes_dissolved_bydirection.trip_headsign.unique():
    line = headsign_shapes_dissolved_bydirection[headsign_shapes_dissolved_bydirection.trip_headsign == headsign].geometry
    buspositions[buspositions.headsign == headsign].apply(lambda row: line.interpolate(line.project(row.geometry)), axis=1)
    # that works and seemingly much quicker
    # but how to add to dataframe? 
    # try to concat within testgdf to see what happens
# could also use buffers, might be quicker, but how big should buffers be?
# TODO@14BewleyM try to use buffers to 

# this way, what happens when there is no matching headsign in the headsigns dataframe? 

#%% calculate distance covered between measurements
# may want to use or consult gtfs_functions' cut_gtfs at some point: https://github.com/Bondify/gtfs_functions/blob/master/gtfs_functions/gtfs_funtions.py
## FIRST calculate length of each route (BUT this is complicated bc not every trip goes the same distance along a route...)
## THEN interpolate from points to points along each line (see here: https://gis.stackexchange.com/questions/306838/snap-points-shapefile-to-line-shapefile-using-shapely)
## (see also here: https://shapely.readthedocs.io/en/stable/manual.html)
## THEN can multiply fraction of total route length by total route length to get cumulative distance to each measured vehicle position
## use that for calculating distance traveled between each measurement

## some more links about methods for calculating distance (both from start of shape and between points)
### from start of shape: https://gis.stackexchange.com/questions/280112/points-layer-distance-from-the-start-of-line-layer-in-qgis
### between points: https://gis.stackexchange.com/questions/184554/measure-distance-between-points-along-a-line-in-qgis
### if possible, probably best to measure between points, because merging the route shapes by headsign and direction may still produce strange routings that don't represent routes real-world vehicles take



#distance_between_measurements = # need to attach to gtfs shapes first

# %% calcuate speeds


# %% some basic summary

first_observation = buspositions.timestamp.min()
last_observation = buspositions.timestamp.max()
logging.info(f"Dataset spans {first_observation} to {last_observation}")

unique_routes = buspositions.route_number.unique()
logging.info(f"Number of unique routes: {len(unique_routes)}")

