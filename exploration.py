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
# incorporating multiprocessing: https://stackoverflow.com/questions/41240067/pandas-and-multiprocessing-memory-management-splitting-a-dataframe-into-multipl

# %% to do

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
import statsmodels.stats.api as sms
import statsmodels.stats.descriptivestats as smsd
import statsmodels.stats.weightstats as smsw
from shapely.geometry import Point, LineString
import logging
import datetime
import gtfs_functions as gtfs
import webbrowser
import fiona
import folium
import os
import re
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, DateTime, Boolean
from geoalchemy2 import Geometry
import qgis.core as qcore
import qgis.processing as qprocessing
import plotly.graph_objs as go
import seaborn as sns
#import chart_studio.plotly as py

# setting some important directories and filepaths
main_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/"
code_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/code/"
gtfs_directory = "C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs/"
data_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/data/"
data_file = "data_deduped.csv"
data_file_bigger = "testdata_bigger.csv"
data_file_smaller = "testdata.csv"
crosswalk_path = data_directory + "crosswalk.xlsx"
data_path_tract = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/nhgis0008_csv/nhgis0008_ds244_20195_2019_tract.csv"
data_path_blockgroup = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/nhgis0008_csv/nhgis0008_ds244_20195_2019_blck_grp.csv"
boundary_path_tract = "C:/Users/Bewle/OneDrive/Documents/data/geographic/boundaries/2019_NJ_tracts/US_tract_2019.shp"
boundary_path_blockgroup = "C:/Users/Bewle/OneDrive/Documents/data/geographic/boundaries/2019_NJ_block_groups/NJ_blck_grp_2019.shp"

# %% initial settings
# projections
NJ_equidistant = "EPSG:3424"
equal_area = "EPSG:5070"
original_projection = "EPSG:4326"

projections_dict = {"equidistant": NJ_equidistant,
                    "equal_area": equal_area,
                    "original_projection": original_projection}

# %% load buspositions data
def load_buspositions_data(projections_dict, data_file=None, data_directory=None): 
    # provide data_path if want to read in data from cvs, data_directory if want to export from postgresql to csv
    # provide projections as a dict of projections containing an original projection, an equal area projection, and an equidistant projection
    # provide just projections if want to load with postgresql but not export to csv
    if projections_dict is None:
        raise TypeError("Must provide a dict of projections containing an original projection, an equal area projection, and an equidistant projection")
    if data_file is None: # read in data from postgresql database
        print(f"Reading in data from postgresql database")

        with open("database_pass.txt", "r") as file:
            lines = file.readlines(0)
            database_pass = lines[0]
        
        engine = sqlalchemy.create_engine(f"postgresql://postgres:{database_pass}@localhost:5432/bustest", echo=True)
        
        # identify busposition table, as per: https://stackoverflow.com/questions/30785892/simple-select-statement-on-existing-table-with-sqlalchemy
        metadata = sqlalchemy.MetaData(bind=None)
        buspositions = sqlalchemy.Table("buspositions",
                                        metadata,
                                        autoload=True,
                                        autoload_with=engine
                                        )
        # sqlalchemy tutorial on querying: https://docs.sqlalchemy.org/en/14/orm/tutorial.html#querying
        # a more helpful sqlalchemy tutorial: https://hackersandslackers.com/database-queries-sqlalchemy-orm/

        # select distinct rows (rows that probably aren't duplicated observations)
        # based on this method: https://stackoverflow.com/questions/54418/how-do-i-or-can-i-select-distinct-on-multiple-columns
        sql = "SELECT * FROM buspositions b \
                                        WHERE NOT EXISTS( \
                                            SELECT FROM buspositions b1 \
                                            WHERE b.run_number = b1.run_number \
                                            AND b.vehicle_id = b1.vehicle_id \
                                            AND b.lon = b1.lon \
                                            AND b.lat = b1.lat \
                                            AND date_trunc('day', b.timestamp)::date = date_trunc('day', b1.timestamp)::date \
                                            AND b.headsign = b1.headsign \
                                            AND b.id <> b1.id \
                                            AND b.has_service = True \
                                            )"
        #distinct_rows = engine.execute(sqlalchemy.text(sql))
        #print(f"Result after deduping with SQL has {distinct_rows.rowcount} rows")
        print(f"Reading in and deduplicating dataset from database")
        buspositions = pd.read_sql_query(sql, engine) # should maybe try to wrap in sqlalchemy.text(), as best practice
        print(f"Dataset has {buspositions.shape[0]} rows")
        print("Creating date column")
        buspositions["date"] = buspositions.timestamp.dt.date
        print(f"Removing {(buspositions.has_service == False).sum()} rows with no service")
        buspositions = buspositions[buspositions.has_service==True]

        # some headsigns had ampersand replaced with &amp;amp;
        # put ampersand back
        print(f"Replacing ampersands in {buspositions.headsign.str.contains('amp;amp;').sum()} rows that had ampersands replaced with amp;amp;")
        buspositions.headsign = buspositions.headsign.str.replace("amp;amp;", "")

        if data_directory is not None:
            print("Exporting buspositions file to csv")
            buspositions.to_csv(data_directory + "/buspositions_deduped.csv")

        print(f"Converting dataframe to geodataframe and converting to equidistant projection")
        buspositions = gpd.GeoDataFrame(buspositions, 
                            geometry=gpd.points_from_xy(buspositions.lon, buspositions.lat))
        buspositions.crs = projections_dict["original_projection"]
        buspositions = buspositions.to_crs(projections_dict["equidistant"])                
    
    
    elif (data_file is not None) and (data_directory is not None): # read in deduped data from path
        data_path = data_directory + data_file
        print(f"Reading in data from {data_path}")
        filetype = data_path.split(".")[1]
        if filetype == ".csv":
            chunk_size = 1.5 * (10 ** 6) # number of rows to read in at a time

            cols = ["timestamp","vehicle_id", "route_number", "run_number", "headsign", "has_service", "lon", "lat"]
            buspositions = pd.DataFrame(cols)

            with pd.read_csv(data_path, chunksize=chunk_size) as reader:
                iterator = 0
                for chunk in reader:
                    print(f"Chunk {iterator}: Reading in chunk {iterator} with chunk size of {chunk_size} max rows")

                    chunk = gpd.GeoDataFrame(chunk, 
                                        geometry=gpd.points_from_xy(chunk.lon, chunk.lat))
                    chunk.crs = projections_dict["original_projection"]

                    # set timestamp column to be datetime and add date column
                    chunk.timestamp = pd.to_datetime(chunk.timestamp)
                    chunk["date"] = chunk.timestamp.dt.date

                    # set other data types
                    chunk.vehicle_id = chunk.vehicle_id.astype(str).str.split(".").str[0]
                    chunk.route_number = chunk.route_number.astype(str).str.split(".").str[0]
                    chunk.astype({"run_number": "str", "headsign": "str", "lon": "str", "lat": "str"})

                    print(f"Chunk {iterator}: Replacing ampersand in {~chunk.headsign.str.contains('amp;amp;').sum()} rows that have bad formatting")
                    chunk.headsign = chunk.headsign.str.replace("amp;amp;", "")

                    # check for weirdly formatted bool and drop observations for routes that returned no service
                    # drop observations for routes that returned no service
                    if chunk.has_service.dtype != bool:
                        chunk.has_service = chunk.has_service.map({"t":True, "f":False})
                    print(f"Chunk {iterator}: Dropping {chunk.shape[0] - chunk.has_service.sum()} observations from routes with no service")
                    chunk = chunk[chunk.has_service==True]

                    # check for duplicates by date, run number, vehicle id, and geometry
                    # records with the same value for all of these probably reflect observations when a vehicle was moving but its position was not updated yet
                    duplicate_count = chunk.duplicated(subset=["date", "run_number", "headsign", "vehicle_id", "geometry"], keep="first").sum()
                    print(f"Chunk {iterator}: There are {duplicate_count} observations with duplicate days, run numbers, headsigns, vehicle ids, and geometries")

                    # drop all except first duplicate for each vehicle
                    print(f"Chunk {iterator}: Dropping {duplicate_count} duplicates")
                    chunk = chunk.drop_duplicates(subset=["date", "run_number", "vehicle_id", "geometry"], keep="first")
                    buspositions = buspositions.append(chunk)
                    print(f"Chunk {iterator}: Dataset has {buspositions.shape[0]} records after dropping") # when run all together, this prints a value different from the actual final row count...

                    iterator += 1        

        if filetype == ".gpkg":
            print("Loading buspositions from geopackage - METHOD INCOMPLETE, WILL NOT LOAD")
            
    # create time elapsed column
    print(f"Creating elapsed time column")
    buspositions["time_elapsed"] = buspositions.sort_values(by=["timestamp", "vehicle_id", "run_number"]).groupby(by=["date", "vehicle_id", "run_number"])["timestamp"].diff().fillna(pd.Timedelta(seconds=0))
    # reformat headsigns
    print("Reformatting headsigns")
    buspositions.headsign = buspositions.headsign.str.replace("  ", " ")
    buspositions.headsign = buspositions.headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)
    # testing for negative times
    negative_times = (buspositions.time_elapsed < datetime.timedelta(0)).sum()
    print(f"{negative_times} records calculated with negative elapsed time")
    # create timedelta columns in seconds
    buspositions["time_elapsed_seconds"] = buspositions.time_elapsed.apply(pd.Timedelta.total_seconds)

    return buspositions

# %% load gtfs
def create_headsign_crosswalk(buspositions, gtfs_directory, crosswalk_path=None, export_directory=None, export_filename=None):
    # load and headsigns and return a merged dataframe that's the union of headsigns from buspositions and gtfs
    # also add a headsign crosswalk for later use in associating buspositions headsigns with shapes from gtfs (and optionally use a csv file to populate crosswalk field)
    # also return trips dataframe for later use

    # routes, stops, stop_times, trips, shapes = gtfs.import_gtfs(gtfs_directory)
    # load routes as dataframe
    print("Loading routes from gtfs")
    routes = gtfs.import_gtfs(gtfs_directory)[0]

    # load trips as dataframe
    print("Loading trips from gtfs")
    trips = pd.read_csv(gtfs_directory + "trips.txt").astype(str)
    trips = pd.merge(trips, routes[["route_id", "route_short_name"]], on="route_id")

    # unique headsigns
    ## from gtfs
    print("Creating dataframe of unique headsigns from gtfs")
    headsigns_gtfs = trips[["trip_headsign", "route_short_name"]].drop_duplicates(subset=["trip_headsign", "route_short_name"]).rename(columns={"trip_headsign": "headsign"})
    headsigns_gtfs.headsign = headsigns_gtfs.headsign.str.replace("  ", " ")
    headsigns_gtfs.headsign = headsigns_gtfs.headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)

    ## from buspositions
    print("Creating dataframe of unique headsigns from buspositions")
    headsigns_buspositions = buspositions[["headsign", "route_number"]].drop_duplicates(subset=["headsign", "route_number"]).rename(columns={"route_number": "route_short_name"})
    headsigns_buspositions.headsign = headsigns_buspositions.headsign.str.replace("  ", " ")
    headsigns_buspositions.headsign = headsigns_buspositions.headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)

    ## merged
    print("Merging gtfs and buspositions headsigns")
    #headsigns_merged = headsigns_buspositions.merge(headsigns_gtfs, how="outer", on="headsign", suffixes=("_buspositions", "_gtfs")).sort_values(by="headsign")
    headsigns_merged = headsigns_gtfs.merge(headsigns_buspositions, how="outer", on="headsign", suffixes=("_gtfs", "_buspositions")).sort_values(by="headsign")
    
    ## create crosswalk field depending on crosswalk argument
    if crosswalk_path is not None:
        print("Reading in crosswalk values from provided crosswalk file")
        filetype = crosswalk_path.split(".")[1]
        if filetype == "xlsx":
            crosswalk_df = pd.read_excel(crosswalk_path, usecols=["headsign", "headsign_crosswalk"])
        elif filetype == "csv":
            crosswalk_df = pd.read_excel(crosswalk_path, usecols=["headsign", "headsign_crosswalk"])
        else: 
            raise TypeError("Must pass either a csv or xlsx file to use for crosswalk")
        headsigns_merged = headsigns_merged.merge(crosswalk_df, how="left", on="headsign")
    elif (export_directory is not None) and (export_filename is not None):
        ### pre-populate crosswalk field if the headsigns match
        print(f"Exporting headsign crosswalk to {export_directory + export_filename}")
        headsigns_merged["headsign_crosswalk"] = headsigns_merged[(headsigns_merged.route_short_name_buspositions.notna()) & (headsigns_merged.route_short_name_gtfs.notna())].headsign
        headsigns_merged.to_csv(export_directory + export_filename)
    else:
        print("No crosswalk path or export information passed. Returning merged headsigns and gtfs trips as dataframes with new crosswalk field.")
        headsigns_merged["headsign_crosswalk"] = headsigns_merged[(headsigns_merged.route_short_name_buspositions.notna()) & (headsigns_merged.route_short_name_gtfs.notna())].headsign

    return headsigns_merged, trips

# %% compare headsigns
def get_matched_routes(headsigns_merged, buspositions):
    # rows where headsign_crosswalk is not null are rows where buspositions headsigns can be matched to gtfs headsigns
    print(f"There are {headsigns_merged.headsign_crosswalk.notna().sum()} matched headsigns")
    # routes that appear in a list of routes with null values for headsign_crosswalk are routes that aren't entirely matched
    # need separate checks for gtfs and buspositions, bc the route names aren't exactly the same
    unmatched_routes_buspositions = headsigns_merged[headsigns_merged.headsign_crosswalk.isna()].route_short_name_buspositions.unique()
    matched_routes_buspositions = headsigns_merged[~headsigns_merged.route_short_name_buspositions.isin(unmatched_routes_buspositions)].route_short_name_buspositions.unique()
    unmatched_routes_gtfs = headsigns_merged[headsigns_merged.headsign_crosswalk.isna()].route_short_name_gtfs.unique()
    matched_routes_gtfs = headsigns_merged[~headsigns_merged.route_short_name_gtfs.isin(unmatched_routes_gtfs)].route_short_name_gtfs.unique()
    print(f"There are {len(matched_routes_buspositions)} routes from buspositions whose headsigns have a matching gtfs headsign shape of {len(headsigns_merged.route_short_name_buspositions.unique())} total routes")

    buspositions_matched = buspositions[buspositions.route_number.isin(matched_routes_buspositions)] # select matching observations using matching headsigns from buspositions column
    print(f"Returning {buspositions_matched.shape[0]} observations from matched routes of {buspositions.shape[0]} total observations")
    return buspositions_matched, matched_routes_buspositions, matched_routes_gtfs

# %% construct final buspositions dataset
def create_shapes(matched_routes_gtfs, headsigns_merged, trips, projections_dict):
    print("Creating dataframe of headsigns associated with shape id")
    # get unique shape ids for each headsign
    headsign_shapeids = trips[trips.route_short_name.isin(matched_routes_gtfs)][["trip_headsign","shape_id"]].drop_duplicates()
    headsign_shapeids.trip_headsign = headsign_shapeids.trip_headsign.str.replace("  ", " ")
    headsign_shapeids.trip_headsign = headsign_shapeids.trip_headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)
    headsign_shapeids.shape_id = headsign_shapeids.shape_id.astype(str)

    # merge shape ids into merged headsigns dataframe
    headsign_shapes = headsigns_merged.merge(headsign_shapeids, left_on="headsign_crosswalk", right_on="trip_headsign", how="right")

    # construct shapes from gtfs
    print("Constructing shapes for each headsign from gtfs")
    shapes = pd.read_csv(gtfs_directory + "shapes.txt")
    shapes.shape_id = shapes.shape_id.astype(str)
    geometry_gtfs = gpd.points_from_xy(x=shapes.shape_pt_lon, y=shapes.shape_pt_lat) # create points from latlon
    shapes = gpd.GeoDataFrame(shapes, geometry=geometry_gtfs) # add points to create geodataframe

    # string points together into lines by shape_id
    shapes = shapes.sort_values(by="shape_pt_sequence").groupby(by=["shape_id"])["geometry"].apply(lambda x: LineString(x.tolist())).reset_index()

    # dissolving shape geometries by headsign and merging shapes into headsign dataframe
    # some headsigns are associated with multiple directions, but on inspection these can be safely merged
    print("Dissolving shape geometries by headsign and merging shapes into headsign dataframe")
    headsign_shapes = headsign_shapes.merge(shapes, on="shape_id", how="left")
    headsign_shapes = gpd.GeoDataFrame(headsign_shapes, geometry="geometry")
    headsign_shapes.crs = projections_dict["original_projection"]
    print(f"Setting projection to be equidistant {projections_dict['equidistant']}")
    headsign_shapes = headsign_shapes.to_crs(projections_dict["equidistant"])
    headsign_shapes = headsign_shapes.rename(columns={"trip_headsign": "headsign_gtfs", "headsign": "headsign_buspositions"})
    headsign_shapes = headsign_shapes.dissolve(by="headsign_buspositions").reset_index()
    # the shapes that result after .dissolve() may not be the best representation of the routes, but for calculating incremental distances over short distances they should be okay

    return headsign_shapes    

    # # string points together into lines by shape_id
    # shapes = shapes.sort_values(by="shape_pt_sequence").groupby(by=["shape_id"])["geometry"].apply(lambda x: LineString(x.tolist())).reset_index()    

    # headsign_shapes = trips[trips.route_short_name.isin(matched_routes_gtfs)][["route_short_name", "trip_headsign", "direction_id", "shape_id"]]
    # headsign_shapes = headsign_shapes.drop_duplicates()
    # # edit headsigns in new dataframes to match buspositions format, using pattern from above
    # headsign_shapes["trip_headsign"] = headsign_shapes["trip_headsign"].apply(lambda x: re.sub(r"[ ]*-[ ]*[Ee]x.*", "", x))
    # headsign_shapes.trip_headsign = headsign_shapes.trip_headsign.str.replace("  ", " ")
    # headsign_shapes.trip_headsign = headsign_shapes.trip_headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)
    # # change dtypes to more appropriate types
    # headsign_shapes[["direction_id", "shape_id"]] = headsign_shapes[["direction_id", "shape_id"]].astype(str)    

    # # construct shapes from gtfs
    # print("Constructing shapes for each headsign from gtfs")
    # shapes = pd.read_csv(gtfs_directory + "shapes.txt")
    # shapes.shape_id = shapes.shape_id.astype(str)
    # geometry_gtfs = gpd.points_from_xy(x=shapes.shape_pt_lon, y=shapes.shape_pt_lat) # create points from latlon
    # shapes = gpd.GeoDataFrame(shapes, geometry=geometry_gtfs) # add points to create geodataframe

    # # string points together into lines by shape_id
    # shapes = shapes.sort_values(by="shape_pt_sequence").groupby(by=["shape_id"])["geometry"].apply(lambda x: LineString(x.tolist())).reset_index()

    # # dissolving shape geometries by headsign and merging shapes into headsign dataframe
    # # some headsigns are associated with multiple directions, but on inspection these can be safely merged
    # print("Dissolving shape geometries by headsign and merging shapes into headsign dataframe")
    # headsign_shapes = headsign_shapes.merge(shapes, on="shape_id")
    # headsign_shapes = gpd.GeoDataFrame(headsign_shapes, geometry="geometry")
    # headsign_shapes.crs = projections_dict["original_projection"]
    # print(f"Setting projection to be equidistant {projections_dict['equidistant']}")
    # headsign_shapes = headsign_shapes.to_crs(projections_dict["equidistant"])
    # headsign_shapes = headsign_shapes.dissolve(by="trip_headsign").reset_index()
    # headsign_shapes = headsigns_merged.merge(headsign_shapes, left_on="headsign_crosswalk", right_on="trip_headsign", how="left").drop(columns="headsign")
    # # headsigns_merged.merge(headsign_shapes, left_on="headsign_crosswalk", right_on="trip_headsign", how="left").drop(columns="headsign").route_short_name.unique()

    # return headsign_shapes
    
# %% interpolate points to lines and calculate cumulative distances, and calculate incremental distance and speed
def interpolate_and_calc(buspositions, headsign_shapes, projections_dict):
    print("Matching bus position observations to route shapes")

    # merge in crosswalk to headsign_shapes dataframe in create_shapes function, then do the merge below left on headsign and right on headsign_crosswalk

    #buspositions = buspositions.merge(headsign_shapes[["geometry", "trip_headsign", "shape_id", "direction_id"]].rename(columns={"trip_headsign": "headsign"}), on="headsign", how="left")
    buspositions = buspositions.merge(headsign_shapes[["geometry", "headsign_crosswalk", "headsign_buspositions"]], left_on="headsign", right_on="headsign_buspositions", how="left")    
    buspositions = buspositions.rename(columns={"geometry_x": "point_original", "geometry_y": "route_shape"})
    # set CRS again
    buspositions = gpd.GeoDataFrame(buspositions, crs=projections_dict["equidistant"], geometry="point_original")

    # drop null values (these are likely observations with unmatched headsigns - maybe due to observing new headsigns not present in crosswalk)
    remaining_nulls = buspositions.isna()
    remaining_nulls_headsign = buspositions.headsign_crosswalk.isna()
    print(f"{remaining_nulls.sum().max()} rows are null for at least one field, and {remaining_nulls_headsign.sum()} rows could not be matched with a crosswalk headsign and shape")
    print("Dropping those rows")
    buspositions = buspositions[~remaining_nulls_headsign]

    # find interpolated points on route shapes for each busposition observation
    print("Interpolating bus position observations to route shapes")
    buspositions["point_interpolated"] = buspositions.apply(lambda row: row.route_shape.interpolate(row.route_shape.project(row.point_original)), axis=1)
    # TODO may want to drop observations where the distance between the interpolated point and the original point is above some value
    # this is effectively a buffer, which may be a faster way of selecting the same or nearly the same set of points

    # calculate cumulative distance along each route shape
    print("Calculating cumulative distance for each bus position observation")
    buspositions["distance_traveled_cumul"] = buspositions.apply(lambda row: row.route_shape.project(row.point_interpolated, normalized=False), axis=1)
    # normalized=True produces distances normalized as a portion of the total shape length

    # TODO it is probably much faster to project the points as a numpy array or list
    # see p 38 of the documentation here: https://buildmedia.readthedocs.org/media/pdf/shapely/latest/shapely.pdf
    # but the much longer operation is interpolation, so that may not speed things up much

    # drop route shape and reproject
    # drop point_original column (this can cause problems when exporting to geopackage and can be reconstructed anyway)
    print("Dropping route shape, timedelta, and point_original columns and reproject")
    buspositions = gpd.GeoDataFrame(buspositions.drop(columns=["route_shape", "point_original", "time_elapsed"]), geometry="point_interpolated", crs=projections_dict["equidistant"])

    # calculate incremental distance
    print("Calculating distance between each bus position observation")
    ## grouping by date, vehicle id, headsign, and run number is a sensible way of isolating particular vehicles on particular runs
    buses_sorted_grouped = buspositions.sort_values(by="timestamp").groupby(by=["date", "vehicle_id", "headsign", "run_number"])
    ## calculate distance (difference between distance associated with each observation and its previous observation)
    buspositions["distance_traveled_prev"] = buses_sorted_grouped.distance_traveled_cumul.diff()

    # drop date column (can cause problems when exporting to geopackage and can be reconstructed anyway)
    buspositions = buspositions.drop(columns="date")

    # calculate speeds
    print("Calculating speed for each bus position observation")
    buspositions["speed"] = (buspositions.distance_traveled_prev/5280) / (buspositions.time_elapsed_seconds/ (60 *60))

    # return dataframe
    print("Returning new bus positions geodataframe with interpolated points, distances, and speed")
    return buspositions

# %% do some cleaning (eg remove negative distance and speed measures, as well as average speed)
def clean_speeds(buspositions, speed_cutoff=None):    
    # how many observations have negative speeds? by route
    neg_speeds_byroute = buspositions[buspositions.speed < 0].groupby(by="route_number").size().reset_index().rename(columns={0: "neg_speeds_count"}) 
    neg_speeds_byroute["total_obs_count"] = buspositions.groupby(by="route_number").size().reset_index()[0]
    neg_speeds_byroute["neg_speed_pct"] = neg_speeds_byroute.neg_speeds_count / neg_speeds_byroute.total_obs_count
    neg_speeds_byroute.sort_values(by="neg_speed_pct", ascending=False)
    print(f"There are {(buspositions.speed < 0).sum()} observations with negative speeds, out of {buspositions.shape[0]} total observations")
    print(f"There are {buspositions.speed.isna().sum()} observations with null speeds, out of {buspositions.shape[0]} total observations")
    print(f"These represent {((((buspositions.speed < 0).sum() + buspositions.speed.isna().sum())/buspositions.shape[0]) * 100):.{3}}% of all observations")

    # convert negative speeds to positive, and exclude all speeds greater than provided cutoff, if provided
    print("Converting all negative speeds to positive (they are due to negative calculated incremental distance values)")
    #buspositions.speed = buspositions.loc[buspositions["speed"] < 0, "speed"] * -1
    buspositions.loc[buspositions["speed"] < 0, "speed"] = buspositions.speed * -1
    print(f"There are {(buspositions.speed < 0).sum()} remaining negative speed values")

    if speed_cutoff is not None:
        if type(speed_cutoff) is not int:
            raise TypeError("Must provide speed cutoff as an int (in mph)")
        fast_speeds = buspositions.speed > speed_cutoff
        print(f"Removing {fast_speeds.sum()} observations with speeds greater than provided cutoff (mph)")
        buspositions = buspositions[~fast_speeds]

    return buspositions

# %% load Census data (with some helper functions)
# function to create columns for values of interest
def create_equity_measures(df, by=None, aggfunc=None):
    """
    Takes in a dataframe, and optionally fields to group by and an aggregation function.
    Calculates the equity measures.
    Returns pandas series with the equity measures (counts and proportions). 
    """
    # relevant codes from ACS 2015-2019:
    # GISJOIN: gis join match code
    # TRACTA: tract code
    # BLKGRPA: block group code
    # ALUBE001: total population
    # ALUCE001: total pop for race table
    # ALUCE003: (race) black alone
    # ALUCE005: (race) Asian alone
    # ALUKE001: total pop for ethnicity table
    # ALUKE003: (ethnicity) non-Hispanic white population <- use for constructing "POC" category by subtracting from or dividing by pop
    # ALUKE012: (ethnicity) Hispanic or Latino (all races)
    # ALWYE001: total households for which poverty status determined (denominator for poverty % calculation)
    # ALWYE002: households with incomes below poverty level (numerator for poverty % calculation)
    # codes for variables are the same are the same for tracts and block groups
    if by is not None and aggfunc is None:
        raise TypeError("Must provide both groupby field(s) and aggfunc")
    elif df is None:
        raise TypeError("Must provide a dataframe")
    elif by is None and aggfunc is None:
        prop_poc = (df.aluke001 - df.aluke003) / df.aluke001
        prop_white = df.aluke003 / df.aluke001
        prop_latino_allraces = df.aluke012 / df.aluke001
        prop_black = df.aluce003 / df.aluce001
        prop_asian = df.aluce005 / df.aluce001
        prop_poverty = df.alwye002 / df.alwye001
    elif by is not None and aggfunc is not None: # condition in which both groupby field and aggfunc are passed, to return aggregations of groups
        prop_poc = (df.groupby(by=by).aluke001.agg(aggfunc) - df.groupby(by=by).aluke003.agg(aggfunc)) / df.groupby(by=by).aluke001.agg(aggfunc)
        prop_white = df.groupby(by=by).aluke003.agg(aggfunc) / df.groupby(by=by).aluke001.agg(aggfunc)
        prop_latino_allraces = df.groupby(by=by).aluke012.agg(aggfunc) / df.groupby(by=by).aluke001.agg(aggfunc)
        prop_black = df.groupby(by=by).aluce003.agg(aggfunc) / df.groupby(by=by).aluce001.agg(aggfunc)
        prop_asian = df.groupby(by=by).aluce005.agg(aggfunc) / df.groupby(by=by).aluce001.agg(aggfunc)
        prop_poverty = df.groupby(by=by).alwye002.agg(aggfunc) / df.groupby(by=by).alwye001.agg(aggfunc)
    else: # condition in which only aggfunc is passed, to return aggregation of entire dataframe
        prop_poc = (df.aluke001.agg(aggfunc) - df.aluke003.agg(aggfunc)) / df.aluke001.agg(aggfunc)
        prop_white = df.aluke003.agg(aggfunc) / df.aluke001.agg(aggfunc)
        prop_latino_allraces = df.aluke012.agg(aggfunc) / df.aluke001.agg(aggfunc)
        prop_black = df.aluce003.agg(aggfunc) / df.aluce001.agg(aggfunc)
        prop_asian = df.aluce005.agg(aggfunc) / df.aluce001.agg(aggfunc)
        prop_poverty = df.alwye002.agg(aggfunc) / df.alwye001.agg(aggfunc)
    return prop_poc, prop_white, prop_latino_allraces, prop_black, prop_asian, prop_poverty

# function to create census dataframes
# select relevant tracts by route and service area
# determine equity block groups
# determine equity routes
def add_equity_measures(headsign_shapes, projections_dict, census_boundary_path, census_data_path, buffer_size=None):
    # set up Census data (ACS files from NHGIS, using GISJOIN field to merge)
    # provide buffer size in meters
    ## read in data
    print("Reading in Census data from provided paths")
    acs_blockgroups = pd.read_csv(census_data_path)
    old_cols = list(acs_blockgroups.columns)
    new_cols = [new_col.lower() for new_col in old_cols] 
    acs_blockgroups = acs_blockgroups.rename(columns={old_cols[i]: new_cols[i] for i in range(len(old_cols))})
    ## read in boundaries
    blockgroups = gpd.read_file(census_boundary_path)
    old_cols = list(blockgroups.columns)
    new_cols = [new_col.lower() for new_col in old_cols] 
    blockgroups = blockgroups.rename(columns={old_cols[i]: new_cols[i] for i in range(len(old_cols))})
    # prefer them lowercase

    # create merged dataframe associating block group shapes with their 2015-2019 ACS estimates
    print("Merging Census shapes with data")
    blockgroups_merged = blockgroups.merge(acs_blockgroups, on="gisjoin")

    # dissolve by route
    # TODO if want to do separate analysis by headsign, can avoid doing the dissolve here and just drop duplicate block groups by route name and gisjoin to get a route-level dataset
    print("Dissolving provided headsign shapes file by route")
    route_shapes = headsign_shapes.dissolve(by="route_short_name_buspositions").reset_index()

    # reproject headsign shapes dataframe to equal area
    route_shapes = route_shapes.to_crs(projections_dict["equal_area"])

    if buffer_size is not None:
        print("Creating buffered headsign shapes")
        # create headsign shapes with bufffer 
        buffer = route_shapes.geometry.buffer(buffer_size)
        route_shapes["buffer"] = buffer

        # join block groups that intersect the headsign shape buffers
        blockgroups_routes_joined_buffer = gpd.sjoin(route_shapes.set_geometry("buffer"), blockgroups_merged.to_crs(projections_dict["equal_area"]), op="intersects", how="left")    
        # .to_crs(projections_dict["equal_area"]
    print("Doing spaital join of non-buffered headsign shapes with Census data")
    headsign_columns = ["headsign_buspositions", "headsign_crosswalk", "headsign_gtfs"] # drop headsign columns in dissolved route dataframe so as not to confuse
    routes_blockgroups_joined = gpd.sjoin(route_shapes, blockgroups_merged.to_crs(projections_dict["equal_area"]), op="intersects", how="left").drop(columns=headsign_columns)

    # calculate equity measures
    # create new columns for values of interest
    print("Adding demographic proportions for relevant groups")
    for df in [routes_blockgroups_joined, blockgroups_merged]: # acs_blockgroups]:
        df["prop_poc"], df["prop_white"], df["prop_latino_allraces"], df["prop_black"], df["prop_asian"], df["prop_poverty"] = create_equity_measures(df)    

    # define new dataframes for each route and for the whole service area, to avoid double-counting
    # for each route (dropping block groups that are associated more than once with a route)
    #blockgroups_routes_joined_buffer = blockgroups_headsigns_joined.drop_duplicates(["route_short_name_buspositions", "gisjoin"])
    # for the whole service area (dropping all duplicated block groups, using gisjoin field as unique id for block groups)
    routes_servicearea_joined = routes_blockgroups_joined.drop_duplicates(["gisjoin"])

    # calculate proportions by service area and route
    # service area averages for each of the measures
    # index is each of the values, access by eg service_area_distribution.prop_poc.loc["mean"]
    service_area_distribution = routes_servicearea_joined[["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]].describe()
    # add statewide mean for each of the variables
    service_area_distribution.loc["state_mean"] = blockgroups_merged[["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]].mean() # acs_blockgroups[["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]].mean()

    # calculating measures for the whole service area and for each route
    route_index = list(routes_blockgroups_joined.route_short_name_buspositions.unique())
    route_index.append("service_area")
    columns = ["count_total", "count_poc", "prop_poc", "count_white", "prop_white", "count_latino_allraces", "prop_latino_allraces", "count_black", "prop_black", "count_asian", "prop_asian", "count_poverty_hholds", "prop_poverty"]
    equity_measures_byroute = pd.DataFrame(index=route_index, columns=columns)
    # calculate proportion measures with function (do this also for the merged blockgroups dataframe, for later use in calculating Title VI measures)
    print("Calculating demographic proportions for relevant groups by route")
    equity_measures_byroute["prop_poc"], equity_measures_byroute["prop_white"], equity_measures_byroute["prop_latino_allraces"], equity_measures_byroute["prop_black"], equity_measures_byroute["prop_asian"], equity_measures_byroute["prop_poverty"] = create_equity_measures(routes_blockgroups_joined, by="route_short_name_buspositions", aggfunc=np.sum)
    blockgroups_merged["prop_poc"], blockgroups_merged["prop_white"], blockgroups_merged["prop_latino_allraces"], blockgroups_merged["prop_black"], blockgroups_merged["prop_asian"], blockgroups_merged["prop_poverty"] = create_equity_measures(blockgroups_merged)
    # calculate proportion measures for the entire service area and statewide, and add to dataframe
    print("Calculating figures for entire service area")
    equity_measures_byroute.loc["service_area", ["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]] = list(create_equity_measures(routes_servicearea_joined, aggfunc=np.sum))
    equity_measures_byroute.loc["statewide", ["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]] = list(create_equity_measures(acs_blockgroups, aggfunc=np.sum))
    # calculate total measures by route, as well as for the entire service area and statewide, and add to dataframe
    print("Calculate demographic totals for relevant groups by route")
    totals_map_dict = {"count_total": "alube001",
                        "count_white": "aluke003",
                        "count_latino_allraces": "aluke012",
                        "count_black": "aluce003", 
                        "count_asian": "aluce005",
                        "count_poverty_hholds": "alwye001"}
    for column in totals_map_dict:
        #print(column)
        #print(original_variable)
        equity_measures_byroute[column] = routes_blockgroups_joined.groupby(by="route_short_name_buspositions")[totals_map_dict[column]].sum()
        equity_measures_byroute.loc["service_area", column] = routes_servicearea_joined[totals_map_dict[column]].sum()
        equity_measures_byroute.loc["statewide", column] = acs_blockgroups[totals_map_dict[column]].sum()
    # calculate POC count separately bc you have to subtract
    equity_measures_byroute["count_poc"] = routes_blockgroups_joined.groupby(by="route_short_name_buspositions").aluke001.sum() - routes_blockgroups_joined.groupby(by="route_short_name_buspositions").aluke003.sum()
    equity_measures_byroute.loc["service_area", "count_poc"] = routes_servicearea_joined.aluke001.sum() - routes_servicearea_joined.aluke003.sum()
    equity_measures_byroute.loc["statewide", "count_poc"] = acs_blockgroups.aluke001.sum() - acs_blockgroups.aluke003.sum()

    return equity_measures_byroute, routes_blockgroups_joined, route_shapes.drop(columns=headsign_columns), blockgroups_merged, service_area_distribution #, blockgroups_merged

# function to assign routes as equity routes
def determine_equity_routes(equity_measures_byroute, route_shapes, blockgroups, projections_dict, columns=None, proportion_cutoffs=None, mileage_cutoff=1/3):
    # determine which routes are equity routes 

    # equity_measures_byroute is a dataframe with overall demographic proportions for routes and the service area
    # joined_routes_census_gdf is a dataframe that associates routes with the census geographies they intersect
    # pass a column or list of columns to determine equity route status for those measures (one of prop_poc, prop_white, prop_latino_allraces, prop_black, prop_asian, or prop_poverty)
    # can pass proportion cutoffs (as a decimal), if not will use value for the whole service area for each measure

    # calculate FTA equity measures
    print("Calculating FTA equity measures")

    # reproject geometry field to equidistant projection
    projection = projections_dict["equidistant"]
    for gdf in [route_shapes, blockgroups]:
        if gdf.crs != projection:
            print(f"Reprojecting to equidistant projection {projection}")
            gdf = gdf.to_crs(projection)
            print(f"Projection is now {gdf.crs}")
    # drop duplicate geometries before reprojecting to save time
    # there should be just one geometry per route
    #print("Dropping duplicate geometries before reprojecting and re-merging - be sure there is only one geometry per route in the joined_census_route_df!")
    #route_shapes = route_shapes.merge(route_shapes.drop_duplicates(subset=["route_short_name_buspositions"]).to_crs(projection), on="route_short_name_buspositions", how="left", suffixes=("", "_reprojected"))
    #route_shapes = route_shapes.set_geometry("geometry_reprojected")

    # add equity measures to blockgroups dataframe
    # to determine associate an equity designation with each blockgroup shape
    equity_columns = ["prop_white", "prop_poc", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]
    if columns is not None:
        print("Using provided column values with service area measureas (FTA standard) to calculate equity mileage")        
        equity_columns = columns
    else:
        print("Using default column values with service area measureas (FTA standard) to calculate equity mileage")
    #final_columns = equity_columns.copy()
    #final_columns.append("gisjoin")
    #blockgroups = blockgroups.merge(route_shapes[final_columns], on="gisjoin", how="left")
    #blockgroups = blockgroups.to_crs(projection)

    # create list of final column names for use in final merge
    final_columns = []

    # assign equity designation to blockgroups for each equity measure
    for column in equity_columns:
        # identify service area proportion for given measure
        service_area_prop = equity_measures_byroute.loc["service_area", column]
        if proportion_cutoffs is not None: # this wouldn't work the way you want it to
            service_area_prop = proportion_cutoffs
        # designate block groups as equity block groups if their proportion for the measure is greater than the service area proportion
        blockgroups["equity_bg_" + column] = False
        blockgroups.loc[blockgroups[column] > service_area_prop, "equity_bg_" + column] = True

    # using overlay:
    #overlay = gpd.overlay(blockgroups, route_shapes, how="intersection")
    #return overlay

    # using clip: https://www.earthdatascience.org/courses/use-data-open-source-python/intro-vector-data-python/vector-data-processing/clip-vector-data-in-python-geopandas-shapely/
    # just doing prop poc for now to see how fast it goes
    #clipped_equity = gpd.clip(route_shapes, blockgroups[blockgroups.equity_bg_prop_poc==True])
    #clipped_total = gpd.clip(route_shapes, blockgroups)
    
    # create length columns
    #clipped_equity["length"] = clipped_equity.length
    #clipped_total["length"] = clipped_total.length
    #equity_miles = clipped_equity.length
    route_shapes["total_length"] = route_shapes.length / 5280
    final_columns.append("total_length")
    for column in equity_columns:
        clipped_equity = gpd.clip(route_shapes, blockgroups[blockgroups["equity_bg_" + column]==True])
        route_shapes["equity_miles_" + column] = clipped_equity.length / 5280
        route_shapes["equity_miles_ratio_" + column] = route_shapes["equity_miles_" + column] / route_shapes["total_length"]
        # determine equity routes
        print(f"Determining equity routes (those with more than {mileage_cutoff} of revenue miles in equity block groups) for {column}")
        route_shapes["equity_route_" + column] = False
        route_shapes.loc[route_shapes["equity_miles_ratio_" + column] > mileage_cutoff, "equity_route_" + column] = True
        # add columns to list of names
        final_columns.append("equity_miles_" + column)
        final_columns.append("equity_miles_ratio_" + column)
        final_columns.append("equity_route_" + column)
    
    # add route number column so can use to merge
    final_columns.append("route_short_name_buspositions")
    equity_measures_byroute = equity_measures_byroute.merge(route_shapes[final_columns], left_index=True, right_on="route_short_name_buspositions", how="outer").set_index("route_short_name_buspositions")

    return route_shapes, equity_measures_byroute


def determine_equity_routes_v1(equity_measures_byroute, joined_routes_census_gdf, blockgroups, projections_dict, columns=None, proportion_cutoffs=None, mileage_cutoff=1/3):
    # determine which routes are equity routes 

    # equity_measures_byroute is a dataframe with overall demographic proportions for routes and the service area
    # joined_routes_census_gdf is a dataframe that associates routes with the census geographies they intersect
    # pass a column or list of columns to determine equity route status for those measures (one of prop_poc, prop_white, prop_latino_allraces, prop_black, prop_asian, or prop_poverty)
    # can pass proportion cutoffs (as a decimal), if not will use value for the whole service area

    # calculate FTA equity measures
    print("Calculating FTA equity measures")

    # reproject geometry field to equidistant projection
    projection = projections_dict["equidistant"]
    print(f"Reprojecting line shapes to equidistant projection {projection}")
    # drop duplicate geometries before reprojecting to save time
    # there should be just one geometry per route
    print("Dropping duplicate geometries before reprojecting and re-merging - be sure there is only one geometry per route in the joined_census_route_df!")
    joined_routes_census_gdf = joined_routes_census_gdf.merge(joined_routes_census_gdf.drop_duplicates(subset=["route_short_name_buspositions"]).to_crs(projection), on="route_short_name_buspositions", how="left", suffixes=("", "_reprojected"))
    joined_routes_census_gdf = joined_routes_census_gdf.set_geometry("geometry_reprojected")
    # merge block group geometry back in for later use in calculating lengths of routes within each block group
    joined_routes_census_gdf = joined_routes_census_gdf.merge(blockgroups[["gisjoin", "geometry"]], on="gisjoin", how="left", suffixes=("", "_blockpoly"))
    joined_routes_census_gdf.geometry_blockpoly = joined_routes_census_gdf.geometry_blockpoly.to_crs(projections_dict["equidistant"])
    # may need to make a pull request and use .overlay that works with shapes other than polygons, see here: https://github.com/geopandas/geopandas/issues/821

    # attempt this with QGIS
    # some resources: 
    # https://gis.stackexchange.com/questions/173303/using-each-feature-of-vector-layer-as-clip-mask-over-separate-point-layer-stori
    # https://gis.stackexchange.com/questions/362979/loading-geodataframe-as-qgis-vector-layer-without-exporting-to-shapefile
    # https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/vector.html
    # https://docs.qgis.org/3.16/en/docs/pyqgis_developer_cookbook/crs.html

    # OR maybe do this by creating two dataframes,
    # one of polygons from equity=True block groups, and one for equity=False
    # and then intersect each route shape with each of those geodataframes

    if (type(columns) is not list) and (columns is not None):
        raise TypeError("Must provide columns as a list (even for single values), or pass nothing to use all equity measures")
    elif columns is None: 
        if proportion_cutoffs is not None:
            raise ValueError("If passing proportion cutoffs, column values must also be provided")
        print("Using default column values with service area measureas (FTA standard) to calculate equity mileage")
        columns = ["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]
        for column in columns:
            # identify service area proportion for given measure
            service_area_prop = equity_measures_byroute.loc["service_area", column]
            # designate block groups as equity block groups if their proportion for the measure is greater than the service area proportion
            joined_routes_census_gdf["equity_block_group"] = False
            joined_routes_census_gdf.loc[joined_routes_census_gdf[column] > service_area_prop, "equity_block_group"] = True
            # calculate total length for route shapes
            joined_routes_census_gdf["length"] = joined_routes_census_gdf.geometry.length
            # calculate length of each route associated with equity and non-equity block groups
            ## ORIGINAL
            rev_miles_total = joined_routes_census_gdf.drop_duplicates(subset="route_short_name_buspositions").groupby(by="route_short_name_buspositions").length.sum() / 5280
            rev_miles_equity = joined_routes_census_gdf[joined_routes_census_gdf.equity_block_group==True].groupby(by="route_short_name_buspositions").length.sum() / 5280
            ## NEW
            #equity_gdf = joined_routes_census_gdf[joined_routes_census_gdf.equity_block_group==True][["route_short_name_buspositions", "gisjoin", "geometry_blockpoly"]]
            ###nonequity_gdf = joined_routes_census_gdf[joined_routes_census_gdf.equity_block_group==False][["route_short_name_buspositions", "gisjoin", "geometry_blockpoly"]]
            #route_shapes = joined_routes_census_gdf.drop_duplicates(subset=["route_short_name_buspositions"]).to_crs(projection)
            #rev_miles_total = joined_routes_census_gdf.drop_duplicates(subset="route_short_name_buspositions").groupby(by="route_short_name_buspositions").length.sum() / 5280            
            #rev_miles_equity = route_shapes.apply(lambda row: gpd.clip(gpd.GeoDataFrame(row, geometry=row.geometry), equity_gdf).geometry.length(), axis=1)
            # add lengths to equity measures dataframe as new columns, rev_miles_total, rev_miles_equity, and rev_miles_equity_ratio
            equity_measures_byroute["rev_miles_total"] = rev_miles_total.astype(float)
            equity_measures_byroute["rev_miles_equity_" + column] = rev_miles_equity.astype(float)
            equity_measures_byroute["rev_miles_equity_ratio_" + column] = rev_miles_equity/rev_miles_total
            equity_measures_byroute["equity_route_" + column] = False
            equity_measures_byroute[equity_measures_byroute["rev_miles_equity_ratio_" + column] > mileage_cutoff]["equity_route_" + column] = True
    else: # use provided columns in this case
        if proportion_cutoffs is not None:
            if type(proportion_cutoffs) is not list:
                raise TypeError("Proportion cutoffs must be passed as list (even for single values)")
            if len(columns) != len(proportion_cutoffs):
                raise ValueError("The number of proportion cutoffs must match the number of columns")
            print("Using provided proportion cutoffs to calculate equity mileage for provided columns")
            column_cutoff_dict = dict(zip(columns, proportion_cutoffs))
            projection = projections_dict["equidistant"]
            for column, proportion_cutoff in column_cutoff_dict.items():
                # identify service area proportion for given measure
                service_area_prop = equity_measures_byroute.loc["service_area", column]
                # designate block groups as equity block groups if their proportion for the measure is greater than the service area proportion
                joined_routes_census_gdf["equity_block_group"] = False
                joined_routes_census_gdf.loc[joined_routes_census_gdf[column] > proportion_cutoff, "equity_block_group"] = True
                # calculate total length for route shapes
                joined_routes_census_gdf["length"] = joined_routes_census_gdf.geometry.length
                # calculate length of each route associated with equity and non-equity block groups
                rev_miles_total = joined_routes_census_gdf.groupby(by="route_short_name_buspositions").length.sum() / 5280
                rev_miles_equity = joined_routes_census_gdf[joined_routes_census_gdf.equity_block_group==True].groupby(by="route_short_name_buspositions").length.sum() / 5280
                # add lengths to equity measures dataframe as new columns, rev_miles_total, rev_miles_equity, and rev_miles_equity_ratio
                equity_measures_byroute["rev_miles_total"] = rev_miles_total
                equity_measures_byroute["rev_miles_equity"] = rev_miles_equity
                equity_measures_byroute["rev_miles_equity_ratio"] = rev_miles_equity/rev_miles_total
                equity_measures_byroute["equity_route"] = False
                equity_measures_byroute[equity_measures_byroute.rev_miles_equity_ratio > mileage_cutoff]["equity_route"] = True                
        else: # use provided column values, with service area measures as the cutoff
            print("Using service area values (FTA standard) to calculate equity mileage for provided columns")
            projection = projections_dict["equidistant"]
            for column in columns:
                # identify service area proportion for given measure
                service_area_prop = equity_measures_byroute.loc["service_area", column]
                # designate block groups as equity block groups if their proportion for the measure is greater than the service area proportion
                joined_routes_census_gdf["equity_block_group"] = False
                joined_routes_census_gdf.loc[joined_routes_census_gdf[column] > service_area_prop, "equity_block_group"] = True
                # calculate total length for route shapes
                joined_routes_census_gdf["length"] = joined_routes_census_gdf.geometry.length
                # calculate length of each route associated with equity and non-equity block groups
                rev_miles_total = joined_routes_census_gdf.groupby(by="route_short_name_buspositions").length.sum() / 5280
                rev_miles_equity = joined_routes_census_gdf[joined_routes_census_gdf.equity_block_group==True].groupby(by="route_short_name_buspositions").length.sum() / 5280
                # add lengths to equity measures dataframe as new columns, rev_miles_total, rev_miles_equity, and rev_miles_equity_ratio
                equity_measures_byroute["rev_miles_total"] = rev_miles_total
                equity_measures_byroute["rev_miles_equity"] = rev_miles_equity
                equity_measures_byroute["rev_miles_equity_ratio"] = rev_miles_equity/rev_miles_total
                equity_measures_byroute["equity_route"] = False
                equity_measures_byroute[equity_measures_byroute.rev_miles_equity_ratio > mileage_cutoff]["equity_route"] = True                       
    return equity_measures_byroute

# speeds by block group (spatial join of speed observations with block groups)

# %% add speed data
def add_speed_data(buspositions, equity_measures_byroute, peak_cutoffs={"am_start": 6, "am_end": 9, "pm_start": 16, "pm_end": 18}):
    # add equity route designation to buspositions
    equity_measure_columns = ["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]    
    equity_columns = [("equity_route_" + column) for column in equity_measure_columns]
    full_columns = equity_columns.append("route_short_name_buspositions")
    buspositions = buspositions.merge(equity_measures_byroute[full_columns], left_on="route_number", right_on="route_short_name_buspositions", suffixes=("", ""), how="left")
    
    # create temp dataframe to load speed summary data
    speed_columns = ["speed_mean", "speed_stdev"]
    equity_measures_temp = pd.DataFrame(columns=])

    equity_measures_byroute["speed_mean"] = buspositions.groupby(by="route_number")["speed"].mean().astype(float)
    equity_measures_byroute.loc["service_area", "speed"] = buspositions.speed.mean().astype(float)
    equity_measures_byroute["speed_stdev"] = buspositions.groupby(by="route_number")["speed"].std().astype(float)
    equity_measures_byroute.loc["service_area", "speed"] = buspositions.speed.std().astype(float)
    equity_measures_byroute["speed_coeff_var"] = equity_measures_byroute.speed_stdev / equity_measures_byroute.speed_mean
    # R package about comparing coefficients of variation and testing for significance: https://cran.r-project.org/web/packages/cvequality/vignettes/how_to_test_CVs.html

    # calculate peak and off-peak speeds
    morning_peak = (buspositions.timestamp.dt.hour > peak_cutoffs["am_start"]) & (buspositions.timestamp.dt.hour < peak_cutoffs["am_end"])
    afternoon_peak = (buspositions.timestamp.dt.hour > peak_cutoffs["pm_start"]) & (buspositions.timestamp.dt.hour < peak_cutoffs["pm_end"])
    equity_measures_byroute["speed_peak_mean"] = buspositions[(morning_peak) | (afternoon_peak)].groupby(by="route_number").speed.mean().astype(float)
    equity_measures_byroute["speed_offpeak_mean"] = buspositions[(~morning_peak) & (~afternoon_peak)].groupby(by="route_number").speed.mean().astype(float)
    equity_measures_byroute["speed_peak_stdev"] = buspositions[(morning_peak) | (afternoon_peak)].groupby(by="route_number").speed.std().astype(float)
    equity_measures_byroute["speed_offpeak_stdev"] = buspositions[(~morning_peak) & (~afternoon_peak)].groupby(by="route_number").speed.std().astype(float)    
    equity_measures_byroute["speed_peak_coeff_var"] = equity_measures_byroute.speed_peak_stdev / equity_measures_byroute.speed_peak_mean
    equity_measures_byroute["speed_offpeak_coeff_var"] =  equity_measures_byroute.speed_offpeak_stdev / equity_measures_byroute.speed_offpeak_mean

    return buspositions, equity_measures_byroute
# %% export geodataframes as layers to single geopackage
# some resources:
# https://gis.stackexchange.com/questions/307749/open-gpkg-embedded-layers-in-python
# https://fiona.readthedocs.io/en/latest/fiona.html#fiona.open
# https://stackoverflow.com/questions/54562069/multi-layer-gdb-files-in-python/54563846#54563846
# https://stackoverflow.com/questions/56165069/can-geopandas-get-a-geopackages-or-other-vector-file-all-layers

# %% analysis

# speed by route

# coefficient of variation

# %% figures
# see the heatmaps and line charts from gtfs_functions readme: https://github.com/Bondify/gtfs_functions

# map of routes that made the cut, maybe over counties or light road network

# then map showing equity measure (pct poc), with service area average

# then map showing which routes are above and below the average
# (i.e., equity and non-equity routes)

# simple comparison of mean speeds equity vs nonequity routes


# build up line plots of average speeds over time (or box plots by hour, which may be more easy to parse)
def line_plots_speed(buspositions):
# define traces for plotting with plotly
    plot = go.Scatter(x=buspositions.timestamp, y=buspositions.speed)
    layout = go.Layout()


    plotting_dict = {}
    for equity_status in [True, False]:
        plotting_dict[equity_status] = buspositions[buspositions["equity_route_prop_poc"]==equity_status]

    i=0
    trace_dict_mean = {}
    trace_dict_upper = {}
    trace_dict_lower = {}
    for route, df in resampled_dict.items():
        color_new = colors_rutgers_byroute[route] # colors_rutgers[i%5]
        color_new_transparent = colors_rutgers_byroute_transparent[route] # colors_rutgers_transparent[i%5]
        i += 1
        
        trace_dict_mean[route] = go.Scatter(
            x=df.index,
            y=df["speed"],
            name="Avg. speed (mph), " + route,
            line_color=color_new,
            opacity=0.8)
        
        # add shaded areas in same color for 95% conf int
        trace_dict_upper[route] = go.Scatter(
            x=df.index,
            y=df["upper_bound_95%"],
            name="95% CI upper bound, " + route,
            #fill="tonexty",
            line=dict(color="rgba(255,255,255,0)"),
            fillcolor=color_new_transparent,
            opacity=0.3,
            showlegend=False)
        
        trace_dict_lower[route] = go.Scatter(
            x=df.index,
            y=df["lower_bound_95%"],
            name="95% conf. int., " + route,
            fill="tonexty",
            line=dict(color="rgba(255,255,255,0)"),
            fillcolor=color_new_transparent,
            opacity=0.3,
            showlegend=True)
    ## first, over time for nonequity routes
    ## then, over time for equity routes
    ## then, overlay



    ## then an example for a couple routes?

# then do the same for variation

# then consider spatial variation

## map showing equity and nonequity routes

# %% main routine
def main():
    read_from_file = 0
    export_dataframes = 0

    if read_from_file==1: # read in data from files
        buspositions = gpd.read_file("buspositions.gpkg", layer="buspositions")
        headsign_shapes = gpd.read_file("buspositions.gpkg", layer="headsign_shapes")

    else: # create dataframes
        # load from postgresql db
        buspositions = load_buspositions_data(projections_dict=projections_dict) 
        # create merged headsigns and trips dataframes (with crosswalk from file)
        headsigns_merged, trips = create_headsign_crosswalk(buspositions, gtfs_directory, crosswalk_path=crosswalk_path) # be on the lookout for weird stuff with GO28/258, not sure you did the crosswalk right
        # create dataframe with only observations from routes that have been entirely matched with crosswalk (ie mappable using gtfs shapes)
        buspositions_matched, matched_routes_buspositions, matched_routes_gtfs = get_matched_routes(headsigns_merged, buspositions)
        # create headsign shapes
        headsign_shapes = create_shapes(matched_routes_gtfs, headsigns_merged, trips, projections_dict=projections_dict)
        # interpolate and calculate cumulative distances for each bus observation
        buspositions = interpolate_and_calc(buspositions_matched, headsign_shapes, projections_dict)
        # remove excessive speeds
        buspositions = clean_speeds(buspositions, speed_cutoff=70)

        # export files as layers to single geopackage
        if export_dataframes==1:
            buspositions.to_file("buspositions.gpkg", driver="GPKG", layer="buspositions")
            headsign_shapes.to_file("buspositions.gpkg", driver="GPKG", layer="headsign_shapes")
    
    # associate equity measures with routes
    equity_measures_byroute, routes_blockgroups_joined, route_shapes, blockgroups_withdata, service_area_distribution = add_equity_measures(headsign_shapes, projections_dict, boundary_path_blockgroup, data_path_blockgroup, buffer_size=None)
    #equity_measures_byroute = determine_equity_routes(equity_measures_byroute, routes_blockgroups_joined, blockgroups, projections_dict)
    #clipped, equity_measures_byroute_v2 = determine_equity_routes_v2(equity_measures_byroute, routes_blockgroups_joined, blockgroups, projections_dict)
    route_shapes, equity_measures_byroute = determine_equity_routes(equity_measures_byroute, route_shapes, blockgroups_withdata, projections_dict)

    # associate speed and variability measures with block groups
    buspositions_new, equity_measures_byroute = add_speed_data(buspositions, equity_measures_byroute)

    print("Done")
if __name__ == "__main__":
    main()

# %% test visualizations
# some seaborn examples: https://seaborn.pydata.org/generated/seaborn.lineplot.html#seaborn.lineplot
sns.lineplot(data=buspositions, x=buspositions.timestamp.dt.hour, y="speed")
sns.lineplot(data=buspositions, # by weekday
            x=buspositions.timestamp.dt.hour, 
            y="speed",
            hue=buspositions.timestamp.dt.weekday)
sns.lineplot(data=buspositions, # hue by weekday vs weekend
            x=buspositions.timestamp.dt.hour, 
            y="speed",
            hue=buspositions.timestamp.dt.weekday.isin([5, 6]))
# includes 95% confidence intervals by default!
routes = pd.read_csv(gtfs_directory + "routes.txt")
routes = routes.merge(pd.read_csv(gtfs_directory + "route_info.csv")[["route_short_name", "interstate_intrastate"]], on="route_short_name")
buspositions = buspositions.merge(routes[["route_short_name", "interstate_intrastate"]], left_on="route_number", right_on="route_short_name", how="left")
sns.lineplot(data=buspositions, 
            x=buspositions.timestamp.dt.hour, 
            y="speed", 
            hue="interstate_intrastate")#,
            #style=buspositions.timestamp.dt.weekday.isin([0, 1, 2, 3, 4]))
# some heatmaps
sns.heatmap(pd.pivot_table(buspositions, values="speed", index="route_number", columns=buspositions.timestamp.dt.hour, aggfunc="mean"))
sns.heatmap(pd.pivot_table(buspositions, values="speed", index="route_number", columns=buspositions.timestamp.dt.weekday, aggfunc="mean"))
# doing facetgrids, faceted by route number
sns.relplot(data=buspositions,
            x=buspositions.timestamp.dt.hour, y="speed", col="route_number",
            hue=buspositions.timestamp.dt.weekday.isin([0, 1, 2, 3, 4]), kind="line"
            )

