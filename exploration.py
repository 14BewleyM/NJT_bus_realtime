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
from shapely.geometry import Point, LineString
import logging
import datetime
import gtfs_functions as gtfs
import webbrowser
import folium
import os
import re
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, DateTime, Boolean
from geoalchemy2 import Geometry

# setting some important directories and filepaths
main_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper"
code_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/code"
gtfs_directory = "C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs"
data_directory = "C:/Users/Bewle/OneDrive/Documents/school/Rutgers_courses_etc/2021_spring/3_transportation_equity/paper/data"
data_path = data_directory + "/data_deduped.csv"
data_path_bigger = data_directory + "/testdata_bigger.csv"
data_path_smaller = data_directory + "/testdata.csv"
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
def load_buspositions_data(projections, data_path=None, data_directory=None): 
    # provide data_path if want to read in data from cvs, data_directory if want to export from postgresql to csv
    # provide projections as a dict of projections containing an original projection, an equal area projection, and an equidistant projection
    # provide just projections if want to load with postgresql but not export to csv
    if projections is None:
        raise TypeError("Must provide a dict of projections containing an original projection, an equal area projection, and an equidistant projection")
    if data_path is None: # read in data from postgresql database
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
    
    
    elif data_path is not None: # read in deduped data from path
        print(f"Reading in data from {data_path}")
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

    # create time elapsed column
    print(f"Creating elapsed time column")
    buspositions["time_elapsed"] = buspositions.sort_values(by=["timestamp", "vehicle_id", "run_number"]).groupby(by=["date", "vehicle_id", "run_number"])["timestamp"].diff().fillna(pd.Timedelta(seconds=0))
    # testing for negative times
    negative_times = (buspositions.time_elapsed < datetime.timedelta(0)).sum()
    print(f"{negative_times} records calculated with negative elapsed time")
    # create timedelta columns in seconds
    buspositions["time_elapsed_seconds"] = buspositions.time_elapsed.apply(pd.Timedelta.total_seconds)

    return buspositions

# %% load gtfs
def load_headsigns(gtfs_directory):
    # load and headsigns and return a merged dataframe that's the union of headsigns from buspositions and gtfs
    # also return trips dataframe for later use

    # routes, stops, stop_times, trips, shapes = gtfs.import_gtfs(gtfs_directory)
    # load routes as dataframe
    print("Loading routes from gtfs")
    routes = gtfs.import_gtfs(gtfs_directory)[0]

    # load trips as dataframe
    print("Loading trips from gtfs")
    trips = pd.read_csv(gtfs_directory + "/trips.txt").astype(str)
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
    headsigns_merged = headsigns_buspositions.merge(headsigns_gtfs, how="outer", on="headsign", suffixes=("_buspositions", "_gtfs")).sort_values(by="headsign")
    ### pre-populate crosswalk field if the headsigns match
    headsigns_merged["headsign_crosswalk"] = headsigns_merged[(headsigns_merged.route_short_name_buspositions.notna()) & (headsigns_merged.route_short_name_gtfs.notna())].headsign
    ### merge in direction_id and shape_id

    return headsigns_merged, trips

# %% compare headsigns
def get_matched_routes(headsigns_merged, buspositions):
    # rows where headsign_crosswalk is not null are rows where buspositions headsigns can be matched to gtfs headsigns
    print(f"There are {headsigns_merged.headsign_crosswalk.notna().sum()} matched headsigns")
    # routes from buspositions all of whose rows have headsign_crosswalk not null are entirely matched
    print(f"There are {headsigns_merged.groupby(by='route_short_name_buspositions').headsign_crosswalk}")
    # routes that don't appear in a list of routes with no null values for headsign_crosswalk are routes that are entirely matched
    matched_routes = headsigns_merged[headsigns_merged.headsign_crosswalk.notna()].route_short_name_buspositions.unique()
    print(f"There are {len(matched_routes)} routes from buspositions whose headsigns have a matching gtfs headsign shape of {len(headsigns_merged.route_short_name_buspositions.unique())} total routes")

    buspositions_matched = buspositions[buspositions.route_number.isin(matched_routes)]
    print(f"Returning {buspositions_matched.shape[0]} observations from matched routes of {buspositions.shape[0]} total observations")
    return buspositions_matched, matched_routes

# %% construct final buspositions dataset
def create_shapes(matched_routes, headsigns_merged): # trips):
    print("L")
    headsign_shapes = trips[trips.route_short_name.isin(matched_routes)][["route_short_name", "trip_headsign", "direction_id"]]
    headsign_shapes = headsign_shapes.drop_duplicates()
    # edit headsigns in new dataframes to match buspositions format, using pattern from above
    headsign_shapes["trip_headsign"] = headsign_shapes["trip_headsign"].apply(lambda x: re.sub(r"[ ]*-[ ]*[Ee]x.*", "", x))
    headsign_shapes.trip_headsign = headsign_shapes.trip_headsign.str.replace("  ", " ")
    headsign_shapes.trip_headsign = headsign_shapes.trip_headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)
    # change dtypes to more appropriate types
    headsign_shapes[["direction_id", "shape_id"]] = headsign_shapes[["direction_id", "shape_id"]].astype(str)    

    # construct shapes from gtfs
    shapes = pd.read_csv(gtfs_directory + "/shapes.txt")
    shapes.shape_id = shapes.shape_id.astype(str)
    geometry_gtfs = gpd.points_from_xy(x=shapes.shape_pt_lon, y=shapes.shape_pt_lat) # create points from latlon
    shapes = gpd.GeoDataFrame(shapes, geometry=geometry_gtfs) # add points to create geodataframe

    # string points together into lines by shape_id
    shapes = shapes.sort_values(by="shape_pt_sequence").groupby(by=["shape_id"])["geometry"].apply(lambda x: LineString(x.tolist())).reset_index()



# including speed and matched headsigns

# %% load Census data

# %% create demographic measures

# %% analysis

# speed by route

# coefficient of variation

# %% figures
# see the heatmaps and line charts from gtfs_functions readme: https://github.com/Bondify/gtfs_functions



# %% main routine
def main():
    buspositions = load_buspositions_data(projections=projections_dict) # data_path=data_path)
    headsigns_merged, trips = load_headsigns(gtfs_directory)
    buspositions_matched, matched_routes = get_matched_routes(headsigns_merged, buspositions)
    print("Done")
if __name__ == "__main__":
    main()

# %% some initial settings
generate_maps = 0
if os.path.exists(data_directory + "/buspositions_deduped.csv"):
    print("Existing deduped csv file detected and will be read in")
    try_sql = 0
else:
    try_sql = 1
read_in_pandas = 0
export_bus_positions = 1
export_route_shapes = 0
generate_folium_maps = 0



# %% set up logging
logging.basicConfig(filename="log.log", 
                    level=logging.INFO) # set to INFO if you want info messages to log

logging.info(f"Logging started at {datetime.datetime.now()}")

# %% load data and setup dataframe
logging.info(f"Reading in data from {data_path}")
# read in in chunks for big file
chunk_size = 1.5 *  (10 ** 6) # number of rows to read in at a time (or once, as the case currently is)
buspositions = pd.DataFrame()

if try_sql == 1:
    # trying to dedupe dataset with sql before reading in, to reduce size 
    DeclarativeBase = declarative_base()

    # sqlalchemy documentation about types: https://docs.sqlalchemy.org/en/13/core/type_basics.html
    # one of these represents one vehicle position
    class Busposition(DeclarativeBase):
        __tablename__ = "buspositions"
        id = Column(Integer, primary_key=True)
        timestamp = sqlalchemy.Column(DateTime)
        vehicle_id = sqlalchemy.Column(String)
        route_number = sqlalchemy.Column(String)
        direction_compass = sqlalchemy.Column(String)
        #direction_compass_abbrev = sqlalchemy.Column(String)
        direction = sqlalchemy.Column(String)
        lonlat = sqlalchemy.Column(Geometry('POINT'))
        run_number = sqlalchemy.Column(String)
        headsign = sqlalchemy.Column(String)
        has_service = sqlalchemy.Column(Boolean)
        request_id = sqlalchemy.Column(Integer)
        lon = sqlalchemy.Column(String)
        lat = sqlalchemy.Column(String)

    # create connection to postgresql database with full dataset
    # create engine
    with open("database_pass.txt", "r") as file:
        lines = file.readlines()
        database_pass = lines[0]
    engine = sqlalchemy.create_engine(f"postgresql://postgres:{database_pass}@localhost:5432/bustest", echo=True)
    # create session
    Session = sqlalchemy.orm.sessionmaker(bind=engine)
    session = Session()

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

    print(f"Converting dataframe to geodataframe and converting to equidistant projection")
    buspositions = gpd.GeoDataFrame(buspositions, 
                        geometry=gpd.points_from_xy(buspositions.lon, buspositions.lat))
    buspositions.crs = projections_dict["original_projection"]
    buspositions = buspositions.to_crs(projections_dict["equidistant"])

elif read_in_pandas == 1:
    # read in with pandas
    cols = ["timestamp","vehicle_id", "route_number", "run_number", "headsign", "has_service", "lon", "lat"]
    with pd.read_csv(data_path, chunksize=chunk_size) as reader:
        iterator = 0
        for chunk in reader:
            print(f"Chunk {iterator}: Reading in chunk {iterator} with chunk size of {chunk_size} max rows")

            chunk = gpd.GeoDataFrame(chunk, 
                                geometry=gpd.points_from_xy(chunk.lon, chunk.lat))
            chunk.crs = projections_dict["original_projection"]

            # set timestamp column to be datetime and add day and month columns   
            chunk.timestamp = pd.to_datetime(chunk.timestamp)
            chunk["month"] = chunk.timestamp.dt.month
            chunk["day"] = chunk.timestamp.dt.day

            # set other data types
            chunk.vehicle_id = chunk.vehicle_id.astype(str).str.split(".").str[0]
            chunk.route_number = chunk.route_number.astype(str).str.split(".").str[0] 
            #chunk.has_service = chunk.has_service.map({"t":True, "f":False})

            # drop observations for routes that returned no service
            if chunk.has_service.dtype != bool:
                chunk.has_service = chunk.has_service.map({"t":True, "f":False})
            print(f"Chunk {iterator}: Dropping {chunk.shape[0] - chunk.has_service.sum()} observations from routes with no service")
            chunk = chunk[chunk.has_service==True]

            # some headsigns had ampersand replaced with &amp;amp;
            # put ampersand back
            print(f"Chunk {iterator}: Replacing ampersand in {~chunk.headsign.contains('amp;amp;').sum()} rows that have bad formatting")
            chunk.headsign = chunk.headsign.str.replace("amp;amp;", "")

            # check for duplicates by month, day, run number, vehicle id, and geometry
            # records with the same value for all of these probably reflect observations when a vehicle was moving but its position was not updated yet
            duplicate_count = chunk.duplicated(subset=["date", "run_number", "vehicle_id", "geometry"], keep="first").sum()
            print(f"Chunk {iterator}: There are {duplicate_count} observations with duplicate days, run numbers, vehicle ids, and geometries")

            # drop all except first duplicate for each vehicle
            print(f"Chunk {iterator}: Dropping {duplicate_count} duplicates")
            chunk = chunk.drop_duplicates(subset=["date", "run_number", "vehicle_id", "geometry"], keep="first")
            buspositions = buspositions.append(chunk)
            print(f"Chunk {iterator}: Dataset has {buspositions.shape[0]} records after dropping") # when run all together, this prints a value different from the actual final row count...

            iterator += 1
            if iterator == 3:
                break # heylisten


    buspositions = gpd.GeoDataFrame(buspositions, 
                        geometry=gpd.points_from_xy(buspositions.lon, buspositions.lat))
    buspositions.crs = projections_dict["original_projection"]
    buspositions = buspositions.to_crs(projections_dict["equidistant"])
    print(f"Initial dataset has {buspositions.shape[0]} records")

    # set timestamp column to be datetime and add day and month columns   
    buspositions.timestamp = pd.to_datetime(buspositions.timestamp)
    buspositions["month"] = buspositions.timestamp.dt.month
    buspositions["day"] = buspositions.timestamp.dt.day

    # set other data types
    buspositions.vehicle_id = buspositions.vehicle_id.astype(str).str.split(".").str[0]
    buspositions.route_number = buspositions.route_number.astype(str).str.split(".").str[0] 
    if buspositions.has_service.dtype != bool:
        buspositions.has_service = buspositions.has_service.map({"t":True, "f":False})

    # drop observations for routes that returned no service
    print(f"Dropping {buspositions.shape[0] - buspositions.has_service.sum()} observations from routes with no service")
    buspositions = buspositions[buspositions.has_service==True]

    # some headsigns had ampersand replaced with &amp;amp;
    # put ampersand back
    buspositions.headsign = buspositions.headsign.str.replace("amp;amp;", "")

    # check for duplicates by month, day, run number, vehicle id, and geometry
    # records with the same value for all of these probably reflect observations when a vehicle was moving but its position was not updated yet
    duplicate_count = buspositions.duplicated(subset=["date", "run_number", "vehicle_id", "geometry"], keep="first").sum()
    print(f"There are {duplicate_count} observations with duplicate days, run numbers, vehicle ids, and geometries")

    # drop all except first duplicate for each vehicle
    print(f"Dropping {duplicate_count} duplicates")
    #buspositions = buspositions[~buspositions.sort_values(by="timestamp").duplicated(subset=["month", "day", "run_number", "vehicle_id", "geometry"], keep="first")]
    buspositions = buspositions.drop_duplicates(subset=["date", "run_number", "vehicle_id", "geometry"], keep="first")
    print(f"Final dataset has {buspositions.shape[0]} records")

else: # if not reading in with sql or anew with pandas, read in from existing csv
    print("Reading in deduped buspositions csv")
    buspositions = pd.read_csv(data_directory + "/buspositions_deduped.csv")

# create time elapsed column
print(f"Creating elapsed time column")
buspositions["time_elapsed"] = buspositions.sort_values(by=["timestamp", "vehicle_id", "run_number"]).groupby(by=["date", "vehicle_id", "run_number"])["timestamp"].diff().fillna(pd.Timedelta(seconds=0))
# testing for negative times
negative_times = (buspositions.time_elapsed < datetime.timedelta(0)).sum()
print(f"{negative_times} records calculated with negative elapsed time")
# create timedelta columns in seconds
buspositions["time_elapsed_seconds"] = buspositions.time_elapsed.apply(pd.Timedelta.total_seconds)

if export_bus_positions == 1:
    # drop timedelta column (bc can only be exported to geopackage as seconds, not as timedelta type)
    os.chdir(data_directory)
    if not os.path.exists(data_directory + "/buspositions_cleaned.gpkg"):
        buspositions.drop(columns="time_elapsed").astype({"date": str}).to_file("buspositions_cleaned.gpkg", driver="GPKG")
    if not os.path.exists(data_directory + "/buspositions_deduped.csv"):
        buspositions.to_csv("buspositions_deduped.csv")
    os.chdir(code_directory)

# %% load gtfs (including route shapes) and construct dataframes of headsigns from both datasets
# routes, stops, stop_times, trips, shapes = gtfs.import_gtfs(gtfs_directory)
routes = gtfs.import_gtfs(gtfs_directory)[0]
# segments_gdf = gtfs.cut_gtfs(stop_times, stops, shapes)

# function to replace within each item in a list
# so that you can apply a it to a whole series of lists
def replace_in_list(patt, repl, list):
    newlist = []
    for element in list: 
        newlist.append(re.sub(patt, repl, element))
        #list[element] = element
    return newlist

# if they have headsigns in common, it'll be much easier to make the vehicle positions to GTFS service patterns
trips_full = pd.read_csv("C:/Users/Bewle/OneDrive/Documents/data/geographic/NJT/NJT_bus_gtfs/trips.txt")
trips_full = trips_full.astype(str)
trips_full = pd.merge(trips_full, routes[["route_id", "route_short_name"]], on="route_id")

# number of unique headsigns per route from gtfs given by:
unique_headsigns_gtfs = trips_full.groupby(by=["route_short_name"])["trip_headsign"].unique().apply(len).reset_index().sort_values(by="route_short_name").rename(columns={"trip_headsign": "headsign_count"})
unique_headsigns_gtfs["headsigns"] = trips_full.groupby(by=["route_short_name"])["trip_headsign"].unique().reset_index()["trip_headsign"] # column to hold lists of headsigns
unique_headsigns_gtfs.headsigns = unique_headsigns_gtfs.headsigns.apply(lambda x: replace_in_list("  ", " ", x)) # replace double spaces with single spaces
trips_full.groupby(by="route_id")["trip_headsign"].unique().apply(len).sum() # total
# by direction and headsign
trips_full.groupby(by=["direction_id", "trip_headsign"]).size().reset_index().rename(columns={0:"count"})

# headsigns from vehicle measurements
unique_headsigns_buspositions = buspositions.groupby(by="route_number")["headsign"].unique().apply(len).reset_index().sort_values(by="route_number").rename(columns={"headsign": "headsign_count"})
unique_headsigns_buspositions["headsigns"] = buspositions.groupby(by=["route_number"])["headsign"].unique().reset_index()["headsign"]
unique_headsigns_buspositions.headsigns = unique_headsigns_buspositions.headsigns.apply(lambda x: replace_in_list("  ", " ", x)) # replace double spaces with single spaces

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

unique_headsigns_gtfs["headsigns"] = unique_headsigns_gtfs["headsigns"].apply(lambda x: replace_in_list(r"[ ]*-[ ]*[Ee][x].*", "", x))
# new dataframe merging headsigns etc from both datasets
unique_headsigns_buspositions = unique_headsigns_buspositions.rename(columns={"route_number": "route_short_name"})
headsigns_merged = pd.merge(unique_headsigns_buspositions, unique_headsigns_gtfs, on="route_short_name", suffixes=("_buspositions", "_gtfs"))


# %% check for headsign crosswalk and if exists use it to replace headsigns in both datasets

# create dataframe where each row is a headsign (rather than a list of headsigns)
headsigns_split_gtfs = trips_full[["trip_headsign", "route_short_name"]].drop_duplicates(subset=["trip_headsign", "route_short_name"]).rename(columns={"trip_headsign":"headsign"})
headsigns_split_gtfs.headsign = headsigns_split_gtfs.headsign.str.replace("  ", " ")
headsigns_split_gtfs.headsign = headsigns_split_gtfs.headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)

headsigns_split_buspositions = buspositions[["headsign", "route_number"]].drop_duplicates(subset=["headsign", "route_number"]).rename(columns={"route_number": "route_short_name"})
headsigns_split_buspositions.headsign = headsigns_split_buspositions.headsign.str.replace("  ", " ")
headsigns_split_buspositions.headsign = headsigns_split_buspositions.headsign.str.replace(r"[ ]*-[ ]*[Ee][x].*", "", regex=True)

# merge them to find which headsigns match and which don't
headsigns_split_merged = headsigns_split_buspositions.merge(headsigns_split_gtfs, how="outer", on="headsign", suffixes=("_buspositions", "_gtfs")).sort_values(by="headsign")
# pre-populate crosswalk field if the headsigns match
headsigns_split_merged["headsign_crosswalk"] = headsigns_split_merged[(headsigns_split_merged.route_short_name_buspositions.notna()) & (headsigns_split_merged.route_short_name_gtfs.notna())].headsign
headsigns_split_merged.to_excel(data_directory + "/headsigns_merged.xlsx")
# manually create the rest of the crosswalk values
# basically, copy matching values from the gtfs to buspositions, so that the busposition headsigns can be matched to the gtfs route shapes

#for index, row in unique_headsigns_gtfs.iterrows():
 #   temp_list.append(list(pd.Series(row.headsigns)))
  #  print(f"Appending series of headsigns for route {row.route_short_name}")
    #headsigns_split_gtfs.headsigns = headsigns_split_gtfs.headsigns.append(temp_series, ignore_index=True)
#headsigns_split_gtfs.headsigns = headsigns_split_gtfs.headsigns.append(pd.Series(temp_list))
#unique_headsigns_gtfs.apply(lambda x: pd.Series(x.headsigns))
#pd.Series(unique_headsigns_buspositions.headsigns.to_list()[0])


# %% compare busposition and gtfs headsigns
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
print(f"Dropping {buspositions.route_number.isin(matched_route_numbers).sum()} observations from buspositions whose headsigns don't exactly match a GTFS headsign")
buspositions = buspositions[buspositions.route_number.isin(matched_route_numbers)]
print(f"Dataset has {buspositions.shape[0]} records remaining after dropping")
headsign_shapes = trips_full[trips_full.route_short_name.isin(matched_route_numbers)][["route_short_name", "trip_headsign", "direction_id", "shape_id"]]
headsign_shapes = headsign_shapes.drop_duplicates()
# edit headsigns in new dataframes to match buspositions format, using pattern from above
headsign_shapes["trip_headsign"] = headsign_shapes["trip_headsign"].apply(lambda x: re.sub(r"[ ]*-[ ]*[Ee]x.*", "", x))
# change dtypes to more appropriate types
headsign_shapes[["direction_id", "shape_id"]] = headsign_shapes[["direction_id", "shape_id"]].astype(str)
# there are still multiple shapes per headsign
# maybe (hopefully) these are just different parts of what is basically a single shape
# in that case, can just merge them to create one shape per headsign (assuming each headsign corresponds to only one direction)

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
headsign_shapes.crs = projections_dict["original_projection"]
headsign_shapes = headsign_shapes.to_crs(projections_dict["equidistant"])
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
if export_route_shapes == 1:
    os.chdir(main_directory)
    headsign_shapes.to_file("headsign_shapes.gpkg", driver="GPKG")
    headsign_shapes_dissolved.to_file("headsign_shapes_dissolved.gpkg", driver="GPKG")
    headsign_shapes_dissolved_bydirection.to_file("headsign_shapes_dissolved_bydirection.gpkg", driver="GPKG")
    os.chdir(code_directory)


# TODO@14BewleyM subset buspositions dataset so that you use only observations from routes in the final headsigns dataframe


# interpolate busposition points to lines based on headsign and shape ids
# https://stackoverflow.com/questions/33769860/pandas-apply-but-only-for-rows-where-a-condition-is-met
# using interpolate: https://gis.stackexchange.com/questions/306838/snap-points-shapefile-to-line-shapefile-using-shapely
buspositions = buspositions.merge(headsign_shapes_dissolved_bydirection[["geometry", "trip_headsign", "shape_id", "direction_id"]].rename(columns={"trip_headsign": "headsign"}), on="headsign", how="left") # try to merge in the line shapes (had hoped to avoid doing for memory reasons, but may be the simplest way)
buspositions = buspositions.rename(columns={"geometry_x": "point_original", "geometry_y": "route_shape"}) # distinguish new geometry columns from each other
# the CRS comes undone after the last two operations, so I'm not sure what it uses for the next few, which would seem to require one
buspositions = gpd.GeoDataFrame(buspositions, crs=projections_dict["equidistant"], geometry="point_original")
#buspositions[["route_shape", "point_original"]] = buspositions[["route_shape", "point_original"]].apply(lambda col: gpd.GeoSeries(col).set_crs(projections_dict["equidistant"]))
buspositions["point_interpolated"] = buspositions.apply(lambda row: row.route_shape.interpolate(row.route_shape.project(row.point_original)), axis=1)
# geometry columns already projected to ESPG 3424 (NJ state plane with units of ft) for distance calculations: https://www.spatialreference.org/ref/epsg/3424/
buspositions["distance_traveled_cumul"] = buspositions.apply(lambda row: row.route_shape.project(row.point_interpolated, normalized=False), axis=1) # gotta be a way to do this in one call to .apply(), create a series or something
# this gives cumulative distance as proportion of the overall route shape
# need distance in feet (multiply by total route shape length)
# wait, no, the normalized argument is false by default, so this should be in feet...
# that the vast majority of measurements are 1 or below makes me think this is normalized, whatever the default arg is supposed to be
# goddammit, run again with normalized=False explicitly stated
buspositions = gpd.GeoDataFrame(buspositions.drop(columns=["route_shape"]), geometry="point_interpolated", crs=projections_dict["equidistant"])


#%% calculate distance covered between measurements
# may want to use or consult gtfs_functions' cut_gtfs at some point: https://github.com/Bondify/gtfs_functions/blob/master/gtfs_functions/gtfs_funtions.py
# shift within each vehicle number (and run number, in case a vehicle goes out of service and then reappears a long time or distance later?)
# link: https://gis.stackexchange.com/questions/347732/finding-distance-of-consecutive-points-in-a-geopandas-data-frame

## some more links about methods for calculating distance (both from start of shape and between points)
### from start of shape: https://gis.stackexchange.com/questions/280112/points-layer-distance-from-the-start-of-line-layer-in-qgis
### between points: https://gis.stackexchange.com/questions/184554/measure-distance-between-points-along-a-line-in-qgis
### if possible, probably best to measure between points, because merging the route shapes by headsign and direction may still produce strange routings that don't represent routes real-world vehicles take

# %% analyze vehicle and run numbers see if vehicles typically go in and out of service in a way that leaves long time or distance gaps

# vehicles with the most observations
buspositions.groupby(by="vehicle_id").size().sort_values()
buspositions[buspositions.vehicle_id=="5946"].sort_values(by="timestamp", ascending=True)
# how many run on more than one route? few
number_vehicles_total = len(buspositions.vehicle_id.unique())
vehicles_byroute_summed_desc = buspositions.groupby(by="vehicle_id").route_number.unique().apply(len).describe()
vehicles_multiple_routes = buspositions.groupby(by="vehicle_id").route_number.unique().apply(len) > 1
print(f"There are {number_vehicles_total} distinct vehicle ids")
print(f"{vehicles_multiple_routes.sum()} vehicles are associated with >1 route")
vehicles_byrun_summed_desc = buspositions.groupby(by="vehicle_id").run_number.unique().apply(len).describe()
vehicles_multiple_runs = buspositions.groupby(by="vehicle_id").run_number.unique().apply(len) > 1
print(f"{vehicles_multiple_runs.sum()} vehicles are associated with >1 run")
print(f"{(vehicles_multiple_runs & vehicles_multiple_routes).sum()} vehicles are associated with >1 route and >1 run")
print(f"{(~vehicles_multiple_runs & ~vehicles_multiple_routes).sum()} vehicles are associated with only 1 route and only 1 run")

# basically, it's safe to say there are a significant number of vehicles that do multiple runs
# so should at least group by vehicle id and run number
# but should also group by headsign, bc there are sometimes gaps in time or distance where the run number doesn't change but the headsign does
# (eg vehicle 5946 on route 1 with run number 653 changes headsign but not run number btw observations 958877 and 983363, with a gap of about 16 minutes)
# this will exclude some relevant observations (bc you won't be able to calculate differences and speeds for some rows that become the first in a group)
# but that's better than having some observations that would be useful but have to be discarded bc of very high speeds, and not in a sensible way

buses_sorted_grouped = buspositions.sort_values(by="timestamp").groupby(by=["date", "vehicle_id", "headsign", "run_number"])
# reproject into ESPG 3424 (NJ state plane with units of ft) for distance calculations: https://www.spatialreference.org/ref/epsg/3424/

# get distance along line of each observation in group (using .project())

# get difference between distance associated with each observation and its previous observation (.diff() on the relevant column should be enough)
buspositions["distance_traveled_prev"] = buses_sorted_grouped.distance_traveled_cumul.diff()


# %% calcuate speeds
# calculate speed values in mph
buspositions["speed"] = (buspositions.distance_traveled_prev/5280) / (buspositions.time_elapsed_seconds/ (60 *60))
# passes smell test, but there are some absurd values

# calculate max distance traveled and time taken by a vehicle for each headsign and run combination
# subtract first observation (or timestamp.min()) from last observation (or timestamp.max()) within each group
# but where to assign this? this may be more an analysis than a data processing step

# identify rows with very low speeds
# (there are some measurements with very large time differences, like several minutes)
# (some even in the thousands of seconds - those have gotta be buses just sitting, right?)
# (if you exclude those measurements, you should probably recalculate all the times, right?)

# how many observations have negative speeds? by route
neg_speeds_byroute = buspositions[buspositions.speed < 0].groupby(by="route_number").size().reset_index().rename(columns={0: "neg_speeds_count"}) 
neg_speeds_byroute["total_obs_count"] = buspositions.groupby(by="route_number").size().reset_index()[0]
neg_speeds_byroute["neg_speed_pct"] = neg_speeds_byroute.neg_speeds_count / neg_speeds_byroute.total_obs_count
neg_speeds_byroute.sort_values(by="neg_speed_pct", ascending=False)
print(f"There are {(buspositions.speed < 0).sum()} observations with negative speeds, out of {buspositions.shape[0]} total observations")
print(f"There are {buspositions.speed.isna().sum()} observations with null speeds, out of {buspositions.shape[0]} total observations")
print(f"These represent {((((buspositions.speed < 0).sum() + buspositions.speed.isna().sum())/buspositions.shape[0]) * 100):.{3}}% of all observations")
# the negative speeds are all due to negative distance measures
(buspositions.time_elapsed_seconds < 0).sum()
(buspositions.distance_traveled_prev < 0).sum()
# this makes me think the negative speeds are legit
# and just the result of the shapes being incorrectly oriented 
buspositions[buspositions.speed < -60] #~2400 observations with speeds less than 60mph

# for now, just exclude the observations that are null or have speeds less than -60mph (judging those to be reasonable measurements)
# and convert remaining negative speeds to positive
# you'll exclude nulls in any case bc those should be the first observation within each group
print(f"Dropping observations with null speed or with negative speeds less than -60mph: {buspositions[(buspositions.speed < -60) | (buspositions.speed.isna())].shape[0]} observations")
buspositions = buspositions[(buspositions.speed >= -60) & (~buspositions.speed.isna())]
print(f"Converting {buspositions[buspositions.speed < 0].shape[0]} remaining negative speed observations")
buspositions.loc[(buspositions.speed < 0), "speed"] = buspositions.speed * -1
print(f"Resulting dataset has {buspositions.shape[0]} records")

# average speeds by route
avg_speeds_by_route = buspositions.groupby(by="route_number").speed.mean().reset_index()
avg_speeds_by_route.to_csv("avg_speeds_byroute.csv")
# add scheduled speeds (calculated using gtfs_functions' .speed_from_gtfs())
# speeds = gtfs.speeds_from_gtfs(routes, stop_times, segments_gdf)
# avg_speeds_by_route["scheduled_speed_avg"] = speeds[speeds.route.isin(matched_route_numbers)].groupby(by="route").speed_mph.mean()
# you can also specify windows of time as a list of breakpoints and pass it to .speeds_from_gtfs() as an arg "cutoffs"

# %% some basic summary

first_observation = buspositions.timestamp.min()
last_observation = buspositions.timestamp.max()
logging.info(f"Dataset spans {first_observation} to {last_observation}")

unique_routes = buspositions.route_number.unique()
logging.info(f"Number of unique routes: {len(unique_routes)}")


# %% setting up Census data

#acs_tract = pd.read_csv(data_path_tract)
acs_blockgroups = pd.read_csv(data_path_blockgroup)
old_cols = list(acs_blockgroups.columns)
new_cols = [new_col.lower() for new_col in old_cols] 
acs_blockgroups = acs_blockgroups.rename(columns={old_cols[i]: new_cols[i] for i in range(len(old_cols))})
#tracts = gpd.read_file(boundary_path_tract)
blockgroups = gpd.read_file(boundary_path_blockgroup)
old_cols = list(blockgroups.columns)
new_cols = [new_col.lower() for new_col in old_cols] 
blockgroups = blockgroups.rename(columns={old_cols[i]: new_cols[i] for i in range(len(old_cols))})
# prefer them lowercase

# create merged dataframe associating block group shapes with their 2015-2019 ACS estimates
blockgroups_merged = blockgroups.merge(acs_blockgroups, on="gisjoin")

# projection is USA_Contiguous_Albers_Equal_Area_Conic (EPSG 102003 or 5070), suitable for buffering bc of equal area

# create headsign shapes for spatial join with Census data
# declare buffer size 
buffersize = 402.336 # 1/4 mile in meters
# create headsign shapes with bufffer 
buffer = headsign_shapes_dissolved_bydirection.to_crs(projections_dict["equal_area"]).geometry.buffer(buffersize)
headsign_shapes_dissolved_bydirection["buffer"] = buffer

# join block groups that intersect the headsign shape buffers
blockgroups_headsigns_joined = gpd.sjoin(headsign_shapes_dissolved_bydirection.set_geometry("buffer").to_crs(projections_dict["equal_area"]), blockgroups_merged.to_crs(projections_dict["equal_area"]), op="intersects", how="left")

# %% calculate equity measures

# see for some examples:
# Boston, p 11: https://d3n8a8pro7vhmx.cloudfront.net/livablestreetsalliance/pages/6582/attachments/original/1569205099/lsa-better-buses-2019-v9-20sep19.pdf?1569205099


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

# function to create columns for values of interest
def create_equity_measures(df, by=None, aggfunc=None):
    """
    Takes in a dataframe, and optionally fields to group by and an aggregation function.
    Calculates the equity measures.
    Returns pandas series with the equity measures (counts and proportions). 
    """
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

# create new columns for values of interest
for df in [blockgroups_headsigns_joined, acs_blockgroups]:
    df["prop_poc"], df["prop_white"], df["prop_latino_allraces"], df["prop_black"], df["prop_asian"], df["prop_poverty"] = create_equity_measures(df)
    #df["prop_poc"] = (df.aluke001 - df.aluke003) / df.aluke001
    #df["prop_white"] = df.aluke003 / df.aluke001
    #df["prop_latino_allraces"] = df.aluke012 / df.aluke001
    #df["prop_black"] = df.aluce003 / df.aluce001
    #df["prop_asian"] = df.aluce005 / df.aluce001
    #df["prop_poverty"] = df.alwye002 / df.alwye001

# define new dataframes for each route and for the whole service area, to avoid double-counting
# for each route (dropping block groups that are associated more than once with a route)
blockgroups_routes_joined_buffer = blockgroups_headsigns_joined.drop_duplicates(["route_short_name", "gisjoin"])
# for the whole service area (dropping all duplicated block groups, using gisjoin field as unique id for tracts)
blockgroups_servicearea_joined = blockgroups_headsigns_joined.drop_duplicates(["gisjoin"])

# service area averages for each of the measures
# index is each of the values, access by eg service_area_distribution.prop_poc.loc["mean"]
service_area_distribution = blockgroups_servicearea_joined[["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]].describe()
# add statewide mean for each of the variables
service_area_distribution.loc["state_mean"] = acs_blockgroups[["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]].mean()

# calculating measures for the whole service area and for each route
route_index = list(blockgroups_servicearea_joined.route_short_name.unique())
route_index.append("service_area")
columns = ["count_total", "count_poc", "prop_poc", "count_white", "prop_white", "count_latino_allraces", "prop_latino_allraces", "count_black", "prop_black", "count_asian", "prop_asian", "count_poverty_hholds", "prop_poverty"]
equity_measures_byroute = pd.DataFrame(index=route_index, columns=columns)
# calculate proportion measures with function (do this also for the merged blockgroups dataframe, for later use in calculating Title VI measures)
equity_measures_byroute["prop_poc"], equity_measures_byroute["prop_white"], equity_measures_byroute["prop_latino_allraces"], equity_measures_byroute["prop_black"], equity_measures_byroute["prop_asian"], equity_measures_byroute["prop_poverty"] = create_equity_measures(blockgroups_routes_joined_buffer, by="route_short_name", aggfunc=np.sum)
blockgroups_merged["prop_poc"], blockgroups_merged["prop_white"], blockgroups_merged["prop_latino_allraces"], blockgroups_merged["prop_black"], blockgroups_merged["prop_asian"], blockgroups_merged["prop_poverty"] = create_equity_measures(blockgroups_merged)
# calculate proportion measures for the entire service area and statewide, and add to dataframe
equity_measures_byroute.loc["service_area", ["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]] = list(create_equity_measures(blockgroups_servicearea_joined, aggfunc=np.sum))
equity_measures_byroute.loc["statewide", ["prop_poc", "prop_white", "prop_latino_allraces", "prop_black", "prop_asian", "prop_poverty"]] = list(create_equity_measures(acs_blockgroups, aggfunc=np.sum))
# calculate total measures by route, as well as for the entire service area and statewide, and add to dataframe
totals_map_dict = {"count_total": "alube001",
                    "count_white": "aluke003",
                    "count_latino_allraces": "aluke012",
                    "count_black": "aluce003", 
                    "count_asian": "aluce005",
                    "count_poverty_hholds": "alwye001"}
for column in totals_map_dict:
    #print(column)
    #print(original_variable)
    equity_measures_byroute[column] = blockgroups_routes_joined_buffer.groupby(by="route_short_name")[totals_map_dict[column]].sum()
    equity_measures_byroute.loc["service_area", column] = blockgroups_servicearea_joined[totals_map_dict[column]].sum()
    equity_measures_byroute.loc["statewide", column] = acs_blockgroups[totals_map_dict[column]].sum()
# calculate POC count separately bc you have to subtract
equity_measures_byroute["count_poc"] = blockgroups_routes_joined_buffer.groupby(by="route_short_name").aluke001.sum() - blockgroups_routes_joined_buffer.groupby(by="route_short_name").aluke003.sum()
equity_measures_byroute.loc["service_area", "count_poc"] = blockgroups_servicearea_joined.aluke001.sum() - blockgroups_servicearea_joined.aluke003.sum()
equity_measures_byroute.loc["statewide", "count_poc"] = acs_blockgroups.aluke001.sum() - acs_blockgroups.aluke003.sum()

# POC proportion is higher than the ~0.44 reported in 2017 Title VI report for the service area
# maybe your 1/4 mile buffer is narrower than the one they use

# calculate FTA equity measures
# routes 1/3 of whose revenue miles run through "equity neighborhoods", which are those where the percentage minority population is greater than the percentage for the service area
service_area_poc_prop = equity_measures_byroute.loc["service_area", "prop_poc"]
service_area_black_prop = equity_measures_byroute.loc["service_area", "prop_black"]
service_area_latino_prop = equity_measures_byroute.loc["service_area", "prop_latino_allraces"]
service_area_asian_prop = equity_measures_byroute.loc["service_area", "prop_asian"]
# designate block groups as equity block groups if their proportion is greater than relevant service area proportion
blockgroups_routes_joined_buffer["equity_block_group"] = False
blockgroups_routes_joined_buffer.loc[blockgroups_routes_joined_buffer.prop_poc > service_area_poc_prop, "equity_block_group"] = True
blockgroups_merged["equity_block_group"] = False
blockgroups_merged.loc[blockgroups_merged.prop_poc > service_area_poc_prop, "equity_block_group"] = True
# redesignate geometry and reproject to equidistant projection for distance measurements
# need to calculate a new intersection using only the line shapes, not their buffer
blockgroups_routes_joined_line = gpd.sjoin(headsign_shapes_dissolved_bydirection.to_crs(projections_dict["equidistant"]), blockgroups_merged.to_crs(projections_dict["equidistant"]), op="intersects", how="left")
# find ratio of route length within equity block groups to route length outside of equity block groups
# for each route, dissolve on equity_block_group field, then calculate length
blockgroups_routes_joined_line = blockgroups_routes_joined_line.dissolve(by=["route_short_name", "equity_block_group"]) # maybe add as_index=False if don't want groupby fields to be new index
blockgroups_routes_joined_line["length"] = blockgroups_routes_joined_line.geometry.length
blockgroups_routes_joined_line.groupby(by=["route_short_name", "equity_block_group"]).length.sum().reset_index()
# add to equity_measures_byroute dataframe as new columns, "rev_miles_total", "rev_miles_equity", "rev_miles_equity_ratio"
equity_measures_byroute.loc[equity_measures_byroute.reset_index().route == blockgroups_routes_joined_line.route_short_name, "rev_miles_total"] = blockgroups_routes_joined_line


# TODO@14BewleyM calculate speed variability, assign to block groups, and assess 