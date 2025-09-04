import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from functions import helper_fcts
import geopandas as gpd

# prepare data
data_path = "D:/data/CAMELS_DE/"

# check if folders exist
results_path = "results/"
if not os.path.isdir(results_path):
    os.makedirs(results_path)
figures_path = "figures/"
if not os.path.isdir(figures_path):
    os.makedirs(figures_path)
figures_timeseries_path = "figures/timeseries/"
if not os.path.isdir(figures_path):
    os.makedirs(figures_path)

# load attributes needed for quality control
df_topo = pd.read_csv(data_path + "CAMELS_DE_topographic_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
df_climate = pd.read_csv(data_path + "CAMELS_DE_climatic_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
df_landcover = pd.read_csv(data_path + "CAMELS_DE_landcover_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
df_humaninfluence = pd.read_csv(data_path + "CAMELS_DE_humaninfluence_attributes.csv", sep=',', skiprows=0, encoding='latin-1')

df_attr = pd.merge(df_topo, df_climate, on='gauge_id')
df_attr = pd.merge(df_attr, df_landcover, on='gauge_id')
df_attr = pd.merge(df_attr, df_humaninfluence, on='gauge_id')

'''
Level 1 Network:
-  ≤10% urbanization in catchment
-  Minimal human impact on flows
-  No major land use changes affecting streamflow
-  Very high-quality daily mean flow data
-  Reliable high and low flow representation
-  At least 40 years of records
-  No data gaps >3 years

Level 2 Network:
-  ≤20% urbanization in catchment
-  Modest impact on monthly and annual flows
-  High to fair quality daily mean flow data
-  Reliable monthly average flow representation
-  At least 20 years of records
-  No specific requirement for data gaps
'''

gauge_id_list = []
mean_P_list = []
mean_PET_list = []
mean_Q_list = []
perc_complete_list = []
record_length_list = []
data_gap_list = []
data_flag_list = []

for id in df_attr["gauge_id"]:
    # print(id)

    # load data
    df_tmp = pd.read_csv(data_path + "timeseries/CAMELS_DE_hydromet_timeseries_" + str(id) + ".csv", sep=',')
    df_tmp["date"] = pd.to_datetime(df_tmp["date"])
    df_sim = pd.read_csv(data_path + "timeseries_simulated/CAMELS_DE_discharge_sim_" + str(id) + ".csv", sep=',')
    df_sim["date"] = pd.to_datetime(df_sim["date"])
    df_tmp = pd.merge(df_tmp, df_sim[["date", "pet_hargreaves"]], on='date')

    # remove NaNs at beginning
    first_valid_index = df_tmp["discharge_spec_obs"].first_valid_index()
    df_tmp = df_tmp.loc[first_valid_index:]

    # check data
    perc_complete = np.sum(~np.isnan(df_tmp["discharge_spec_obs"].values)) / len(df_tmp["discharge_spec_obs"].values)
    record_length = (df_tmp["date"].max() - df_tmp["date"].min()).days / 365.25
    df_copy = df_tmp.copy().dropna(subset=["discharge_spec_obs"])
    data_gap = df_copy["date"].diff().max()
    if df_tmp["discharge_spec_obs"].any() < 0:
        print("Negative values in streamflow data: ", id)
    if df_tmp["discharge_spec_obs"].any() > 1000:
        print("Unrealistic high values in streamflow data (> 1000 mm/d): ", id)

    # add data
    gauge_id_list.append(id)
    perc_complete_list.append(perc_complete)
    record_length_list.append(record_length)
    data_gap_list.append(data_gap)
    mean_P_list.append(np.nanmean(df_tmp["precipitation_mean"]))
    mean_Q_list.append(np.nanmean(df_tmp["discharge_spec_obs"]))
    mean_PET_list.append(np.nanmean(df_tmp["pet_hargreaves"]))

# create dataframe
df = pd.DataFrame()
df["gauge_id"] = gauge_id_list  # as a check
df = pd.merge(df_attr, df, on='gauge_id')
df["mean_P"] = mean_P_list
df["mean_Q"] = mean_Q_list
df["mean_PET"] = mean_PET_list
df["runoff_ratio"] = df["mean_Q"] / df["mean_P"]
df["aridity"] = df["mean_PET"] / df["mean_P"]
df["perc_complete"] = perc_complete_list
df["record_length"] = record_length_list
df["data_gap"] = data_gap_list
df["area_error"] = np.abs((df["area_metadata"] - df["area"]) / df["area_metadata"])

# quality control
# flags: level 1 = 1, level 2 = 2, did not pass checks = 0

# human impacts
df["humanimpact_flag"] = 0
df.loc[(df["artificial_surfaces_perc"] >= 0) & (df["artificial_surfaces_perc"] <= 10) & (
        df["dams_num"] < 1), "humanimpact_flag"] = 1
df.loc[(df["artificial_surfaces_perc"] > 10) & (df["artificial_surfaces_perc"] <= 20) & (
        df["dams_num"] < 1), "humanimpact_flag"] = 2
# land use change is currently not possible to check

# data quality
df["recordlength_flag"] = 0
df.loc[(df["record_length"] >= 40) & (df["perc_complete"] >= 0.9) & (
        df["data_gap"] < pd.Timedelta(days=1095)), "recordlength_flag"] = 1
df.loc[(df["record_length"] >= 20) & (df["record_length"] < 40) & (df["perc_complete"] >= 0.9) & (
        df["data_gap"] < pd.Timedelta(days=1095)), "recordlength_flag"] = 2
# NOTE: information on individual gauging stations, rating curve uncertainties, etc. not available

# other checks
df["dataquality_flag"] = 0
df.loc[(df["area_error"] < 0.1) & (df["mean_P"] > df["mean_Q"]), "dataquality_flag"] = 1
# check for strange Q values: nan, negative, many 0s, consecutive days with same flow, step changes -> ROBIN CHECKS?

# ROBIN list
df_ROBIN = pd.read_csv("camel_robin_overlap.csv", sep=',', skiprows=0, encoding='latin-1')
df["ROBIN_flag"] = 0
df.loc[df["gauge_id"].isin(df_ROBIN["gauge_id"]), "ROBIN_flag"] = 1

# final data flag
df["data_flag"] = 0
df.loc[((df["humanimpact_flag"] == 1) | (df["humanimpact_flag"] == 2)) &
       ((df["recordlength_flag"] == 1) | (df["recordlength_flag"] == 2)) &
       ((df["dataquality_flag"] == 1) | (df["dataquality_flag"] == 2)) &
       (df["ROBIN_flag"] == 1), "data_flag"] = 2
df.loc[(df["humanimpact_flag"] == 1) & (df["recordlength_flag"] == 1) & (df["dataquality_flag"] == 1) & (
        df["ROBIN_flag"] == 1), "data_flag"] = 1

# create dataframe that only consists of "data_flag" not equal to 0
df_checked = df[df["data_flag"] != 0]

# after initial screening, save all timeseries as plots and check manually again
'''
for id in df_checked["gauge_id"]:
    # load data
    df_tmp = pd.read_csv(data_path + "timeseries/CAMELS_DE_hydromet_timeseries_" + str(id) + ".csv", sep=',')
    df_tmp["date"] = pd.to_datetime(df_tmp["date"])

    # plot data
    fig, ax = plt.subplots(figsize=(12, 4), tight_layout=True)
    im = ax.plot(df_tmp["date"], df_tmp["discharge_spec_obs"], alpha=0.9)
    ax.set_ylabel("Streamflow (mm/d)")
    fig.savefig(figures_timeseries_path + "CAMELS_DE_" + id + ".png", dpi=600, bbox_inches='tight')
    plt.close()
'''

# extract attributes that are necessary and save results
'''
    ID
    Name
    Longitude
    Latitude
    Catchment Area (km²)
    Measuring Organisation
    Level 1 or Level 2 (Catchment Development)
    Level 1 or Level 2 (Data Quality)
    Record Length
    Missing Data Criteria Met? (Yes/No)
'''

df_checked["ID"] = df_checked["gauge_id"]
df_checked["Name"] = df_checked["gauge_name"]
df_checked["Longitude"] = df_checked["gauge_lon"]
df_checked["Latitude"] = df_checked["gauge_lat"]
df_checked["Catchment Area (km²)"] = df_checked["area"]
df_checked["Measuring Organisation"] = df_checked["federal_state"]
df_checked["Level 1 or Level 2 (Catchment Development)"] = df_checked["humanimpact_flag"]
df_checked["Level 1 or Level 2 (Data Quality)"] = df_checked["data_flag"]
df_checked["Record Length"] = np.round(df_checked["record_length"], 2)
df_checked["Missing Data Criteria Met? (Yes/No)"] = df_checked["data_gap"].apply(
    lambda x: "Yes" if x < pd.Timedelta(days=1095) else "No")

# count how many in df_final Level 1 or Level 2 are Level 1 or Level 2
print(df_checked["Level 1 or Level 2 (Catchment Development)"].value_counts())
print(df_checked["Level 1 or Level 2 (Data Quality)"].value_counts())

df_final = df_checked[["ID", "Name", "Longitude", "Latitude", "Catchment Area (km²)", "Measuring Organisation",
                       "Level 1 or Level 2 (Catchment Development)", "Level 1 or Level 2 (Data Quality)",
                       "Record Length", "Missing Data Criteria Met? (Yes/No)",
                       "gauge_name", "water_body_name", "artificial_surfaces_perc", "dams_num", "perc_complete",
                       "data_gap",
                       "aridity", "runoff_ratio"]]
df_final.to_csv(results_path + 'camels_de_ROBIN.csv', index=False)
print("Finished saving data.")

# plot standard Budyko plot
fig = plt.figure(figsize=(4, 3), constrained_layout=True)
axes = plt.axes()
im = axes.scatter(df_final["aridity"], 1 - df_final["runoff_ratio"], s=10, c="tab:blue", alpha=0.8, lw=0)
axes.set_xlabel("Aridity [-]")
axes.set_ylabel("1 - Runoff ratio [-]")
axes.set_xlim([0, 2])
axes.set_ylim([-0.25, 1.25])
helper_fcts.plot_Budyko_limits(df_final["aridity"], 1 - df_final["runoff_ratio"], axes)
helper_fcts.plot_Budyko_curve(np.linspace(0, 10, 100), axes)
fig.savefig(figures_path + "Budyko_plot_ROBIN_CAMELS_DE.png", dpi=600, bbox_inches='tight')
plt.close()

# plot map
fig, ax = plt.subplots(figsize=(8, 4))
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
scatter_df = gpd.GeoDataFrame(df_final, geometry=gpd.points_from_xy(df_final.Longitude, df_final.Latitude))
world.boundary.plot(ax=ax, linewidth=0.5, color='black')
world.plot(ax=ax, color='lightgrey', edgecolor='black', )
scatter_df.plot(ax=ax, markersize=10, color='tab:blue')
plt.title('CAMELS-DE catchments')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim(4, 16)
plt.ylim(46, 56)
plt.gca().set_aspect('equal', adjustable='box')
fig.savefig(figures_path + "map_ROBIN_CAMELS_DE" + ".png", dpi=600, bbox_inches='tight')
plt.close()

### alternative selection of all original ROBIN files

# original ROBIN list
df_ROBIN_orig = df[df["ROBIN_flag"] == 1]

df_ROBIN_orig["ID"] = df_ROBIN_orig["gauge_id"]
df_ROBIN_orig["Name"] = df_ROBIN_orig["gauge_name"]
df_ROBIN_orig["Longitude"] = df_ROBIN_orig["gauge_lon"]
df_ROBIN_orig["Latitude"] = df_ROBIN_orig["gauge_lat"]
df_ROBIN_orig["Catchment Area (km²)"] = df_ROBIN_orig["area"]
df_ROBIN_orig["Measuring Organisation"] = df_ROBIN_orig["federal_state"]
df_ROBIN_orig["Level 1 or Level 2 (Catchment Development)"] = df_ROBIN_orig["humanimpact_flag"]
df_ROBIN_orig["Level 1 or Level 2 (Data Quality)"] = df_ROBIN_orig["data_flag"]
df_ROBIN_orig["Record Length"] = np.round(df_ROBIN_orig["record_length"], 2)
df_ROBIN_orig["Missing Data Criteria Met? (Yes/No)"] = df_ROBIN_orig["data_gap"].apply(
    lambda x: "Yes" if x < pd.Timedelta(days=1095) else "No")

# count how many in df_final Level 1 or Level 2 are Level 1 or Level 2
print(df_ROBIN_orig["Level 1 or Level 2 (Catchment Development)"].value_counts())
print(df_ROBIN_orig["Level 1 or Level 2 (Data Quality)"].value_counts())

df_ROBIN_final = df_ROBIN_orig[["ID", "Name", "Longitude", "Latitude", "Catchment Area (km²)", "Measuring Organisation",
                                "Level 1 or Level 2 (Catchment Development)", "Level 1 or Level 2 (Data Quality)",
                                "Record Length", "Missing Data Criteria Met? (Yes/No)",
                                "gauge_name", "water_body_name", "artificial_surfaces_perc", "dams_num",
                                "perc_complete", "data_gap", "area_error",
                                "aridity", "runoff_ratio", "data_flag"]]
df_ROBIN_final.to_csv(results_path + 'camels_de_ROBIN_orig.csv', index=False)
print("Finished saving data.")

# plot standard Budyko plot
fig = plt.figure(figsize=(4, 3), constrained_layout=True)
axes = plt.axes()
im = axes.scatter(df_ROBIN_final["aridity"], 1 - df_ROBIN_final["runoff_ratio"], s=10, c="tab:blue", alpha=0.8, lw=0)
axes.set_xlabel("Aridity [-]")
axes.set_ylabel("1 - Runoff ratio [-]")
axes.set_xlim([0, 2])
axes.set_ylim([-0.25, 1.25])
helper_fcts.plot_Budyko_limits(df_ROBIN_final["aridity"], 1 - df_ROBIN_final["runoff_ratio"], axes)
helper_fcts.plot_Budyko_curve(np.linspace(0, 10, 100), axes)
fig.savefig(figures_path + "Budyko_plot_ROBIN_CAMELS_DE_orig.png", dpi=600, bbox_inches='tight')
plt.close()

# plot map
fig, ax = plt.subplots(figsize=(8, 4))
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
scatter_df = gpd.GeoDataFrame(df_ROBIN_final, geometry=gpd.points_from_xy(df_ROBIN_final.Longitude, df_ROBIN_final.Latitude))
world.boundary.plot(ax=ax, linewidth=0.5, color='black')
world.plot(ax=ax, color='lightgrey', edgecolor='black', )
scatter_df.plot(ax=ax, markersize=10, color='tab:blue')
plt.title('CAMELS-DE catchments')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim(4, 16)
plt.ylim(46, 56)
plt.gca().set_aspect('equal', adjustable='box')
fig.savefig(figures_path + "map_ROBIN_CAMELS_DE_orig" + ".png", dpi=600, bbox_inches='tight')
plt.close()

# check reference network (https://www.umweltbundesamt.de/en/monitoring-on-das/cluster/water-supply/ww-i-3/indicator)
df_UBA = pd.read_csv(results_path + "reference_gauges_Germany_UBA.csv", sep=',', skiprows=0, encoding='latin-1')
df_ROBIN_final = pd.read_csv(results_path + 'camels_de_ROBIN_orig.csv', sep=',', skiprows=0, encoding='latin-1')
df_UBA["in_CAMELS_DE"] = df_UBA["Pegelname"].isin(df_attr["gauge_name"])
df_UBA["in_ROBIN_DE"] = df_UBA["Pegelname"].isin(df_ROBIN_final["gauge_name"])
print(sum(df_UBA["in_CAMELS_DE"]))
print(sum(df_UBA["in_ROBIN_DE"]))
# some problems due to Umlaute, but overall it seems like there is not much overlap
