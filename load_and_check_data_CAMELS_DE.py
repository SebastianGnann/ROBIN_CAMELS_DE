import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# from mpl_toolkits.basemap import Basemap

# prepare data
data_path = "C:/Users/Sebastian/Documents/Data/camels_de/"

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
df_topo = pd.read_csv(data_path + "CAMELS_DE_topographic_attributes.csv",
                      sep=',', skiprows=0, encoding='latin-1')
df_climate = pd.read_csv(data_path + "CAMELS_DE_climatic_attributes.csv",
                         sep=',', skiprows=0, encoding='latin-1')
df_landcover = pd.read_csv(data_path + "CAMELS_DE_landcover_attributes.csv",
                           sep=',', skiprows=0, encoding='latin-1')
df_humaninfluence = pd.read_csv(data_path + "CAMELS_DE_humaninfluence_attributes.csv",
                                sep=',', skiprows=0, encoding='latin-1')

df_attr = pd.merge(df_topo, df_climate, on='gauge_id')
df_attr = pd.merge(df_attr, df_landcover, on='gauge_id')
df_attr = pd.merge(df_attr, df_humaninfluence, on='gauge_id')

### checks and additional attributes
# record length
# missing value fraction
# strange Q values: nan, negative, many 0s, consecutive days with same flow, step changes
# plot map
# plot timeseries

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

# NOTE: information on individual gauging stations, rating curve uncertainties, etc.
# not possible without considerable efforts, partly because it is different for every state
# flags: level 1 = 1, level 2 = 2, not in network = 0
# todo: add data source in metadata: CAMELS-DE or more specific

gauge_id_list = []
mean_P_list = []  # precip average
mean_Q_list = []  # Q average
perc_complete_list = []
record_length_list = []
data_gap_list = []
data_flag_list = []

for id in df_attr["gauge_id"]:
    print(id)

    # load data
    df_tmp = pd.read_csv(data_path + "timeseries/CAMELS_DE_hydromet_timeseries_" + str(id) + ".csv", sep=',')
    df_tmp["date"] = pd.to_datetime(df_tmp["date"])

    # check data
    perc_complete = np.sum(~np.isnan(df_tmp["discharge_spec_obs"].values)) / len(df_tmp["discharge_spec_obs"].values)
    record_length = len(df_tmp["date"])
    data_gaps = 0  # check longest gap (<3y)

    # add data
    gauge_id_list.append(id)
    perc_complete_list.append(perc_complete)
    record_length_list.append(record_length)
    data_gap_list.append(data_gaps)
    mean_P_list.append(np.nanmean(df_tmp["precipitation_mean"]))
    mean_Q_list.append(np.nanmean(df_tmp["discharge_spec_obs"]))

    # plot data
    fig, ax = plt.subplots(figsize=(12, 4), tight_layout=True)
    im = ax.plot(df_tmp["date"], df_tmp["discharge_spec_obs"], alpha=0.9)
    ax.set_ylabel("Streamflow (mm/d)")
    # ax.set_xlim([0., 5.])
    # ax.set_ylim([0., 5.])
    # log y axis
    fig.savefig(figures_timeseries_path + "CAMELS_DE_" + id + ".png", dpi=600, bbox_inches='tight')
    plt.close()

# create dataframe
df = pd.DataFrame()
df["gauge_id"] = gauge_id_list  # to check
df = pd.merge(df_attr, df, on='gauge_id')
df["mean_P"] = mean_P_list
df["mean_Q"] = mean_Q_list
df["perc_complete"] = perc_complete_list
df["record_length"] = record_length_list
df["data_gap"] = data_gap_list

# human impacts
# artificial_surfaces_perc dams_num
if df["artificial_surfaces_perc"][id] < 10:
    urban_flag = 1
elif df["artificial_surfaces_perc"][id] < 20:
    urban_flag = 2
else:
    urban_flag = 0

# no dams etc.

# todo: landuse change?

# data quality
if df["record_length"][id] < 20:
    length_flag = 0

# area check
if np.abs(df["area_metadata"] - df["area"] / df["area_metadata"]) > 0.1:
    area_flag = 0

# extract attributes that are necessary
# gauge_id provider_id gauge_name water_body_name gauge_lon gauge_lat area_metadata area
#

# save results
df.to_csv(results_path + 'camels_de_ROBIN.csv', index=False)
print("Finished saving data.")

# scatter
fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)
line_x = np.linspace(-10, 10, 100)
line_y = line_x
ax.plot(line_x, line_y, 'k--')
im = ax.scatter(df["P_mean"], df["Q_mean"], s=20, alpha=0.9)
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xlim([0., 5.])
ax.set_ylim([0., 5.])
plt.show()

# map
fig, ax = plt.subplots(figsize=(4, 4))
m = Basemap(projection='robin', resolution='l', area_thresh=1000.0, lat_0=0, lon_0=0)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='lightgrey', lake_color='white')
m.drawmapboundary(fill_color='white')
x, y = m(df["gauge_lon"].values, df["gauge_lat"].values)
scatter = m.scatter(x, y, s=20, c=df["mean_Q"], alpha=0.9, vmin=0.0, vmax=5.0, cmap='viridis')  # invert colormap
cbar = plt.colorbar(scatter, ax=ax, pad=0.02, shrink=0.3, aspect=20)
ax.set_xlim(np.min(x) * 0.99, np.max(x) * 1.01)
ax.set_ylim(np.min(y) * 0.99, np.max(y) * 1.01)
cbar.set_label('Q mean [-]', rotation=270, labelpad=15)
plt.tight_layout()
plt.show()
