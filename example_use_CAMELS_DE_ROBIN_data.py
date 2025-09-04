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

# load attributes
df_topo = pd.read_csv(data_path + "CAMELS_DE_topographic_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
df_climate = pd.read_csv(data_path + "CAMELS_DE_climatic_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
df_landcover = pd.read_csv(data_path + "CAMELS_DE_landcover_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
df_humaninfluence = pd.read_csv(data_path + "CAMELS_DE_humaninfluence_attributes.csv", sep=',', skiprows=0, encoding='latin-1')
# there are more attributes available...

df_attr = pd.merge(df_topo, df_climate, on='gauge_id')
df_attr = pd.merge(df_attr, df_landcover, on='gauge_id')
df_attr = pd.merge(df_attr, df_humaninfluence, on='gauge_id')

# ROBIN catchments
df_ROBIN_list = pd.read_csv("D:/Python/ROBIN_CAMELS_DE/results/camels_de_ROBIN.csv")
df_CAMELS_ROBIN = df_attr[df_attr["gauge_id"].isin(df_ROBIN_list["ID"].values)]
df_CAMELS_ROBIN = df_CAMELS_ROBIN.reset_index(drop=True)

# plot map
fig, ax = plt.subplots(figsize=(8, 4))
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
scatter_df = gpd.GeoDataFrame(df_attr, geometry=gpd.points_from_xy(df_attr.gauge_lon, df_attr.gauge_lat))
scatter_df_ROBIN = gpd.GeoDataFrame(df_CAMELS_ROBIN, geometry=gpd.points_from_xy(df_CAMELS_ROBIN.gauge_lon, df_CAMELS_ROBIN.gauge_lat))
world.boundary.plot(ax=ax, linewidth=0.5, color='black')
world.plot(ax=ax, color='lightgrey', edgecolor='black', )
scatter_df.plot(ax=ax, markersize=10, color='grey')
scatter_df_ROBIN.plot(ax=ax, markersize=10, color='tab:blue')
plt.title('CAMELS-DE catchments')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim(4, 16)
plt.ylim(46, 56)
plt.gca().set_aspect('equal', adjustable='box')
fig.savefig(figures_path + "map_ROBIN_CAMELS_DE_comparison" + ".png", dpi=600, bbox_inches='tight')
plt.show()
#plt.close()

# plot map with attribute
fig, ax = plt.subplots(figsize=(8, 4))
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
scatter_df_ROBIN = gpd.GeoDataFrame(df_CAMELS_ROBIN, geometry=gpd.points_from_xy(df_CAMELS_ROBIN.gauge_lon, df_CAMELS_ROBIN.gauge_lat))
world.boundary.plot(ax=ax, linewidth=0.5, color='black')
world.plot(ax=ax, color='lightgrey', edgecolor='black', )
scatter_df_ROBIN.plot(ax=ax, markersize=10, column=scatter_df_ROBIN['elev_mean'], cmap='viridis', vmin=0, vmax=1000)
plt.title('ROBIN catchments')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim(4, 16)
plt.ylim(46, 56)
plt.gca().set_aspect('equal', adjustable='box')
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1000)), ax=ax, pad=0.02, shrink=0.5, aspect=20)
cbar.set_label('Elevation [m]', rotation=270, labelpad=15)
fig.savefig(figures_path + "map_ROBIN_CAMELS_DE_elevation" + ".png", dpi=600, bbox_inches='tight')
plt.show()
#plt.close()

# load and plot individual catchment
id = 'DE111270'
df_obs = pd.read_csv(data_path + "timeseries/CAMELS_DE_hydromet_timeseries_" + str(id) + ".csv", sep=',')
df_obs["date"] = pd.to_datetime(df_obs["date"])
df_sim = pd.read_csv(data_path + "timeseries_simulated/CAMELS_DE_discharge_sim_" + str(id) + ".csv", sep=',')
df_sim["date"] = pd.to_datetime(df_sim["date"])
df_obs = pd.merge(df_obs, df_sim[["date", "pet_hargreaves"]], on='date')
# remove NaNs at beginning
first_valid_index = df_obs["discharge_spec_obs"].first_valid_index()
df_obs = df_obs.loc[first_valid_index:]

# plot time series
fig, ax = plt.subplots(figsize=(8, 3))
ax2 = ax.twinx()
ax.plot(df_obs["date"], df_obs["discharge_spec_obs"], color='tab:blue', label='Observed Discharge', lw=1.5)
ax2.plot(df_obs["date"], df_obs["precipitation_mean"], color='grey', label='Precipitation', lw=1)
ax2.invert_yaxis()
#ax.set_yscale('log')
ax.set_ylabel('Discharge [mm/d]')
ax2.set_ylabel('Precipitation [mm/d]')
plt.title(f'CAMELS-DE catchment {id}')
plt.tight_layout()
ax.set_xlim(pd.to_datetime('2001-01-11'), pd.to_datetime('2020-10-31'))
fig.savefig(figures_timeseries_path + "timeseries_CAMELS_DE_catchment_" + str(id) + ".png", dpi=600, bbox_inches='tight')
plt.show()
#plt.close()

# check reference network (https://www.umweltbundesamt.de/en/monitoring-on-das/cluster/water-supply/ww-i-3/indicator)
df_UBA = pd.read_csv(results_path + "reference_gauges_Germany_UBA.csv", sep=',', skiprows=0, encoding='latin-1')
df_ROBIN_final = pd.read_csv(results_path + 'camels_de_ROBIN_orig.csv', sep=',', skiprows=0, encoding='latin-1')
df_UBA["in_CAMELS_DE"] = df_UBA["Pegelname"].isin(df_attr["gauge_name"])
df_UBA["in_ROBIN_DE"] = df_UBA["Pegelname"].isin(df_ROBIN_final["gauge_name"])
print(sum(df_UBA["in_CAMELS_DE"]))
print(sum(df_UBA["in_ROBIN_DE"]))
# some problems due to Umlaute, but overall it seems like there is not much overlap

