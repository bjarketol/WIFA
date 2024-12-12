import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# Load the datasets
dat_10min = xr.load_dataset("A1_profiling_lidar_10min.nc")
dat_10min["time"] = dat_10min.time - np.timedelta64(5, "m")
tempdat = xr.load_dataset("A_1-d03_2023-08-23.nc")
dat_1s = xr.load_dataset("A1_profiling_lidar_1s.nc")

# Interpolate temperature data to match the height of 'dat_10min'
tempdat_interp = tempdat.interp(height=dat_10min.height, method="nearest")

# Calculate the Turbulence Intensity (TI) from 'dat_1s'
u = np.sqrt(dat_1s["u"] ** 2 + dat_1s["v"] ** 2)
u_resampled = u.resample(time="10T").mean()  # Resample to 10-minute intervals
u_std = u.resample(
    time="10T"
).std()  # Calculate the standard deviation over the intervals

# Calculate TI as the ratio of standard deviation to mean, keeping the time and height dimensions
TI = u_std / u_resampled

# Extract relevant data variables from the 10-minute dataset
wind_speed = dat_10min["U"]
wind_direction = dat_10min["WD"]
turbulence_intensity = dat_10min["w"]
potential_temperature = tempdat_interp["TH"]

# Create an xarray Dataset to store all these variables
ds_combined = xr.Dataset(
    {
        "wind_speed": wind_speed,
        "wind_direction": wind_direction,
        "turbulence_intensity": turbulence_intensity,
        "potential_temperature": potential_temperature,
        "turbulence_intensity": TI,
    }
)

# Optionally, you can add metadata to each variable to describe it
ds_combined["wind_speed"].attrs["units"] = "m/s"
ds_combined["wind_direction"].attrs["units"] = "degrees"
ds_combined["turbulence_intensity"].attrs["units"] = "ratio"
ds_combined["potential_temperature"].attrs["units"] = "K"
ds_combined["turbulence_intensity"].attrs["units"] = "ratio"


# Plot function for visualization
def plot_variable(data, title, units, filename):
    plt.figure(figsize=(10, 6))

    # Create a colormap that assigns NaN values to black
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("black")

    # Convert DataArray to a masked array that handles NaNs
    data_array = data.to_masked_array()

    # Plot the data with pcolormesh, ensuring matching dimensions
    time, height = np.meshgrid(data.time, data.height)
    pcm = plt.pcolormesh(time, height, data_array.T, shading="auto", cmap=cmap)

    # Add color bar and labels
    plt.colorbar(pcm, label=f"{title} ({units})")
    plt.title(f"{title} vs Height and Time")
    plt.xlabel("Time")
    plt.ylabel("Height (m)")

    # Save the plot to a file
    plt.savefig(f"{filename}.png")
    plt.close()


# Plot all variables from the dataset
plot_variable(ds_combined["wind_speed"], "Wind Speed", "m/s", "wind_speed_plot")
plot_variable(
    ds_combined["wind_direction"], "Wind Direction", "degrees", "wind_direction_plot"
)
plot_variable(
    ds_combined["turbulence_intensity"],
    "Turbulence Intensity",
    "ratio",
    "turbulence_intensity_plot",
)
plot_variable(
    ds_combined["potential_temperature"],
    "Potential Temperature",
    "K",
    "potential_temperature_plot",
)
plot_variable(
    ds_combined["turbulence_intensity"], "Turbulence Intensity (TI)", "ratio", "TI_plot"
)

ds_combined = ds_combined.sel(time=slice("2023-08-24", "2023-08-24T23:59:59"))
ds_combined["time"] = ds_combined.time.astype(str)

# Print out the dataset to verify the structure
print(ds_combined)


# Save the combined dataset to a new NetCDF file
ds_combined.to_netcdf("combined_signals.nc")
