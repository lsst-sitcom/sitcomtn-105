---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: LSST
  language: python
  name: lsst
---

# Compare TMA balancing data from previous events with the next events

```{abstract}
We are getting ready to balance the telescope twice in the next weeks. First, we will balance the telescope with ComCam and M2 Glass. The M2 Glass and M2 Surrogate have similar weights, with a small difference. We expect the torques applied by the elevation drives will be very close to the previous balancing event(s). A couple of weeks later, we will repeat the procedure with ComCam, M2 Glass, and M1M3 Glass. The M1M3 Glass and M1M3 Cell assembly is much heavier than the M1M3 Mass Simulator (yellow cross) and hundreds of kilograms heavier than the M1M3 Surrogate and M1M3 Cell configuration. This procedure will be much more delicate due to the size and mass of the mirror.

We want to establish a baseline before we start the procedure, and we need someone to review the data to determine whether we can proceed quickly.

The links below point to old night logs that might contain useful information. Feel free to unlink them if they are not useful.

Here is an approximate timeline of different integration phases where we needed to re-balance the telescope. We do not necessarily need the whole process. We need the torques once the telescope is already balanced as a baseline.

May to Aug 2023 - M1M3 Surrogate and M1M3 Cell on the TMA

Nov 2023 to Jan 2024 - M1M3 Surrogate and Cell, M2 Surrogate and Cell on the TMA

Feb to Apr 2024 - M2 Surrogate and Cell on the TMA
```

+++

```{code-cell} ipython3
:tags: [hide-cell]
:mystnb:
:  code_prompt_show: "My show prompt for {type}"
:  code_prompt_hide: "My hide prompt for {type}"
```

```{code-cell} ipython3
# Notebook extensions for formatting and auto-reload libraries
%matplotlib inline
%load_ext lab_black
%load_ext autoreload
%autoreload 2

# Standard Python Libraries
import os
import sys

from astropy.time import Time

from lsst_efd_client import EfdClient

module_path = os.path.abspath(os.path.join("."))
sys.path.append(module_path + "/python")

# Import utility functions from the python directory
from utils import *
```

## Identify time periods with telemetry and where the TMA is still

In order to do that we will check the "lsst.sal.MTMount.logevent_elevationMotionState" and "lsst.sal.MTMount.logevent_azimuthMotionState" topics and select time ranges where
the `state` is equal to 1, meaning that the TMA is `stopped`

Once we have the time ranges in both azimuth and elevation we select the overlaps between both sets.

```{code-cell} ipython3
# Create a directory to save plots
plot_dir = "./plots"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
```

## Dates with interesting data

The most interesting data for this analysis are those from BLOCK-177 https://rubinobs.atlassian.net/browse/BLOCK-177
initially designed for laser tracker tests.
The TMA is still on the azimuth axis and the elevation increase and decrease by steps of 5 (or 10 ?) degrees

May to Aug 2023 - M1M3 Surrogate and M1M3 Cell on the TMA

* 2023-06-22

Nov 2023 to Jan 2024 - M1M3 Surrogate and Cell, M2 Surrogate and Cell on the TMA

* 2024-01-06 - Az: 59.6 deg - Block 177 between 09:00 and 10:00
* 2024-01-09 - Az: 59.6 deg
* 2024-01-10 - Az: 60.3 deg
* 2024-01-12 - Az: 0 deg

Feb to Apr 2024 - M2 Surrogate and Cell on the TMA

 * 2024-03-28

Coarse balance with Yellow cross - M2+Cell and ComCam on TMA 

* 2024-09-03 14:30 - 16:00

```{code-cell} ipython3
# Define the time period that we are going to investigate

# date_dict = {
#    "date": "2024-01-06",
#    "start_time": "09:00:00.00",
#    "end_time": "10:00:00.00",
# }
date_dict = {
    "date": "2024-09-03",
    "start_time": "14:30:00.00",
    "end_time": "16:00:00.00",
}
# date_dict = {
#    "date": "2024-01-03",
#    "start_time": "08:30:00.00",
#    "end_time": "09:30:00.00",
# }

date = date_dict["date"]
start_time = Time(f"{date} {date_dict['start_time']}")
end_time = Time(f"{date} {date_dict['end_time']}")

# Select the EFD server and initialize EFD client
# summit_efd contains the most recent data
# usdf_efd contains older / archived data

# client = EfdClient("summit_efd")
client = EfdClient("usdf_efd")
```

```{code-cell} ipython3
# Get all time ranges where the TMA is still in azimuth and in elevation
t_range_azi = await get_time_range_by_axis(client, "azimuth", start_time, end_time)
t_range_ele = await get_time_range_by_axis(client, "elevation", start_time, end_time)

# Find the overlaps between the 2 sets of time ranges
# We want data spanning over a minimum amount of time during the overlap period
# min_delta is in seconds
min_delta = 20
overlaps = get_overlaps(t_range_azi, t_range_ele, min_delta)
print(
    f"We found {len(overlaps)} time periods where the TMA is still in both azimuth and elevation"
)
```

# Torque versus Elevation Angle

```{code-cell} ipython3
# We are going to plot the Torque as a function of the elevation angle for a given position of the TMA in azimuth
# So we will first identify what is the most common TMA position in azimuth over all the overlap periods and we will make the analysis for this
# position only
az_min, az_max = get_common_azimuth(overlaps, client)

# Create pandas dataframe with the relevant EFD data
df = make_dataframe(overlaps, client, az_min, az_max)

# Plot torque versus elevation angle
plot_torque_versus_elevation(df, date, 0.5 * (az_max + az_min), plot_dir)
```

## Plot data in chronological order

It is useful to plot a few quantities in chronological order in order to be able to understand the exact sequence of events

```{code-cell} ipython3
# We will retrieve data back t_back seconds in the past
t_back = 1000
plot_history(
    df, client, t_back, date, 0.5 * (az_max + az_min), plot_dir, save_plot=True
)
```

```{code-cell} ipython3

```
