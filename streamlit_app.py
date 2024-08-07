import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.constants import c, k, h, pi, parsec
import matplotlib.colors as mcolors
from scipy.special import j1, jn_zeros


st.markdown(
    """
    # Stars available for Intensity Interferometry

    On this webpage one can explore and plot the HR-diagram of available stars that could be observed using intensity interferometry.
    """
    )


# Initialize a dictionary to store data from each column
data = {}
plt.style.use('dark_background')

# Open the CSV file
with open('10000stars_data.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)  # Read the header row

    # Create an empty list for each column
    for column_name in header:
        data[column_name] = []

    # Iterate over each row in the CSV file
    for row in reader:
        # Iterate over each item in the row and append it to the respective column list
        for idx, item in enumerate(row):
            column_name = header[idx]  # Get the corresponding column name
            data[column_name].append(item)

# Create colormap that corresponds to temperatures of the stars (color table taken from http://www.vendian.org/mncharity/dir3/blackbody/)
def parse_colormap(file_path):
    temperatures = []
    colors = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            temperature = int(parts[0])
            color = parts[2]
            temperatures.append(temperature)
            colors.append(color)
    return temperatures, colors


def create_custom_colormap(temperatures, colors):
    norm = mcolors.Normalize(vmin=min(temperatures), vmax=max(temperatures))
    tuples = list(zip(map(norm, temperatures), colors))
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", tuples)
    return cmap


temperatures, bb_colors = parse_colormap('blackbody_colors')
bb_cmap = create_custom_colormap(temperatures, bb_colors)

norm = plt.Normalize(vmin=0, vmax=29800)


# function converting the entries retrieved from the .csv from string to float
def convert_strings_to_floats(input_array):
    output_array = []
    for element in input_array:
        try:
            converted_float = float(element)
        except ValueError:
            converted_float = np.nan  # Use NaN for invalid or empty strings
        output_array.append(converted_float)
    return output_array

# Extracting entries and converting

Phi_V = np.array(data['Phi_V'])
Phi_V = convert_strings_to_floats(Phi_V)
Phi_V = np.array(Phi_V)

diameter_V = np.array(data['Diameter_V'])
diameter_V = convert_strings_to_floats(diameter_V)
diameter_V = np.array(diameter_V)

temps = np.array(data['Temp'])
temps = convert_strings_to_floats(temps)
temps = np.array(temps)

dist = convert_strings_to_floats(data['Distance'])
dist = np.array(dist)

Vmag = convert_strings_to_floats(data['Vmag'])
Vmag = np.array(Vmag)

Dec_decimal = convert_strings_to_floats(data['Dec_decimal'])
Dec_decimal = np.array(Dec_decimal)

RA_decimal = convert_strings_to_floats(data['RA_decimal'])
Ra_decimal = np.array(RA_decimal)

sii_analyzed = np.array(data['SII'])
sii_analyzed = [True if entry.lower() == 'x' else False for entry in sii_analyzed]

# Remove entries where distance is NaN
valid_indices = ~np.isnan(dist) & (dist > 0)

dist = dist[valid_indices]
Phi_V = Phi_V[valid_indices]
diameter_V = diameter_V[valid_indices]
temps = temps[valid_indices]
Vmag = Vmag[valid_indices]
sii_analyzed = np.array(sii_analyzed)[valid_indices]
Dec_decimal = np.array(Dec_decimal)[valid_indices]
RA_decimal = np.array(RA_decimal)[valid_indices]


def visibility(b, theta, lambda_=540 * 1e-9):
    """The squared visibility, often denoted in papers as |V_12|^2 and equals g**(2)-1"""
    input = pi * b * theta / lambda_
    if b == 0:
        I = 1
    else:
        I = (2 * j1(input) / input) ** 2
    return I

def mas_to_ang(theta_mas):
    return theta_mas/ 1000 * pi / (3600 * 180)

def baseline_needed(theta, lambda_=540e-9):
    """Function to calculate the minimum baseline needed for a given theta."""
    theta = mas_to_ang(theta)
    j1_root = jn_zeros(1, 1)
    return float(j1_root/(pi*theta / (540 * 1e-9)))

# Determine baselines needed for each star
baselines_needed = []

for i in range(len(diameter_V)):
    baselines_needed.append(baseline_needed(diameter_V[i]))


baselines_needed_array = np.array(baselines_needed)


# determine the inverse of the diameter and the diameter in angular degrees
inverse_diameter = []
diameter_in_rad = []

for i in range(len(diameter_V)):
    diameter_in_rad.append(mas_to_ang(diameter_V[i]))
    inverse_diameter.append(1 / diameter_V[i])
inverse_diameter = np.array(inverse_diameter)
print(dist[7])

# Determining Phi from the magnitude of the star
def mag_from_phi(Phi, wavelength=540 * 1e-9):
    nu = c / wavelength
    magnitude = -2.5 * (22.44 + np.log10(2 * nu * h * Phi))
    return magnitude

# Plot using the colormap based on temperatures
fig, ax1 = plt.subplots()

# Plot the first dataset with the colormap based on temperatures
sc = ax1.scatter(inverse_diameter, Phi_V, c=temps, cmap=bb_cmap, marker='.', norm=norm)
sc2 = ax1.scatter(inverse_diameter[sii_analyzed], Phi_V[sii_analyzed], c=temps[sii_analyzed], cmap=bb_cmap, marker='*', label='SII Analyzed Stars', norm=norm)
cbar = plt.colorbar(sc, ax=ax1, label='Temperature (K)', pad=0.15)
ax1.set_yscale('log')
ax1.set_xlabel('1/θ [mas$^{-1}$]')
ax1.set_ylabel(r'Φ [photons m$^{-2}$ s$^{-1}$ Hz$^{-1}$]')
ax1.set_title('Φ vs θ')
#ax1.set_xlim(0, 25)

ax2 = ax1.twinx()
ax3 = ax1.twiny()

# Set the limits for the second y-axis based on the transformation
phi_min, phi_max = ax1.get_ylim()
thetainverse_min, thetainverse_max = ax1.get_xlim()
if thetainverse_min < 0:
    thetainverse_min = 0

ax2.set_ylim(mag_from_phi(phi_min), mag_from_phi(phi_max))
ax2.set_ylabel('magnitude')

ax3.scatter(baselines_needed, Phi_V, c=temps, cmap=bb_cmap, marker='')
ax3.set_xlabel('baseline needed [m]')

ax2.grid(True)
ax3.grid(True)

st.markdown(
    """
    The stars used are from the Yale Bright Star Catalog which contains 9110 of the brightest stars (found at http://tdc-www.harvard.edu/catalogs/bsc5.html). 
    The cataloge is in ASCII format and was converted to a .JSON format in the repository https://github.com/brettonw/YaleBrightStarCatalog. 
    The data file used from that repository is bsc5-all.json.
    """
    )

st.write("Of the 9110 stars available, ",len(Vmag)," had enough information to create suitable HR-diagrams.")
st.write("The ones marked with a '★' had their diameters already measured by Hanbury Brown.")

st.pyplot(plt)


st.markdown(
    """
    ## Input

    Here you can select the available baseline, the minimum magnitude of the stars you want to observe and if desired the RA and Dec. Alternatively you can specify the angular diameter instead of the baseline.
    
    ### Filtering Options

    Use the options below to filter the stars based on your requirements.

    #### Baseline or Angular Diameter
    Choose whether to filter by baseline or angular diameter. If you select "baseline", you will specify the range of baseline lengths in meters. If you select "angular diameter", you will specify the range of angular diameters in milliarcseconds.
    """
)

filter_b_d = st.radio(
    "Filter by baseline or angular diameter",
    ["baseline", "angular diameter"]
)
print(filter_b_d)

if filter_b_d == "baseline":
    st.markdown("#### Baseline Selection")
    baseline_available = st.slider("Baseline available in meters", 0, 3500, (0, 3500))
    baseline_min, baseline_max = baseline_available
    indices_baseline = np.argwhere((baselines_needed_array >= baseline_min) & (baselines_needed_array <= baseline_max))
    baselines_needed_ = baselines_needed_array[indices_baseline].reshape(-1)
    Vmag_ = Vmag[indices_baseline].reshape(-1)
    Phi_V_ = Phi_V[indices_baseline].reshape(-1)
    diameter_V_ = diameter_V[indices_baseline].reshape(-1)
    temps_ = temps[indices_baseline].reshape(-1)
    inverse_diameter_ = inverse_diameter[indices_baseline].reshape(-1)
    dist_ = dist[indices_baseline].reshape(-1)
    sii_analyzed_ = sii_analyzed[indices_baseline].reshape(-1)
    Dec_decimal_ = Dec_decimal[indices_baseline].reshape(-1)
    RA_decimal_ = RA_decimal[indices_baseline].reshape(-1)

else:
    st.markdown("#### Angular Diameter Selection")
    desired_angular_diameter = st.slider("Angular diameter of star desired in milliarcseconds", 0.00, 50.00, (0.00, 50.00))
    diameter_min, diameter_max = desired_angular_diameter
    indices_angular_diameter = np.argwhere((diameter_V >= diameter_min) & (diameter_V <= diameter_max))
    baselines_needed_ = baselines_needed_array[indices_angular_diameter].reshape(-1)
    Vmag_ = Vmag[indices_angular_diameter].reshape(-1)
    Phi_V_ = Phi_V[indices_angular_diameter].reshape(-1)
    diameter_V_ = diameter_V[indices_angular_diameter].reshape(-1)
    temps_ = temps[indices_angular_diameter].reshape(-1)
    inverse_diameter_ = inverse_diameter[indices_angular_diameter].reshape(-1)
    dist_ = dist[indices_angular_diameter].reshape(-1)
    sii_analyzed_ = sii_analyzed[indices_angular_diameter].reshape(-1)
    Dec_decimal_ = Dec_decimal[indices_angular_diameter].reshape(-1)
    RA_decimal_ = RA_decimal[indices_angular_diameter].reshape(-1)

st.markdown(
    """
    #### Minimum Magnitude

    Set the minimum magnitude of the star you want to observe. This helps filter out stars that are too faint.
    """
)
magnitude_min = st.number_input("Minimum magnitude of star", -2, 8, value=8)
indices_magnitude = np.argwhere(Vmag_ < magnitude_min)

baselines_needed__ = baselines_needed_[indices_magnitude].reshape(-1)
Phi_V__= Phi_V_[indices_magnitude].reshape(-1)
diameter_V__ = diameter_V_[indices_magnitude].reshape(-1)
temps__ = temps_[indices_magnitude].reshape(-1)
inverse_diameter__ = inverse_diameter_[indices_magnitude].reshape(-1)
dist__ = dist_[indices_magnitude].reshape(-1)
Vmag__ = Vmag_[indices_magnitude].reshape(-1)
sii_analyzed__ = sii_analyzed_[indices_magnitude].reshape(-1)
Dec_decimal__ = Dec_decimal_[indices_magnitude].reshape(-1)
RA_decimal__ = RA_decimal_[indices_magnitude].reshape(-1)

st.markdown(
    """
    #### Declination Range

    Toggle the option below if you want to specify a range for declination.
    """
)
on_dec = st.toggle("Specify Declination")

if on_dec:
    declination_range = st.slider("Declination range in degrees", -90, 90, (-90, 90))
else:
    declination_range = [-90, 90]

declination_min, declination_max = declination_range
indices_declination = np.argwhere((Dec_decimal__ >= declination_min) & (Dec_decimal__ <= declination_max))
baselines_needed___ = baselines_needed__[indices_declination].reshape(-1)
Vmag___ = Vmag__[indices_declination].reshape(-1)
Phi_V___ = Phi_V__[indices_declination].reshape(-1)
diameter_V___ = diameter_V__[indices_declination].reshape(-1)
temps___ = temps__[indices_declination].reshape(-1)
inverse_diameter___ = inverse_diameter__[indices_declination].reshape(-1)
dist___ = dist__[indices_declination].reshape(-1)
sii_analyzed___ = sii_analyzed__[indices_declination].reshape(-1)
Dec_decimal___ = Dec_decimal__[indices_declination].reshape(-1)
RA_decimal___ = RA_decimal__[indices_declination].reshape(-1)

st.markdown(
    """
    #### Right Ascension Range

    Toggle the option below if you want to specify a range for right ascension.
    """
)
on_ra = st.toggle("Specify Right Ascension Range")

if on_ra:
    # Define the RA range in hours
    ra_min, ra_max = 0, 24
    
    # RA slider
    ra_start = st.number_input("Minimum Right Ascension in hours", 0, 24, value=0)
    ra_end = st.number_input("Maximum Right Ascension in hours", 0, 24, value=24)
    
    st.write(f"RA range: {ra_start}h to {ra_end}h")
    indices_ra = np.argwhere((RA_decimal___ >= ra_start) & (RA_decimal___ < ra_end))

else:
    ra_start = 0
    ra_end = 24
    indices_ra = np.argwhere((RA_decimal___ >= ra_start) & (RA_decimal___ <= ra_end))

baselines_needed____ = baselines_needed___[indices_ra].reshape(-1)
Vmag____ = Vmag___[indices_ra].reshape(-1)
Phi_V____ = Phi_V___[indices_ra].reshape(-1)
diameter_V____ = diameter_V___[indices_ra].reshape(-1)
temps____ = temps___[indices_ra].reshape(-1)
inverse_diameter____ = inverse_diameter___[indices_ra].reshape(-1)
dist____ = dist___[indices_ra].reshape(-1)
sii_analyzed____ = sii_analyzed___[indices_ra].reshape(-1)
Dec_decimal____ = Dec_decimal___[indices_ra].reshape(-1)
RA_decimal____ = RA_decimal___[indices_ra].reshape(-1)


# Plot using the colormap based on temperatures
fig, ax1 = plt.subplots()

# Plot the first dataset with the colormap based on temperatures
sc = ax1.scatter(inverse_diameter____, Phi_V____, c=temps____, cmap=bb_cmap, marker='.', norm=norm)
sc2 = ax1.scatter(inverse_diameter____[sii_analyzed____], Phi_V____[sii_analyzed____], c=temps____[sii_analyzed____], cmap=bb_cmap, marker='*', label='SII Analyzed Stars', norm=norm)
cbar = plt.colorbar(sc, ax=ax1, label='Temperature (K)', pad=0.15)
ax1.set_yscale('log')
ax1.set_xlabel('1/θ [mas$^{-1}$]')
ax1.set_ylabel(r'Φ [photons m$^{-2}$ s$^{-1}$ Hz$^{-1}$]')
ax1.set_title('Φ vs θ')
#ax1.set_xlim(0, 25)

ax2 = ax1.twinx()
ax3 = ax1.twiny()

# Set the limits for the second y-axis based on the transformation
phi_min, phi_max = ax1.get_ylim()
thetainverse_min, thetainverse_max = ax1.get_xlim()
if thetainverse_min < 0:
    thetainverse_min = 0

ax2.set_ylim(mag_from_phi(phi_min), mag_from_phi(phi_max))
ax2.set_ylabel('magnitude')

ax3.scatter(baselines_needed____, Phi_V____, c=temps____, cmap=bb_cmap, marker='')
ax3.set_xlabel('baseline needed [m]')

ax2.grid(True)
ax3.grid(True)

st.markdown(
    """
    #### Filtered Stars
    """
)
st.write("Total number of stars taken into account with the above set filters:", len(Phi_V____))
st.markdown(
    """
    ## Plots of the filtered stars
    The stars that meet the specified criteria are shown below.
    """
)
st.pyplot(plt)




def relmag_to_absmag(rel_magnitude, distance): #distance in parsec
    return rel_magnitude + 5 - 5 * np.log10(distance)

def luminosity(absmag):
    return 10**(0.4*(4.74 - absmag))

luminosities = []
abs_Vmag = []

for i in range(len(Vmag____)):
    abs_Vmag.append(relmag_to_absmag(Vmag____[i], dist____[i]))
    
for i in range(len(Vmag____)):
    luminosities.append(luminosity(abs_Vmag[i]))


fig, ax1 = plt.subplots()
luminosities_array = np.array(luminosities)
sc = ax1.scatter(temps____, luminosities_array, c=temps____, cmap=bb_cmap, marker='.', norm=norm)
sc2 = ax1.scatter(temps____[sii_analyzed____], luminosities_array[sii_analyzed____], c=temps____[sii_analyzed____], cmap=bb_cmap, marker='*', norm=norm)


# Set labels
ax1.set_ylabel(r'Luminosity [L$_\odot$]')
ax1.set_xlabel('Temperatures [K]')
ax1.set_title('H-R Diagram')
cbar = plt.colorbar(sc, ax=ax1, label='Temperature (K)', pad=0.1)
ax1.set_yscale('log')
ax1.set_xscale('log')
plt.gca().invert_xaxis()


ax1.grid(True)

st.markdown(
    """
    ## H-R Diagram
    The corresponding H-R diagram can be seen here.
    """
)

st.pyplot(plt)

# other properties of the stars put into arrays
BayerF = np.array(data['BayerF'])
Common = np.array(data['Common'])
Parallax = np.array(data['Parallax'])
RA_decimal = np.array(data['RA_decimal'])
Dec_decimal = np.array(data['Dec_decimal'])
RA = np.array(data['RA'])
Dec = np.array(data['Dec'])
Diameter_U = np.array(data['Diameter_U'])
Diameter_B = np.array(data['Diameter_B'])
Phi_U = np.array(data['Phi_U'])
Phi_B = np.array(data['Phi_B'])
Umag = np.array(data['Umag'])
Bmag = np.array(data['Bmag'])


# filter arrays as the other properties before by baseline and magnitude
if filter_b_d == "baseline":
    BayerF_ = BayerF[indices_baseline].reshape(-1)
    Common_ = Common[indices_baseline].reshape(-1)
    Parallax_ = Parallax[indices_baseline].reshape(-1)
    RA_ = RA[indices_baseline].reshape(-1)
    Dec_ = Dec[indices_baseline].reshape(-1)
    Diameter_U_ = Diameter_U[indices_baseline].reshape(-1)
    Diameter_B_ = Diameter_B[indices_baseline].reshape(-1)
    Phi_U_ = Phi_U[indices_baseline].reshape(-1)
    Phi_B_ = Phi_B[indices_baseline].reshape(-1)
    Umag_ = Umag[indices_baseline].reshape(-1)
    Bmag_ = Bmag[indices_baseline].reshape(-1)

else:
    BayerF_ = BayerF[indices_angular_diameter].reshape(-1)
    Common_ = Common[indices_angular_diameter].reshape(-1)
    Parallax_ = Parallax[indices_angular_diameter].reshape(-1)
    RA_ = RA[indices_angular_diameter].reshape(-1)
    Dec_ = Dec[indices_angular_diameter].reshape(-1)
    Diameter_U_ = Diameter_U[indices_angular_diameter].reshape(-1)
    Diameter_B_ = Diameter_B[indices_angular_diameter].reshape(-1)
    Phi_U_ = Phi_U[indices_angular_diameter].reshape(-1)
    Phi_B_ = Phi_B[indices_angular_diameter].reshape(-1)
    Umag_ = Umag[indices_angular_diameter].reshape(-1)
    Bmag_ = Bmag[indices_angular_diameter].reshape(-1)


BayerF__ = BayerF_[indices_magnitude].reshape(-1)
Common__ = Common_[indices_magnitude].reshape(-1)
Parallax__ = Parallax_[indices_magnitude].reshape(-1)
RA__ = RA_[indices_magnitude].reshape(-1)
Dec__ = Dec_[indices_magnitude].reshape(-1)
Diameter_U__ = Diameter_U_[indices_magnitude].reshape(-1)
Diameter_B__ = Diameter_B_[indices_magnitude].reshape(-1)
Phi_U__ = Phi_U_[indices_magnitude].reshape(-1)
Phi_B__ = Phi_B_[indices_magnitude].reshape(-1)
Umag__ = Umag_[indices_magnitude].reshape(-1)
Bmag__ = Bmag_[indices_magnitude].reshape(-1)

BayerF___ = BayerF__[indices_declination].reshape(-1)
Common___ = Common__[indices_declination].reshape(-1)
Parallax___ = Parallax__[indices_declination].reshape(-1)
RA___ = RA__[indices_declination].reshape(-1)
Dec___ = Dec__[indices_declination].reshape(-1)
Diameter_U___ = Diameter_U__[indices_declination].reshape(-1)
Diameter_B___ = Diameter_B__[indices_declination].reshape(-1)
Phi_U___ = Phi_U__[indices_declination].reshape(-1)
Phi_B___ = Phi_B__[indices_declination].reshape(-1)
Umag___ = Umag__[indices_declination].reshape(-1)
Bmag___ = Bmag__[indices_declination].reshape(-1)

BayerF____ = BayerF___[indices_ra].reshape(-1)
Common____ = Common___[indices_ra].reshape(-1)
Parallax____ = Parallax___[indices_ra].reshape(-1)
RA____ = RA___[indices_ra].reshape(-1)
Dec____ = Dec___[indices_ra].reshape(-1)
Diameter_U____ = Diameter_U___[indices_ra].reshape(-1)
Diameter_B____ = Diameter_B___[indices_ra].reshape(-1)
Phi_U____ = Phi_U___[indices_ra].reshape(-1)
Phi_B____ = Phi_B___[indices_ra].reshape(-1)
Umag____ = Umag___[indices_ra].reshape(-1)
Bmag____ = Bmag___[indices_ra].reshape(-1)

# save as dataframe
filtered_data = {
    'BayerF': BayerF____,
    'Common': Common____,
    'Parallax': Parallax____,
    'Distance': dist____,
    'Umag': Umag____,  # Assuming 'Umag' was not used and therefore is not available
    'Vmag': Vmag____,
    'Bmag': Bmag____,  # Assuming 'Bmag' was not used and therefore is not available
    'Temp': temps____,
    'RA_decimal': RA_decimal____,
    'Dec_decimal': Dec_decimal____,
    'RA': RA____,
    'Dec': Dec____,
    'Diameter_U': Diameter_U____,
    'Diameter_V': diameter_V____,
    'Diameter_B': Diameter_B____,
    'Phi_U': Phi_U____,
    'Phi_V': Phi_V____,
    'Phi_B': Phi_B____,
    'SII': sii_analyzed____,
    'Baseline_Needed': baselines_needed____,
    'Inverse_Diameter': inverse_diameter____
}

df_filtered = pd.DataFrame(filtered_data)

# convert to csv
csv = df_filtered.to_csv(index=False)


st.markdown(
    """
    #### Download
    To download a CSV file of the filtered stars with its Bayer designation, common name, RA, Dec as well as other properties, click the button below:
    """
)

# Add the download button
st.download_button(
    label="Download",
    data=csv,
    file_name='filtered_stars.csv',
    mime='text/csv'
)
