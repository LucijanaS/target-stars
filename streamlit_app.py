import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.constants import c, h, pi
from scipy.special import j1, jn_zeros

st.markdown(
    """
    # Stars available for Intensity Interferometry

    On this webpage one can explore and plot the HR-diagram of available stars that could be observed using intensity interferometry.
    """
    )

plt.style.use('dark_background')


@st.cache_data
def load_catalog():
    return pd.read_csv('combined_stars.csv', low_memory=False)


df = load_catalog()

# Only the rows a star needs to appear in the Phi-vs-theta plot and the baseline/diameter
# filters below -- distance is checked separately, only where it's actually needed (the H-R
# diagram), so a missing/unmeasured distance doesn't drop an otherwise-usable star.
core_valid = df['theta_mas'].notna() & df['temp_K'].notna() & df['mag'].notna() & (df['theta_mas'] > 0)
df = df[core_valid].reset_index(drop=True)


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

# The blackbody_colors table only tabulates up to 29800 K (hotter stars are all much the same
# saturated blue-white), so the color norm stops there too rather than stretching to this
# catalog's hottest star (~45000 K) -- anything past 29800 K just clips to that end color.
norm = plt.Normalize(vmin=0, vmax=29800)


def mas_to_rad(theta_mas):
    """Convert an angle in milliarcseconds to radians."""
    return theta_mas / 1000 * pi / (3600 * 180)


def baseline_needed(theta_mas, wavelength_m):
    """Minimum interferometric baseline (m) to resolve the first visibility null for a uniform
    disk of angular diameter theta_mas (mas) at the given wavelength (m). Both arguments may be
    arrays -- used here with each star's own catalogue wavelength (V-band for BSC, G-band for
    Gaia) rather than a single hardcoded wavelength, since the combined catalog mixes both."""
    theta_rad = mas_to_rad(theta_mas)
    j1_root = jn_zeros(1, 1)[0]  # jn_zeros returns an array -- [0] to get the scalar root
    return j1_root / (pi * theta_rad / wavelength_m)


def mag_from_phi(phi, wavelength_m):
    """Inverse of Phi(): recover an apparent magnitude from a spectral photon flux density."""
    nu = c / wavelength_m
    return -2.5 * (22.44 + np.log10(2 * nu * h * phi))


def relmag_to_absmag(rel_magnitude, distance_pc):
    """Converts an apparent magnitude to absolute magnitude given a distance in parsec."""
    return rel_magnitude + 5 - 5 * np.log10(distance_pc)


def luminosity_from_absmag(absmag):
    """Bolometric-ish luminosity (L_sun) from an absolute magnitude, calibrated against the Sun (M_V=4.74)."""
    return 10 ** (0.4 * (4.74 - absmag))


wavelength_m = df['wavelength_nm'].to_numpy() * 1e-9
df['baseline_needed_m'] = baseline_needed(df['theta_mas'].to_numpy(), wavelength_m)
df['inverse_theta_mas'] = 1 / df['theta_mas']


def plot_phi_vs_theta(data):
    # Wider/taller than matplotlib's default (6.4x4.8in) -- this plot carries four axis label
    # sets (bottom/left/right/top) plus a colorbar, all sharing one figure, so the default size
    # leaves too little room for the actual data area and every label ends up looking oversized.
    fig, ax1 = plt.subplots(figsize=(9, 6.5))
    ax1.scatter(data['inverse_theta_mas'], data['phi'], c=data['temp_K'], cmap=bb_cmap, marker='.', norm=norm)
    marked = data['sii_observed'].to_numpy(dtype=bool)
    ax1.scatter(data['inverse_theta_mas'][marked], data['phi'][marked], c=data['temp_K'][marked],
                cmap=bb_cmap, marker='*', label='SII Observed Stars', norm=norm)
    sc = ax1.scatter([], [], c=[], cmap=bb_cmap, norm=norm)  # dummy mappable for the colorbar
    plt.colorbar(sc, ax=ax1, label='Temperature (K)', pad=0.15)
    ax1.set_yscale('log')
    ax1.set_xlabel('1/θ [mas$^{-1}$]')
    ax1.set_ylabel(r'Φ [photons m$^{-2}$ s$^{-1}$ Hz$^{-1}$]')
    ax1.set_title('Φ vs θ')

    ax2 = ax1.twinx()
    ax3 = ax1.twiny()

    phi_min, phi_max = ax1.get_ylim()
    # The G-band (Gaia, 622 nm) accounts for ~98% of this catalog, so the secondary magnitude
    # axis uses that wavelength -- it's only approximate for the V-band (BSC) points mixed in.
    ax2.set_ylim(mag_from_phi(phi_min, 622e-9), mag_from_phi(phi_max, 622e-9))
    ax2.set_ylabel('magnitude (approx., G-band)', fontsize=9)
    ax2.tick_params(labelsize=9)

    ax3.scatter(data['baseline_needed_m'], data['phi'], c=data['temp_K'], cmap=bb_cmap, marker='')
    ax3.set_xlabel('baseline needed [m]', fontsize=9)
    ax3.tick_params(labelsize=9)

    ax2.grid(True)
    ax3.grid(True)
    fig.tight_layout()
    return fig


st.markdown(
    """
    The stars used are combined from two catalogues -- see [DATA_SOURCES.md](https://github.com/LucijanaS/brightstar/blob/main/DATA_SOURCES.md)
    in the companion [brightstar](https://github.com/LucijanaS/brightstar) repository for the exact queries and how they were merged:

    - The **Yale Bright Star Catalogue** (9110 of the brightest stars, http://tdc-www.harvard.edu/catalogs/bsc5.html,
      converted to JSON by https://github.com/brettonw/YaleBrightStarCatalog). Angular diameters for these stars are
      reconstructed from V-band magnitude and effective temperature via the blackbody flux relation.
    - **Gaia DR3** (ESA), queried for stars brighter than G=9 with a resolved GSP-Phot radius. Angular diameters for
      these stars come directly from Gaia's own radius and distance fit, not a flux reconstruction -- this is the
      large majority of the catalog and is generally more reliable than the BSC's flux-based diameters.

    Duplicate stars (matched by sky position) are kept once, preferring the Bright Star Catalogue entry since it
    carries real Johnson photometry.
    """
    )

st.write(
    f"Of the stars in the combined catalog, {len(df)} had enough information (angular diameter, "
    f"temperature, and magnitude) to plot below."
)
st.write("The ones marked with a '★' have already been observed via intensity interferometry "
         "-- currently the 32 stars from Hanbury Brown, Davis & Allen (1974).")

st.pyplot(plot_phi_vs_theta(df))
plt.close('all')


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

if filter_b_d == "baseline":
    st.markdown("#### Baseline Selection")
    baseline_available = st.slider("Baseline available in meters", 0, 3500, (0, 3500))
    baseline_min, baseline_max = baseline_available
    mask = (df['baseline_needed_m'] >= baseline_min) & (df['baseline_needed_m'] <= baseline_max)
else:
    st.markdown("#### Angular Diameter Selection")
    desired_angular_diameter = st.slider("Angular diameter of star desired in milliarcseconds", 0.00, 50.00, (0.00, 50.00))
    diameter_min, diameter_max = desired_angular_diameter
    mask = (df['theta_mas'] >= diameter_min) & (df['theta_mas'] <= diameter_max)

st.markdown(
    """
    #### Minimum Magnitude

    Set the minimum magnitude of the star you want to observe. This helps filter out stars that are too faint.
    Note: this is V-band magnitude for Bright Star Catalogue stars and G-band magnitude for Gaia stars -- the two
    aren't perfectly comparable, but both track observability well enough for filtering.
    """
)
magnitude_min = st.number_input("Minimum magnitude of star", -2, 8, value=8)
mask &= df['mag'] < magnitude_min

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
    declination_range = (-90, 90)

declination_min, declination_max = declination_range
mask &= (df['dec_deg'] >= declination_min) & (df['dec_deg'] <= declination_max)

st.markdown(
    """
    #### Right Ascension Range

    Toggle the option below if you want to specify a range for right ascension.
    """
)
on_ra = st.toggle("Specify Right Ascension Range")

if on_ra:
    ra_start = st.number_input("Minimum Right Ascension in hours", 0, 24, value=0)
    ra_end = st.number_input("Maximum Right Ascension in hours", 0, 24, value=24)
    st.write(f"RA range: {ra_start}h to {ra_end}h")
else:
    ra_start, ra_end = 0, 24

ra_deg_min, ra_deg_max = ra_start * 15.0, ra_end * 15.0
mask &= (df['ra_deg'] >= ra_deg_min) & (df['ra_deg'] <= ra_deg_max)

filtered = df[mask].reset_index(drop=True)

st.markdown(
    """
    #### Filtered Stars
    """
)
st.write("Total number of stars taken into account with the above set filters:", len(filtered))
st.markdown(
    """
    ## Plots of the filtered stars
    The stars that meet the specified criteria are shown below.
    """
)
st.pyplot(plot_phi_vs_theta(filtered))
plt.close('all')


st.markdown(
    """
    ## H-R Diagram
    The corresponding H-R diagram can be seen here. Only stars with a known distance are included, since
    luminosity requires converting apparent to absolute magnitude.
    """
)

hr_valid = filtered['distance_pc'].notna() & (filtered['distance_pc'] > 0)
hr = filtered[hr_valid].copy()
hr['abs_mag'] = relmag_to_absmag(hr['mag'], hr['distance_pc'])
hr['luminosity_Lsun'] = luminosity_from_absmag(hr['abs_mag'])

fig, ax1 = plt.subplots()
ax1.scatter(hr['temp_K'], hr['luminosity_Lsun'], c=hr['temp_K'], cmap=bb_cmap, marker='.', norm=norm)
hr_marked = hr['sii_observed'].to_numpy(dtype=bool)
ax1.scatter(hr['temp_K'][hr_marked], hr['luminosity_Lsun'][hr_marked], c=hr['temp_K'][hr_marked],
            cmap=bb_cmap, marker='*', norm=norm)
sc = ax1.scatter([], [], c=[], cmap=bb_cmap, norm=norm)
plt.colorbar(sc, ax=ax1, label='Temperature (K)', pad=0.1)

ax1.set_ylabel(r'Luminosity [L$_\odot$]')
ax1.set_xlabel('Temperatures [K]')
ax1.set_title('H-R Diagram')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.invert_xaxis()
ax1.grid(True)

st.pyplot(fig)
plt.close('all')

st.markdown(
    """
    #### Download
    To download a CSV file of the filtered stars with its name, RA, Dec as well as other properties, click the button below:
    """
)

st.download_button(
    label="Download",
    data=filtered.to_csv(index=False),
    file_name='filtered_stars.csv',
    mime='text/csv'
)
