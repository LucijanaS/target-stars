# target-stars

A [Streamlit](https://streamlit.io) app for exploring which stars are bright enough and large
enough (angular diameter) to be observed with intensity interferometry. Filter by baseline
length or angular diameter, magnitude, and sky position, and download the resulting star list.

## Data

`combined_stars.csv` merges two catalogues:

- The **Yale Bright Star Catalogue**, with angular diameters reconstructed from V-band
  magnitude and effective temperature via a blackbody flux relation.
- **Gaia DR3** (stars brighter than G=9 with a resolved GSP-Phot radius), whose angular
  diameters come directly from Gaia's own radius + distance fit rather than a flux
  reconstruction. This is ~98% of the catalog.

Duplicate stars (matched by sky position) are kept once, preferring the Bright Star Catalogue
entry for its real Johnson photometry. `observation_log.csv` records which stars have already
been observed via intensity interferometry (currently the 32 stars from Hanbury Brown, Davis &
Allen 1974) and is already merged into `combined_stars.csv`'s `sii_observed` column -- it's kept
here for reference, not read directly by the app.

Both files are built in the companion [brightstar](https://github.com/LucijanaS/brightstar)
repository; see its
[DATA_SOURCES.md](https://github.com/LucijanaS/brightstar/blob/main/DATA_SOURCES.md) for the
exact ADQL/SIMBAD queries and `build_combined_catalog.py` for the build script. Re-copy
`combined_stars.csv` (and `observation_log.csv`) here whenever that pipeline is re-run.

`blackbody_colors` is a star-color-by-temperature lookup table (source:
http://www.vendian.org/mncharity/dir3/blackbody/), used to color plot points by temperature.

## Running locally

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```
