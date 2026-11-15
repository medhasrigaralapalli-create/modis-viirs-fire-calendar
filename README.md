# MODIS–VIIRS Harmonized Fire Calendar

**Team:** The Cosmic Coders (Hyderabad) · NASA Space Apps Challenge 2026
**Challenge:** Harmonization of MODIS and VIIRS Hot Spots

**Live dashboard:** [paste your GitHub Pages link here]
**Demo slides:** https://docs.google.com/presentation/d/1uX2DfzrS0NdAHcLakdiZ8SIofN8ZkOht/edit?usp=sharing&ouid=111759685376896790355&rtpof=true&sd=true

## The problem
MODIS has recorded active fires from space since 2000 and is scheduled to end around 2026–27. VIIRS is its successor for fire detection, but it is not a drop-in replacement: its pixels are smaller (375 m vs 1 km), its overpass times differ, and its confidence values use a different scale (low / nominal / high vs 0–100). Without a bridge between the two, a fire trend could change only because the sensor changed.

## What this project does
1. Reads NASA FIRMS active-fire detections from MODIS and from two VIIRS satellites (NOAA-20, NOAA-21).
2. Converts everything to one schema: UTC time, confidence on a 0–1 scale, fire radiative power (FRP), pixel size.
3. Merges the two VIIRS satellites into one table.
4. Places detections on a ~2 km grid (0.02°) for each UTC day.
5. Labels each cell **MODIS only**, **VIIRS only** or **both**, with a unified confidence of `1 − (1 − c_MODIS)(1 − c_VIIRS)`. This formula is our own simple rule, not an official one.
6. Measures how often one sensor's fire cells have a fire from the other sensor within 1 cell (~2 km) the same day, with a shifted-data control for chance level, for all passes, day only and night only.
7. Builds an interactive dashboard (map, day selector, sensor toggles, agreement table, benchmark and fire-intensity charts).

## Results (central India, 1–7 Nov 2026; agreement statistics use full days, 2–6 Nov)
| Measure | Value |
|---|---|
| Detections | VIIRS 4,158 (NOAA-20 1,969 + NOAA-21 2,189) vs MODIS 314 |
| MODIS cells with a VIIRS fire within ~2 km, same day | 42.2% overall, 46.9% day only |
| Control (VIIRS shifted 0.5°, chance level) | 4.6% overall, 4.1% day only |
| VIIRS cells with a MODIS fire within ~2 km, same day | 8.5% overall, 9.0% day only |
| Reference level: NOAA-21 cells seen by NOAA-20 / NOAA-20 cells seen by NOAA-21 | 36.8% / 30.1% (control 4.1%) |
| VIIRS cells seen by MODIS, by fire intensity (FRP, MW) | 3.8% (<2), 7.1% (2–4), 8.2% (4–8), 14.9% (8–16), 19.5% (16–32), 22.7% (32–64, 22 cells) |

Reading the results: VIIRS finds about 13 times more detections. The sensors agree far more than chance, but only partly: even two satellites carrying the same VIIRS instrument agree on only about a third of fire cells, probably because they pass at different times and many fires are short-lived (we did not test this explanation). The calendar keeps cells seen by only one sensor, so neither sensor's gaps are mistaken for fire-free areas.

## How to run it
1. Install Python 3.9 or newer and the libraries: `pip install pandas numpy`
2. Put the FIRMS files (`.csv`, or the downloaded `.zip` files) in the `data/` folder. MODIS vs VIIRS is detected automatically, so file names do not matter.
3. Run: `python run.py`
4. This writes `output/fire_calendar.csv` and `output/index.html` and opens the dashboard. The data is embedded in `index.html`, so that single file can be shared.

Options: `python run.py --cell 0.05` (coarser grid), `--keep-edge-days`, `--no-open`, `--demo` (synthetic data to test that everything runs; results are not real).

The map uses Leaflet and Esri basemap tiles, which need an internet connection. Use the layer button on the map to switch basemaps or choose "No basemap".

## Data
- NASA FIRMS (Fire Information for Resource Management System), near-real-time active-fire data: MODIS (Terra and Aqua, Collection 6.1), VIIRS NOAA-20 and NOAA-21 (375 m). https://firms.modaps.eosdis.nasa.gov/
- Files in `data/` are the unmodified downloads for the area about 76–84°E, 17–24°N, 1–7 Nov 2026.

## Limitations
- One week and one region, so results are indicative.
- "Confirmed by the other sensor" measures sensor agreement. It is not validated against reference fire records.
- Matching is by UTC date, not exact overpass time.
- MODIS has only 52 night detections in this download, so night results are indicative.
- The first and last days (1 and 7 Nov) are partial for one sensor and are left out of the statistics.
- The 2 km grid, the 1-cell neighbourhood and the unified-confidence formula are design choices.
- Near-real-time data can be slightly revised in the later standard-processing version.

## Next steps
Longer periods and other regions, validation against independent fire records, and a model that predicts whether a detection is confirmed by the other sensor.

## Credits
NASA FIRMS data, Leaflet, Esri basemap tiles, Python (pandas, NumPy).
