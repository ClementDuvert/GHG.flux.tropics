# Dataset Documentation
This folder contains raw and processed datasets used for the analysis of GHG fluxes from (sub)tropical lakes and rivers, including supporting metadata.

---

## References

| File name | Description |
|-----------|-------------|
| `reference_list.csv` | Bibliographic references linked to original data sources for measurements |

---

## Datasets overview

| File name | Description |
|-----------|-------------|
| `dataset_lakes_raw.csv` | Raw extraction of lake data |
| `dataset_rivers_raw.csv` | Raw extraction of river data |
| `dataset_lakes_reduced.csv` | Processed and reduced version of lake dataset with means per system |
| `dataset_rivers_reduced.csv` | Processed and reduced version of river dataset with means per site |

---

## River metadata

| File name | Description |
|-----------|-------------|
| `river_properties_BasinATLAS_2025.csv` | Catchment properties from BasinATLAS |
| `river_properties_GRADES_2025.csv` | GRADES-derived river network properties |
| `river_properties_H90_2025.csv` | H90 stream order and channel properties |
| `river_properties_npp_peat_2025.csv` | NPP and peatland extent |
| `river_surface_areas.csv` | Surface areas of rivers per climate zone and stream order |

---

## Lake metadata

| File name | Description |
|-----------|-------------|
| `lake_HydroLAKES.2025.txt` | Lake attribute data from HydroLAKES (area, depth, etc.) |
| `lake_manual_2025.csv` | Manually curated lake metadata |
| `lake_surface_areas_natural.csv` | Surface areas of natural lakes and ponds per climate zone and size class |
| `lake_surface_areas_manmade.csv` | Surface areas of reservoirs and man-made ponds per climate zone and size class |
| `lake_wind_2025.csv` | Mean annual wind speed |

---

## Column descriptions and units

### raw datasets (`dataset_lakes_raw.csv` and `dataset_rivers_raw.csv`)

| Column name             | Units / format      | Description |
|-------------------------|---------------------|-------------|
| `Ref`                   | —                   | Reference number as per `reference_list.csv` |
| `Lake_idx`              | —                   | Internal lake identifier |
| `LakeType`              | —                   | Type of waterbody: "lake" or "reservoir" |
| `Site name`             | —                   | Name of the sampling site |
| `day`                   | DD                  | Day of measurement |
| `month`                 | MM                  | Month of measurement |
| `year`                  | YYYY                | Year of measurement |
| `time`                  | HH:MM:SS            | Time of measurement |
| `measurement_frequency` | —                   | Frequency of measurement as per main text in study (from "interannual" to "HF"; see Methods) |
| `reported_frequency`    | —                   | Frequency of available data (from "interannual" to "HF"; see Methods) |
| `run_on_river`          | YES/NO              | Whether the site is part of a run-on-river survey |
| `Latitude`              | decimal degrees     | Latitude of sampling location |
| `Longitude`             | decimal degrees     | Longitude of sampling location |
| `Elevation`             | masl                | Elevation above mean sea level |
| `Climate`               | categorical (1–5)   | 1: humid tropics; 2: wet-dry tropics; 3: (semi)arid (sub)tropics; 4: humid subtropics; 5: highland (sub)tropics |
| `Slope`                 | m/m                 | Local slope of the channel |
| `Watershed Area`        | km2                 | Total watershed (catchment) area |
| `River width`           | m                   | Width of the river channel |
| `Surface area_m2`       | m2                  | Surface area of lake or reservoir |
| `Volume`                | m3                  | Volume of lake or reservoir |
| `Mean depth`            | m                   | Mean depth of lake or reservoir |
| `Max depth`             | m                   | Maximum depth of lake or reservoir |
| `Secchi depth`          | m                   | Secchi disk transparency depth |
| `Sample depth`          | m                   | Depth at which sample was collected |
| `CO2`                   | —                   | CO2 concentration |
| `CO2_units`             | —                   | Units used for CO2 |
| `CO2_method_cat`        | —                   | Method used to measure CO2: direct (sensor, headspace) or indirect (alkalinity + pH) |
| `CH4`                   | —                   | CH4 concentration |
| `CH4_units`             | —                   | Units used for CH4 |
| `N2O`                   | —                   | N2O concentration |
| `N2O_units`             | —                   | Units used for N2O |
| `F_method_cat`          | —                   | Method used to estimate GHG fluxes: measured or modelled |
| `F_CO2`                 | —                   | Diffusive CO2 flux |
| `F_CO2_units`           | —                   | Units used for CO2 flux |
| `FCH4_diff_avg_val`     | —                   | Diffusive CH4 flux |
| `F_CH4_units`           | —                   | Units used for CH4 flux |
| `FCH4_ebb_avg_val`      | —                   | Ebullitive CH4 flux |
| `F_CH4_ebb_units`       | —                   | Units used for ebullitive CH4 flux |
| `FN2O_avg_val`          | —                   | Diffusive N2O flux |
| `FN2O_UNITS`            | —                   | Units used for N2O flux |
| `Strahler order`        | —                   | Stream order based on Strahler classification |
| `Discharge_mean`        | —                   | River discharge |
| `Discharge_mean_units`  | —                   | Units used for discharge |
| `Velocity`              | —                   | Flow velocity |
| `Velocity_units`        | —                   | Units used for velocity |
| `Temperature`           | °C                  | Water temperature |
| `DO`                    | —                   | Dissolved oxygen |
| `DO_units`              | —                   | Units used for DO |
| `pH`                    | —                   | pH |
| `Conductivity`          | µS/cm               | Electrical conductivity |
| `DOC_value`             | —                   | Dissolved organic carbon |
| `DOC_units`             | —                   | Units used for DOC |
| `DIC_value`             | —                   | Dissolved inorganic carbon |
| `DIC_units`             | —                   | Units used for DIC |
| `k`                     | —                   | Gas transfer velocity |
| `k600`                  | —                   | Standardised gas transfer velocity at 600 Schmidt |
| `k_units`               | —                   | Units used for `k` and `k600` |
| `Comments`              | —                   | Any notes on data extraction, averaging, or corrections |


### reduced datasets (`dataset_lakes_reduced.csv` and `dataset_rivers_reduced.csv`)

| Column name                     | Units / format        | Description |
|--------------------------------|------------------------|-------------|
| `Refs`                           | —                      | Aggregated reference numbers as per `reference_list.csv` |
| `Site`                           | —                      | Aggregated site names |
| `Latitude`                       | decimal degrees        | Latitude of site |
| `Longitude`                      | decimal degrees        | Longitude of site |
| `Elevation_m`                    | m                      | Elevation based on H90 |
| `CatchmentArea_km2`              | km2                    | Catchment area based on H90 |
| `StreamOrder`                    | —                      | Consolidated stream order (combined from reported and H90-derived values) |
| `LakeIndex`                      | —                      | Lake index |
| `LakeType`                       | —                      | Type of waterbody: "lake" or "reservoir" |
| `Slope`                          | m/m                    | Consolidated slope (combined from reported and H90-derived values) |
| `Climate`                        | categorical (1–5)      | 1: humid tropics; 2: wet-dry tropics; 3: (semi)arid (sub)tropics; 4: humid subtropics; 5: highland (sub)tropics |
| `k600_md`                        | m/d                    | Mean consolidated gas transfer velocity (k600), measured or modelled |
| `MeasFreq`                       | —                      | Highest frequency of measurement across studies at this site (from "interannual" to "HF"; see Methods) |
| `ReportedFreq`                   | —                      | Highest frequency of available data at this site (from "interannual" to "HF"; see Methods) |
| `DirectIndirect_CO2`             | —                      | CO2 data source: "direct" (sensor, headspace) or "indirect" (alkalinity+pH) |
| `FluxMethod`                     | —                      | Method used to estimate GHG fluxes: "measured" or "modelled" |
| `Discharge_m3s`                  | m3/s                   | Mean annual runoff from GRADES |
| `Velocity_ms`                    | m/s                    | Mean flow velocity |
| `CO2_uM`                         | uM                     | Mean CO2 concentration |
| `n_CO2`                          | —                      | Number of individual measurements used for site-level averaging |
| `F_CO2_mmol_m2_d`                | mmol m-2 d-1           | Mean diffusive CO2 flux |
| `n_FCO2`                         | —                      | Number of individual measurements used for site-level averaging |
| `F_CO2_measmod_mmol_m2_d`        | mmol m-2 d-1           | Mean CO2 flux estimate (combined measured and modelled) |
| `CH4_uM`                         | uM                     | Mean CH4 concentration |
| `n_CH4`                          | —                      | Number of individual measurements used for site-level averaging |
| `F_CH4_diff_mmol_m2_d`           | mmol m-2 d-1           | Mean diffusive CH4 flux |
| `n_FCH4diff`                     | —                      | Number of individual measurements used for site-level averaging |
| `F_CH4_measmod_mmol_m2_d`        | mmol m-2 d-1           | Mean diffusive CH4 flux estimate (combined measured and modelled) |
| `F_CH4_eb_mmol_m2_d`             | mmol m-2 d-1           | Mean ebullitive CH4 flux |
| `n_FCH4eb`                       | —                      | Number of individual measurements used for site-level averaging |
| `N2O_uM`                         | uM                     | Mean N2O concentration |
| `n_N2O`                          | —                      | Number of individual measurements used for site-level averaging |
| `F_N2O_mmol_m2_d`                | mmol m-2 d-1           | Mean diffusive N2O flux |
| `n_FN2O`                         | —                      | Number of individual measurements used for site-level averaging |
| `F_N2O_measmod_mmol_m2_d`        | mmol m-2 d-1           | Mean N2O flux estimate (combined measured and modelled) |
| `Temp_C`                         | °C                     | Mean water temperature |
| `DO_mgL`                         | mg/L                   | Mean dissolved oxygen |
| `pH`                             | —                      | Mean pH |
| `Cond_uScm`                      | uS/cm                  | Mean electrical conductivity |
| `DOC_mgL`                        | mg/L                   | Mean dissolved organic carbon |
| `DIC_mgL`                        | mg/L                   | Mean dissolved inorganic carbon |

---

## License and usage

These datasets are released under the **[CC-BY 4.0 License](https://creativecommons.org/licenses/by/4.0/)**. You are free to use, share, and adapt the data with proper attribution.

If using this dataset in a publication, please cite:
Duvert et al. (in revision). Hydroclimate and landscape diversity drive highly variable greenhouse gas emissions from (sub)tropical inland waters. Nature Water. 

---

## Contact
For any inquiry, please contact:
**Clément Duvert**  clem.duvert@cdu.edu.au 
Charles Darwin University  
ORCID: [0000-0002-9873-6846](https://orcid.org/0000-0002-9873-6846)

