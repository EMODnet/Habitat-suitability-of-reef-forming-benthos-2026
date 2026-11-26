# Habitat suitability of reef forming benthos - version 2026

## Introduction

This product is an update of the earlier product "Habitat suitability of reef-forming species in the North Sea" (https://github.com/EMODnet/EMODnet-Biology-Reef-Forming-Species-NorthSea). For Lanice conchilega and Sabellaria spinulosa, the prediction was extended to cover a wide geographic scope, encompassing the North Sea, English Channel, Irish and Celtic Seas.

## Directory structure

```
{{directory_name}}/
├── analysis
├── data/
│   ├── derived_data/
│   └── raw_data/
├── docs/
├── product/
└── scripts/
```

* **analysis** - Markdown or Jupyter notebooks
* **data** - Raw and derived data
* **docs** - Rendered reports
* **product** - Output product files
* **scripts** - Reusable code

## Data series

Habitat suitability models combine species distribution data (Presence/Absence) and environmental data (rasters). A complete description of the required data is given in the Readme file of "Habitat suitability of reef-forming species in the North Sea" (https://github.com/EMODnet/EMODnet-Biology-Reef-Forming-Species-NorthSea).

Bathymetry (2024) was extracted from EMODnet and used for the calculation of BPI values (Calculate BPI.R). Temperature, salinity and mixed layer depth (as a measure for stratification) were downloaded from Copernicus. Bottom shear stress and chlorophyll a concentrations were extracted from the Deltares 3D DCSM-FM model. Predicted gravel, sand and mud fractions from Mitchell et al. (2019) were used (https://doi.org/10.14466/CefasDataHub.63). All environmental data were saved as GeoTIFF rasters.

Species distribution data were extracted from EMODnet (EMODnet_presabs.R) and saved as a point shapefile (Presence/Absence). OBIS records were added to account for missing data in the German EEZ.

```
{{data_wfs_request}}
```

## Data product

The updated product contains habitat suitability maps (generated with Random Forest) for Ostrea edulis and Modiolus modiolus in the central and southern North Sea, and for Sabellaria spinulosa and Lanice conchilega in the entire North Sea, English Channel and Irish and Celtic Seas.

## More information:

### References

* P.J. Mitchell, J. Aldridge, and M. Diesing (2019). Legacy data: How decades of seabed sampling can produce robust predictions and Versatile Products. Geosciences 9(4): 182. https://doi.org/10.3390/geosciences9040182
* P.M.J. Herman, and F.F. van Rees (2022). Mapping reef forming North Sea species. Deltares report, Deltares, Delft, the Netherlands. https://pub.kennisbank.deltares.nl/Details/fullCatalogue/1000020826.

### Code and methodology

Code can be found in the "scripts" folder. Further explanation on the methodology is given in Herman and van Rees (2022) and the associated data product "Habitat suitability of reef-forming species in the North Sea" (https://github.com/EMODnet/EMODnet-Biology-Reef-Forming-Species-NorthSea).

### Citation and download link

This product should be cited as:

S. Mestdagh, and W. Stolte (2026). Habitat suitability for reef-forming macrobenthic species in the North Sea. Integrated EMODnet Biology data product. Deltares, Delft, the Netherlands.

Available to download in:

{{link_download}}

### Authors

Sebastiaan Mestdagh & Willem Stolte

contact: sebastiaan.mestdagh@deltares.nl
