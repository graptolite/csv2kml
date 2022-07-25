# Usage
## Data Input
Store data in `data.csv` in table format:
| id | grid ref | img | desc |

- *id*: \[required\] Display name of the site.
- *grid ref*: \[required\] Grid reference (location) of the site. Note that the grid reference must contain two letters at the start and no more than 5 numbers each numerical easting or northing component (e.g. TQ 10000 10000 is fine but TQ 100001 100002 and TQ 10000.1 10000.2 aren't).
- *img*: \[optional\] Image to accompany the site.
- *desc*: \[optional\] Description to accompany the site.

## Configuation
Edit config.py as necessary. All of these are required but have presets values already.
 
- *CSV_FILE*: name of csv file with the necessary column headers
- *KML_FOLDER_NAME*: name of the KML folder that will hold all the points of interest.
- *KML_FILE_NAME*: name of the file the kml is written to.
- *IMG_DIR*: folder that holds all the images - absolute paths are preferred for this.

## Dependencies
- os
- string
- pandas
- exifread
- numpy