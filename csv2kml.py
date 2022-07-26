#!/usr/bin/env python3

import os
import pandas as pd
import exifread

# Copyright (C) 2022  Yingbo Li
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from kml_funcs import *
from coord_funcs import *
import config

def img_taken_time(img_path):
    with open(img_path,"rb") as infile:
        tags = exifread.process_file(infile, stop_tag="EXIF DateTimeOriginal")

    try:
        times = [string.split(":") for string in str(tags["EXIF DateTimeOriginal"]).split(" ")]
        times = times[0] + times[1]
        year,month,day,hour,minute,second = times
        time_str = "/".join([day,month,year]) + " " + ":".join([hour,minute,second])
    except:
        time_str = ""
    return time_str

df = pd.read_csv(config.CSV_FILE).fillna("")

all_placemarks = ""
for _,row in df.iterrows():
    name = row["id"]
    grid_ref = row["grid ref"]
    E,N = grid_ref2osgb(grid_ref)
    lat,lon,_ = osgbEN2wgsLatLon(E,N)
    img_name = row["img"]
    desc = row["desc"].strip()

    if img_name != "":
        img_path = os.path.join(config.IMG_DIR,img_name)
        time_taken = img_taken_time(img_path)
        contents = '<img src="%s" width="250px"></br></br>Time Taken: %s</br>%s' % (img_path,time_taken,desc)
    else:
        contents = str(desc)


    all_placemarks += makePlacemark(name,
                                    "%.15f,%.15f" % (lon,lat),
                                    contents)

kml = initKml(kmlFolder(config.KML_FOLDER_NAME,all_placemarks))

with open(config.KML_FILE_NAME,"w") as outfile:
    outfile.write(kml)
