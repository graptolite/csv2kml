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

def make_xml(tag,contents=""):
    ''' if options are needed for the tag, separate them with a space from the tag name e.g. tag='div style="width:100%;"'
    '''
    def tag_str(tag,contents):
        tag_only = tag.split(" ")[0]
        return "<%s>\n%s\n</%s>" % (tag,contents,tag_only)

    if isinstance(tag,str):
        out = tag_str(tag,contents)
    elif isinstance(tag,list):
        if len(tag) > 1:
            out = make_xml(tag[:-1],tag_str(tag[-1],contents))
        else:
            out = tag_str(tag[0],contents)
    else:
        print("Warning: Tag doesn't seem to be a normal ML tag (should be a string or list of strings)")
    return out

def make_cdata(contents=""):
    ## declare contents as character data e.g. not parsed as the "X"ML
    return "<![CDATA[\n%s\n]]>" % contents

def init_kml(contents=""):
    ## shove some contents into configured kml and document tags
    return make_xml(['kml xmlns="http://www.opengis.net/kml/2.2"',"Document"],contents)

def kml_folder(name,contents=""):
    return make_xml("Folder",make_xml("name",name)+contents)

def make_placemark(name,coordinates,description="",icon_scale=1.0,icon_marker="http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png"):
    ''' coordinates in format:
        lon,lat e.g. -2.334,49.294
    '''
    return make_xml("Placemark",
                   make_xml("name",name)+
                   make_xml("description",make_cdata(description))+
                   make_xml(["Style","IconStyle"],
                           make_xml("scale",icon_scale)+
                           make_xml(["Icon","href"],icon_marker))+
                   make_xml(["Point","coordinates"],coordinates))
