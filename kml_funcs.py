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

def makeXml(tag,contents=""):
    ''' if options are needed for the tag, make the tag into a dictionary:
        {"tag_name":"tag-options"}
        This can be applied to mix-and-match lists as well:
        [{"tag_name":"tag-options"},"tag_name",{"tag_name":"tag-options"}]
    '''
    def tag_str(tag,contents):
        if type(tag) == str:
            tag_only = tag
            option = ""
        elif type(tag) == dict:
            tag_only = list(tag.keys())[0]
            option = " %s" % list(tag.values())[0]

        return "<%s%s>\n%s\n</%s>\n" % (tag_only,option,contents,tag_only)

    if type(tag) == str:
        out = tag_str(tag,contents)
    elif type(tag) == list:
        if len(tag) > 1:
            out = makeXml(tag[:-1],tag_str(tag[-1],contents))
        else:
            out = tag_str(tag[0],contents)
    else:
        print("Warning: Tag doesn't seem to be a normal ML tag (should be a string or list of strings)")
    return out

def makeCdata(contents=""):
    ## declare contents as character data e.g. not parsed as the "X"ML
    return "<![CDATA[\n%s\n]]>" % contents

def initKml(contents=""):
    ## shove some contents into configured kml and document tags
    return makeXml([{"kml":'xmlns="http://www.opengis.net/kml/2.2"'},"Document"],contents)

def kmlFolder(name,contents=""):
    return makeXml("Folder",makeXml("name",name)+contents)

def makePlacemark(name,coordinates,description="",icon_scale=1.0,icon_marker="http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png"):
    ''' coordinates in format:
        lon,lat e.g. -2.334,49.294
    '''
    return makeXml("Placemark",
                   makeXml("name",name)+
                   makeXml("description",makeCdata(description))+
                   makeXml(["Style","IconStyle"],
                           makeXml("scale",icon_scale)+
                           makeXml(["Icon","href"],icon_marker))+
                   makeXml(["Point","coordinates"],coordinates))
