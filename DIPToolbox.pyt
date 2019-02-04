# vim: filetype=python tabstop=4 shiftwidth=4 expandtab smarttab autoindent
############################################################################
# Name:    Digital Image Processing Toolbox for ArcGIS Pro
# Purpose: Provides various digital image processing tools.
# Author:  Huidae Cho
# Since:   February 3, 2019
#
# Copyright (C) 2019, Huidae Cho <https://idea.isnew.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#############################################################################

# path to custom modules
# copy agpy.py and dippy.py to this folder
modules_path = 'T:\\'
# path to output files
output_path = 'T:\\'

import arcpy
import os
import time

import sys
sys.path.append(modules_path)

# https://github.com/HuidaeCho/agpy
import agpy

# https://github.com/HuidaeCho/dippy
import dippy

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Digital Image Processing Toolbox"
        self.alias = "DIP"

        # List of tool classes associated with this toolbox
        self.tools = [
            GrayscaleTransform,
            NegativeTransform,
            LinearTransform,
            LogTransform,
            InverseLogTransform,
            PowerTransform,
            GrayLevelSlice,
            BitPlaneSlice,
            HistogramEqualize,
        ]

class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return

def get_map(suffix):
    tif_path = output_path + 'dip-map-' + suffix + '.tif'
    if os.path.exists(tif_path):
        os.remove(tif_path)
    agpy.export_map_to_tiff(tif_path, 1024, 1024)
    img = arcpy.Raster(tif_path)

    tif_path = output_path + 'dip-old-' + suffix + '.tif'
    if os.path.exists(tif_path):
        os.remove(tif_path)
    img_a = agpy.convert_raster_to_grayscale_array(img)
    img = agpy.convert_array_to_raster(img_a, img)
    arcpy.CopyRaster_management(img, tif_path)
    arcpy.mp.ArcGISProject('CURRENT').activeMap.addDataFromPath(tif_path)
    return img, img_a

def save_map_array(suffix, img_a, ref_img):
    tif_path = output_path + 'dip-new-' + suffix +'.tif'
    if os.path.exists(tif_path):
        os.remove(tif_path)
    img = agpy.convert_array_to_raster(img_a, ref_img)
    arcpy.CopyRaster_management(img, tif_path)
    arcpy.mp.ArcGISProject('CURRENT').activeMap.addDataFromPath(tif_path)

class GrayscaleTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Grayscale transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        new_L = arcpy.Parameter(
            name='new_L',
            displayName='New L',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        params = [new_L]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        new_L = parameters[0].value

        suffix = 'grayscale-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.grayscale_transform(img_a, new_L)

        save_map_array(suffix, new_img_a, img)
        return

class NegativeTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Negative transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        suffix = 'negative-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.negative_transform(img_a)

        save_map_array(suffix, new_img_a, img)
        return

class LinearTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Linear transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        cp1r = arcpy.Parameter(
            name='cp1r',
            displayName='Control point 1 r',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        cp1s = arcpy.Parameter(
            name='cp1s',
            displayName='Control point 1 s',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        cp2r = arcpy.Parameter(
            name='cp2r',
            displayName='Control point 2 r',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        cp1s = arcpy.Parameter(
            name='cp2s',
            displayName='Control point 2 s',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        params = [cp1r, cp1s, cp2r, cp2s]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        cp1r = parameters[0].value
        cp1s = parameters[1].value
        cp2r = parameters[2].value
        cp2s = parameters[3].value

        suffix = 'linear-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.linear_transform(img_a, (cp1r, cp1s), (cp2r, cp2s))

        save_map_array(suffix, new_img_a, img)
        return

class LogTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Log transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        suffix = 'log-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.log_transform(img_a)

        save_map_array(suffix, new_img_a, img)
        return

class InverseLogTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Inverse-log transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        suffix = 'inverselog-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.inverse_log_transform(img_a)

        save_map_array(suffix, new_img_a, img)
        return

class PowerTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Power transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        gamma = arcpy.Parameter(
            name='gamma',
            displayName='gamma',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        params = [gamma]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        gamma = parameters[0].value

        suffix = 'power-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.power_transform(img_a, gamma)

        save_map_array(suffix, new_img_a, img)
        return

class GrayLevelSlice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gray-level slice"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        lower_gray = arcpy.Parameter(
            name='lower_gray',
            displayName='Lower gray level',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        upper_gray = arcpy.Parameter(
            name='upper_gray',
            displayName='Upper gray level',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        new_gray = arcpy.Parameter(
            name='new_gray',
            displayName='New gray level',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        binary = arcpy.Parameter(
            name='binary',
            displayName='Binary',
            direction='Input',
            datatype='GPBoolean',
            parameterType='Required',
        )
        params = [lower_gray, upper_gray, new_gray, binary]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        lower_gray = parameters[0].value
        upper_gray = parameters[1].value
        new_gray = parameters[2].value
        binary = parameters[3].value

        suffix = 'slice-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.gray_level_slice(img_a, (lower_gray, upper_gray), new_gray, binary)

        save_map_array(suffix, new_img_a, img)
        return

class BitPlaneSlice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bit-plane slice"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        bit_plane = arcpy.Parameter(
            name='bit_plane',
            displayName='Bit plane',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        params = [bit_plane]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        bit_plane = parameters[0].value

        suffix = 'bitplane-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.bit_plane_slice(img_a, bit_plane)

        save_map_array(suffix, new_img_a, img)
        return

class HistogramEqualize(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram equalize"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        suffix = 'histogram-' + str(int(time.time()))
        img, img_a = get_map(suffix)

        new_img_a = dippy.histogram_equalize(img_a)

        save_map_array(suffix, new_img_a, img)
        return
