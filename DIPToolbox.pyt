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

import arcpy
import os
import time

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
            ConvertMapToGrayscale,
            GrayscaleTransform,
            NegativeTransform,
            LinearTransform,
            LogTransform,
            InverseLogTransform,
            PowerTransform,
            GrayLevelSlice,
            BitPlaneSlice,
            Histogram,
            HistogramEqualize,
            BitwiseNot,
            BitwiseAnd,
            BitwiseOr,
            BitwiseXor,
            Add,
            Subtract,
            AddNoise,
            Average,
            LocalStatistics,
            LocalEnhance,
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

def export_map(tiff_path):
    if os.path.exists(tiff_path):
        os.remove(tiff_path)
    agpy.export_map_to_tiff(tiff_path, 1024, 1024)

def add_tiff_to_map(tiff_path):
    arcpy.mp.ArcGISProject('CURRENT').activeMap.addDataFromPath(tiff_path)

def convert_tiff_to_grayscale(tiff_path, gray_tiff_path):
    img = arcpy.Raster(tiff_path)
    if os.path.exists(gray_tiff_path):
        os.remove(gray_tiff_path)
    img_a = agpy.convert_raster_to_grayscale_array(img)
    img = agpy.convert_array_to_raster(img_a, img)
    arcpy.CopyRaster_management(img, gray_tiff_path)

def save_raster_array(tiff_path, img_a, ref_img):
    if os.path.exists(tiff_path):
        os.remove(tiff_path)
    img = agpy.convert_array_to_raster(img_a, ref_img)
    arcpy.CopyRaster_management(img, tiff_path)

def get_raster_array(raster_layer):
    img = arcpy.Raster(raster_layer)
    img_a = arcpy.RasterToNumPyArray(img)
    return img, img_a

def save_add_raster_array(tiff_path, img_a, ref_img):
    save_raster_array(tiff_path, img_a, ref_img)
    add_tiff_to_map(tiff_path)

class ConvertMapToGrayscale(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert the map to grayscale"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [output_tiff]
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
        output_tiff = parameters[0].value

        map_tiff = output_tiff.replace('.', '-color.')
        export_map(map_tiff)
        convert_tiff_to_grayscale(map_tiff, output_tiff)
        add_tiff_to_map(output_tiff)
        return

class GrayscaleTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Grayscale transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        gray_levels = arcpy.Parameter(
            name='gray_levels',
            displayName='Gray Levels',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, gray_levels, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        gray_levels = parameters[1].value
        output_tiff = parameters[2].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.grayscale_transform(img_a, gray_levels)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class NegativeTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Negative transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        output_tiff = parameters[1].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.negative_transform(img_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class LinearTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Linear transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
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
        cp2s = arcpy.Parameter(
            name='cp2s',
            displayName='Control point 2 s',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, cp1r, cp1s, cp2r, cp2s, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        cp1r = parameters[1].value
        cp1s = parameters[2].value
        cp2r = parameters[3].value
        cp2s = parameters[4].value
        output_tiff = parameters[5].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.linear_transform(img_a, (cp1r, cp1s), (cp2r, cp2s))
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class LogTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Log transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        output_tiff = parameters[1].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.log_transform(img_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class InverseLogTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Inverse-log transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        output_tiff = parameters[1].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.inverse_log_transform(img_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class PowerTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Power transform"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        gamma = arcpy.Parameter(
            name='gamma',
            displayName='gamma',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, gamma, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        gamma = parameters[1].value
        output_tiff = parameters[2].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.power_transform(img_a, gamma)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class GrayLevelSlice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gray-level slice"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
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
            parameterType='Optional',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, lower_gray, upper_gray, new_gray, binary, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        lower_gray = parameters[1].value
        upper_gray = parameters[2].value
        new_gray = parameters[3].value
        binary = parameters[4].value
        output_tiff = parameters[5].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.gray_level_slice(img_a, (lower_gray, upper_gray), new_gray, binary)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class BitPlaneSlice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bit-plane slice"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        bit_plane = arcpy.Parameter(
            name='bit_plane',
            displayName='Bit plane',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, bit_plane, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        bit_plane = parameters[1].value
        output_tiff = parameters[2].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.bit_plane_slice(img_a, bit_plane)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class Histogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        params = [raster_layer]
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
        raster_layer = parameters[0].valueAsText

        img, img_a = get_raster_array(raster_layer)
        dippy.histogram(img_a)
        return

class HistogramEqualize(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram equalize"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        output_tiff = parameters[1].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.histogram_equalize(img_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class BitwiseNot(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bitwise not"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        output_tiff = parameters[1].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.bitwise_not(img_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class BitwiseAnd(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bitwise and"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer_1 = arcpy.Parameter(
            name='raster_layer_1',
            displayName='Raster Layer 1',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        raster_layer_2 = arcpy.Parameter(
            name='raster_layer_2',
            displayName='Raster Layer 2',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_tiff]
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
        raster_layer_1 = parameters[0].valueAsText
        raster_layer_2 = parameters[1].valueAsText
        output_tiff = parameters[2].value

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.bitwise_and(img1_a, img2_a)
        save_add_raster_array(output_tiff, new_img_a, img1)
        return

class BitwiseOr(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bitwise or"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer_1 = arcpy.Parameter(
            name='raster_layer_1',
            displayName='Raster Layer 1',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        raster_layer_2 = arcpy.Parameter(
            name='raster_layer_2',
            displayName='Raster Layer 2',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_tiff]
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
        raster_layer_1 = parameters[0].valueAsText
        raster_layer_2 = parameters[1].valueAsText
        output_tiff = parameters[2].value

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.bitwise_or(img1_a, img2_a)
        save_add_raster_array(output_tiff, new_img_a, img1)
        return

class BitwiseXor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bitwise xor"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer_1 = arcpy.Parameter(
            name='raster_layer_1',
            displayName='Raster Layer 1',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        raster_layer_2 = arcpy.Parameter(
            name='raster_layer_2',
            displayName='Raster Layer 2',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_tiff]
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
        raster_layer_1 = parameters[0].valueAsText
        raster_layer_2 = parameters[1].valueAsText
        output_tiff = parameters[2].value

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.bitwise_xor(img1_a, img2_a)
        save_add_raster_array(output_tiff, new_img_a, img1)
        return

class Add(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layers = arcpy.Parameter(
            name='raster_layers',
            displayName='Raster Layers',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
            multiValue=True,
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layers, output_tiff]
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
        raster_layers = parameters[0].valueAsText
        output_tiff = parameters[1].value

        imgs_a = []
        for raster_layer in raster_layers.split(';'):
            img, img_a = get_raster_array(raster_layer)
            imgs_a.append(img_a)
        new_img_a = dippy.add(imgs_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class Subtract(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Subtract"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer_1 = arcpy.Parameter(
            name='raster_layer_1',
            displayName='Raster Layer 1',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        raster_layer_2 = arcpy.Parameter(
            name='raster_layer_2',
            displayName='Raster Layer 2',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_tiff]
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
        raster_layer_1 = parameters[0].valueAsText
        raster_layer_2 = parameters[1].valueAsText
        output_tiff = parameters[2].value

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.subtract(img1_a, img2_a)
        save_add_raster_array(output_tiff, new_img_a, img1)
        return

class AddNoise(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add noise"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        prob = arcpy.Parameter(
            name='prob',
            displayName='Probability of Noise',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        max = arcpy.Parameter(
            name='max',
            displayName='Max Noise',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, prob, max, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        prob = parameters[1].value
        max = parameters[2].value
        output_tiff = parameters[3].value

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.add_noise(img_a, prob, max)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class Average(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layers = arcpy.Parameter(
            name='raster_layers',
            displayName='Raster Layers',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
            multiValue=True,
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layers, output_tiff]
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
        raster_layers = parameters[0].valueAsText
        output_tiff = parameters[1].value

        imgs_a = []
        for raster_layer in raster_layers.split(';'):
            img, img_a = get_raster_array(raster_layer)
            imgs_a.append(img_a)
        new_img_a = dippy.average(imgs_a)
        save_add_raster_array(output_tiff, new_img_a, img)
        return

class LocalStatistics(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Local statistics"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        width = arcpy.Parameter(
            name='width',
            displayName='Neighborhood Width',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        height = arcpy.Parameter(
            name='height',
            displayName='Neighborhood Height',
            direction='Input',
            datatype='GPLong',
            parameterType='Required',
        )
        local_mean_tiff = arcpy.Parameter(
            name='local_mean_tiff',
            displayName='Output Local Mean TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Optional',
        )
        local_std_tiff = arcpy.Parameter(
            name='local_std_tiff',
            displayName='Output Local Standard Deviation TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Optional',
        )
        local_median_tiff = arcpy.Parameter(
            name='local_median_tiff',
            displayName='Output Local Median TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Optional',
        )
        local_min_tiff = arcpy.Parameter(
            name='local_min_tiff',
            displayName='Output Local Minimum TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Optional',
        )
        local_max_tiff = arcpy.Parameter(
            name='local_max_tiff',
            displayName='Output Local Maximum TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Optional',
        )
        params = [raster_layer, width, height, local_mean_tiff, local_std_tiff, local_median_tiff, local_min_tiff, local_max_tiff]
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
        raster_layer = parameters[0].valueAsText
        width = parameters[1].value
        height = parameters[2].value
        local_mean_tiff = parameters[3].value
        local_std_tiff = parameters[4].value
        local_median_tiff = parameters[5].value
        local_min_tiff = parameters[6].value
        local_max_tiff = parameters[7].value

        stats = 0
        if local_mean_tiff:
            stats |= 1
        if local_std_tiff:
            stats |= 2
        if local_median_tiff:
            stats |= 4
        if local_min_tiff:
            stats |= 8
        if local_max_tiff:
            stats |= 16

        if stats:
            img, img_a = get_raster_array(raster_layer)
            local_a = dippy.local_statistics(img_a, (height, width), stats)
            if local_mean_tiff:
                save_add_raster_array(local_mean_tiff, local_a[0], img)
            if local_std_tiff:
                save_add_raster_array(local_std_tiff, local_a[1], img)
            if local_median_tiff:
                save_add_raster_array(local_median_tiff, local_a[2], img)
            if local_min_tiff:
                save_add_raster_array(local_min_tiff, local_a[3], img)
            if local_max_tiff:
                save_add_raster_array(local_max_tiff, local_a[4], img)
        return

class LocalEnhance(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Local enhance"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        local_mean_layer = arcpy.Parameter(
            name='local_mean_layer',
            displayName='Local Mean Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        local_std_layer = arcpy.Parameter(
            name='local_std_layer',
            displayName='Local Standard Deviation Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        multi = arcpy.Parameter(
            name='multi',
            displayName='Gray Level Multipler',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        k0 = arcpy.Parameter(
            name='k0',
            displayName='Mean Parameter',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        k1 = arcpy.Parameter(
            name='k1',
            displayName='Lower Standard Deviation Parameter',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        k2 = arcpy.Parameter(
            name='k2',
            displayName='Upper Standard Deviation Parameter',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        output_tiff = arcpy.Parameter(
            name='output_tiff',
            displayName='Output TIFF',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        params = [raster_layer, local_mean_layer, local_std_layer, multi, k0, k1, k2, output_tiff]
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
        raster_layer = parameters[0].valueAsText
        local_mean_layer = parameters[1].valueAsText
        local_std_layer = parameters[2].valueAsText
        multi = parameters[3].value
        k0 = parameters[4].value
        k1 = parameters[5].value
        k2 = parameters[6].value
        output_tiff = parameters[7].value

        img, img_a = get_raster_array(raster_layer)
        local_mean_a = get_raster_array(local_mean_layer)[1]
        local_std_a = get_raster_array(local_std_layer)[1]
        new_img_a = dippy.local_enhance(img_a, local_mean_a, local_std_a, multi, (k0, k1, k2))
        save_add_raster_array(output_tiff, new_img_a, img)
        return
