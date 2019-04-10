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
import numpy as np

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
            RescaleGrayLevels,
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
            Convolute,
            WeightedAverage,
            FirstDerivative,
            SecondDerivative,
            Sharpen,
            HighBoostFilter,
            ConvertRGBToCMY,
            ConvertCMYToRGB,
            ConvertRGBToHSI,
            ConvertHSIToRGB,
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

def add_data_to_map(data_path):
    arcpy.mp.ArcGISProject('CURRENT').activeMap.addDataFromPath(data_path)

def convert_raster_to_grayscale(raster_path, gray_raster_path):
    raster = arcpy.Raster(raster_path)
    if os.path.exists(gray_raster_path):
        os.remove(gray_raster_path)
    raster_a = agpy.convert_raster_to_grayscale_array(raster)
    raster = agpy.convert_array_to_raster(raster_a, raster)
    arcpy.CopyRaster_management(raster, gray_raster_path)

def save_raster_array(raster_path, raster_a, ref_raster):
    if os.path.exists(raster_path):
        os.remove(raster_path)
    raster = agpy.convert_array_to_raster(raster_a, ref_raster)
    arcpy.CopyRaster_management(raster, raster_path)

def get_raster_array(raster_layer):
    raster = arcpy.Raster(raster_layer)
    if raster.isInteger:
        raster_a = arcpy.RasterToNumPyArray(raster)
    else:
        raster_a = arcpy.RasterToNumPyArray(raster, nodata_to_value=np.nan)
    return raster, raster_a

def save_add_raster_array(raster_path, raster_a, ref_raster):
    save_raster_array(raster_path, raster_a, ref_raster)
    add_data_to_map(raster_path)

class ConvertMapToGrayscale(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert the map to grayscale"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        color_output_tiff = arcpy.Parameter(
            name='color_output_tiff',
            displayName='Color Output TIFF',
            direction='Output',
            datatype='DEFile',
            parameterType='Required',
        )
        grayscale_output_raster = arcpy.Parameter(
            name='grayscale_output_raster',
            displayName='Grayscale Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [color_output_tiff, grayscale_output_raster]
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
        color_output_tiff = parameters[0].valueAsText
        grayscale_output_raster = parameters[1].valueAsText

        export_map(color_output_tiff)
        convert_raster_to_grayscale(color_output_tiff, grayscale_output_raster)
        return

class RescaleGrayLevels(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rescale gray levels"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer with 256 Gray Levels',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        lower_gray_level = arcpy.Parameter(
            name='lower_gray_level',
            displayName='Lower Gray Level',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        upper_gray_level = arcpy.Parameter(
            name='upper_gray_level',
            displayName='Upper Gray Level',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, lower_gray_level, upper_gray_level, output_raster]
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
        lower_gray_level = parameters[1].value
        upper_gray_level = parameters[2].value
        output_raster = parameters[3].valueAsText

        img, img_a = get_raster_array(raster_layer)
        gray_level_range = (lower_gray_level, upper_gray_level)
        new_img_a = dippy.rescale_gray_levels(img_a, gray_level_range)
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, gray_levels, output_raster]
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
        output_raster = parameters[2].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.grayscale_transform(img_a, gray_levels)
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, output_raster]
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
        output_raster = parameters[1].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.negative_transform(img_a)
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, cp1r, cp1s, cp2r, cp2s, output_raster]
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
        output_raster = parameters[5].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.linear_transform(img_a, (cp1r, cp1s), (cp2r, cp2s))
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, output_raster]
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
        output_raster = parameters[1].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.log_transform(img_a)
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, output_raster]
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
        output_raster = parameters[1].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.inverse_log_transform(img_a)
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, gamma, output_raster]
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
        output_raster = parameters[2].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.power_transform(img_a, gamma)
        save_raster_array(output_raster, new_img_a, img)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, lower_gray, upper_gray, new_gray, binary, output_raster]
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
        output_raster = parameters[5].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.gray_level_slice(img_a, (lower_gray, upper_gray), new_gray, binary)
        save_raster_array(output_raster, new_img_a, img)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, bit_plane, output_raster]
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
        output_raster = parameters[2].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.bit_plane_slice(img_a, bit_plane)
        save_raster_array(output_raster, new_img_a, img)
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
            displayName='Raster Layer with 256 Gray Levels',
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
            displayName='Raster Layer with 256 Gray Levels',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, output_raster]
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
        output_raster = parameters[1].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.histogram_equalize(img_a)
        save_raster_array(output_raster, new_img_a, img)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, output_raster]
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
        output_raster = parameters[1].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.bitwise_not(img_a)
        save_raster_array(output_raster, new_img_a, img)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_raster]
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
        output_raster = parameters[2].valueAsText

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.bitwise_and(img1_a, img2_a)
        save_add_raster_array(output_raster, new_img_a, img1)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_raster]
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
        output_raster = parameters[2].valueAsText

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.bitwise_or(img1_a, img2_a)
        save_add_raster_array(output_raster, new_img_a, img1)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_raster]
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
        output_raster = parameters[2].valueAsText

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.bitwise_xor(img1_a, img2_a)
        save_add_raster_array(output_raster, new_img_a, img1)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layers, output_raster]
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
        output_raster = parameters[1].valueAsText

        imgs_a = []
        for raster_layer in raster_layers.split(';'):
            img, img_a = get_raster_array(raster_layer)
            imgs_a.append(img_a)
        new_img_a = dippy.add(imgs_a)
        save_add_raster_array(output_raster, new_img_a, img)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer_1, raster_layer_2, output_raster]
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
        output_raster = parameters[2].valueAsText

        img1, img1_a = get_raster_array(raster_layer_1)
        img2, img2_a = get_raster_array(raster_layer_2)
        new_img_a = dippy.subtract(img1_a, img2_a)
        save_add_raster_array(output_raster, new_img_a, img1)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, prob, max, output_raster]
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
        output_raster = parameters[3].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.add_noise(img_a, prob, max)
        save_add_raster_array(output_raster, new_img_a, img)
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layers, output_raster]
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
        output_raster = parameters[1].valueAsText

        imgs_a = []
        for raster_layer in raster_layers.split(';'):
            img, img_a = get_raster_array(raster_layer)
            imgs_a.append(img_a)
        new_img_a = dippy.average(imgs_a)
        save_add_raster_array(output_raster, new_img_a, img)
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
        local_mean_raster = arcpy.Parameter(
            name='local_mean_raster',
            displayName='Output Local Mean Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Optional',
        )
        local_std_raster = arcpy.Parameter(
            name='local_std_raster',
            displayName='Output Local Standard Deviation Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Optional',
        )
        local_median_raster = arcpy.Parameter(
            name='local_median_raster',
            displayName='Output Local Median Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Optional',
        )
        local_min_raster = arcpy.Parameter(
            name='local_min_raster',
            displayName='Output Local Minimum Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Optional',
        )
        local_max_raster = arcpy.Parameter(
            name='local_max_raster',
            displayName='Output Local Maximum Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Optional',
        )
        params = [raster_layer, width, height, local_mean_raster, local_std_raster, local_median_raster, local_min_raster, local_max_raster]
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
        local_mean_raster = parameters[3].valueAsText
        local_std_raster = parameters[4].valueAsText
        local_median_raster = parameters[5].valueAsText
        local_min_raster = parameters[6].valueAsText
        local_max_raster = parameters[7].valueAsText

        stats = 0
        if local_mean_raster:
            stats |= 1
        if local_std_raster:
            stats |= 2
        if local_median_raster:
            stats |= 4
        if local_min_raster:
            stats |= 8
        if local_max_raster:
            stats |= 16

        if stats:
            img, img_a = get_raster_array(raster_layer)
            local_a = dippy.local_statistics(img_a, (height, width), stats)
            if local_mean_raster:
                save_raster_array(local_mean_raster, local_a[0], img)
            if local_std_raster:
                save_raster_array(local_std_raster, local_a[1], img)
            if local_median_raster:
                save_raster_array(local_median_raster, local_a[2], img)
            if local_min_raster:
                save_raster_array(local_min_raster, local_a[3], img)
            if local_max_raster:
                save_raster_array(local_max_raster, local_a[4], img)
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
            displayName='Gray-level Multiplier',
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, local_mean_layer, local_std_layer, multi, k0, k1, k2, output_raster]
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
        output_raster = parameters[7].valueAsText

        img, img_a = get_raster_array(raster_layer)
        local_mean_a = get_raster_array(local_mean_layer)[1]
        local_std_a = get_raster_array(local_std_layer)[1]
        new_img_a = dippy.local_enhance(img_a, local_mean_a, local_std_a, multi, (k0, k1, k2))
        save_add_raster_array(output_raster, new_img_a, img)
        return

class Convolute(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convolute"
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
        mask = arcpy.Parameter(
            name='mask',
            displayName='Filter mask values separated by a space',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        average = arcpy.Parameter(
            name='average',
            displayName='Weighted Average',
            direction='Input',
            datatype='GPBoolean',
            parameterType='Optional',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, width, height, mask, average, output_raster]
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
        mask = parameters[3].value
        average = parameters[4].value
        output_raster = parameters[5].valueAsText

        img, img_a = get_raster_array(raster_layer)
        mask = np.fromstring(mask, sep=' ').reshape(height, width)
        new_img_a = dippy.convolute(img_a, mask, average)
        save_add_raster_array(output_raster, new_img_a, img)
        return

class WeightedAverage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Weighted average"
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
        weights = arcpy.Parameter(
            name='weights',
            displayName='Weights separated by a space',
            direction='Input',
            datatype='GPString',
            parameterType='Required',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, width, height, weights, output_raster]
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
        weights = parameters[3].value
        output_raster = parameters[4].valueAsText

        img, img_a = get_raster_array(raster_layer)
        weights = np.fromstring(weights, sep=' ').reshape(height, width)
        new_img_a = dippy.weighted_average(img_a, weights)
        save_add_raster_array(output_raster, new_img_a, img)
        return

class FirstDerivative(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "First derivative using the Sobel filter"
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
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, output_raster]
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
        output_raster = parameters[1].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.first_derivative(img_a)
        save_add_raster_array(output_raster, new_img_a, img)
        return

class SecondDerivative(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Second derivative using the Laplacian filter"
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
        diag = arcpy.Parameter(
            name='diag',
            displayName='Consider Diagonal Directions',
            direction='Input',
            datatype='GPBoolean',
            parameterType='Optional',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, diag, output_raster]
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
        diag = parameters[1].value
        output_raster = parameters[2].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.second_derivative(img_a, diag)
        save_add_raster_array(output_raster, new_img_a, img)
        return

class Sharpen(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sharpen using the Laplacian"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        raster_layer = arcpy.Parameter(
            name='raster_layer',
            displayName='Raster Layer with 256 Gray Levels',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        diag = arcpy.Parameter(
            name='diag',
            displayName='Consider Diagonal Directions',
            direction='Input',
            datatype='GPBoolean',
            parameterType='Optional',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, diag, output_raster]
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
        diag = parameters[1].value
        output_raster = parameters[2].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.sharpen(img_a, diag)
        save_add_raster_array(output_raster, new_img_a, img)
        return

class HighBoostFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "High-boost filter using the Laplacian"
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
        A = arcpy.Parameter(
            name='A',
            displayName='High-boost filter parameter',
            direction='Input',
            datatype='GPDouble',
            parameterType='Required',
        )
        diag = arcpy.Parameter(
            name='diag',
            displayName='Consider Diagonal Directions',
            direction='Input',
            datatype='GPBoolean',
            parameterType='Optional',
        )
        output_raster = arcpy.Parameter(
            name='output_raster',
            displayName='Output Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [raster_layer, A, diag, output_raster]
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
        A = parameters[1].value
        diag = parameters[2].value
        output_raster = parameters[3].valueAsText

        img, img_a = get_raster_array(raster_layer)
        new_img_a = dippy.high_boost_filter(img_a, A, diag)
        save_add_raster_array(output_raster, new_img_a, img)
        return

class ConvertRGBToCMY(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert RGB to CMY"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        R_layer = arcpy.Parameter(
            name='R_layer',
            displayName='R Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        G_layer = arcpy.Parameter(
            name='G_layer',
            displayName='G Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        B_layer = arcpy.Parameter(
            name='B_layer',
            displayName='B Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        C_raster = arcpy.Parameter(
            name='C_raster',
            displayName='Output C Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        M_raster = arcpy.Parameter(
            name='M_raster',
            displayName='Output M Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        Y_raster = arcpy.Parameter(
            name='Y_raster',
            displayName='Output Y Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [R_layer, G_layer, B_layer, C_raster, M_raster, Y_raster]
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
        R_layer = parameters[0].valueAsText
        G_layer = parameters[1].valueAsText
        B_layer = parameters[2].valueAsText
        C_raster = parameters[3].valueAsText
        M_raster = parameters[4].valueAsText
        Y_raster = parameters[5].valueAsText

        R, R_a = get_raster_array(R_layer)
        G, G_a = get_raster_array(G_layer)
        B, B_a = get_raster_array(B_layer)
        C_a, M_a, Y_a = dippy.rgb_to_cmy(R_a, G_a, B_a)
        save_raster_array(C_raster, C_a, R)
        save_raster_array(M_raster, M_a, G)
        save_raster_array(Y_raster, Y_a, B)
        return

class ConvertCMYToRGB(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert CMY to RGB"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        C_layer = arcpy.Parameter(
            name='C_layer',
            displayName='C Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        M_layer = arcpy.Parameter(
            name='M_layer',
            displayName='M Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        Y_layer = arcpy.Parameter(
            name='Y_layer',
            displayName='Y Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        R_raster = arcpy.Parameter(
            name='R_raster',
            displayName='Output R Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        G_raster = arcpy.Parameter(
            name='G_raster',
            displayName='Output G Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        B_raster = arcpy.Parameter(
            name='B_raster',
            displayName='Output B Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [C_layer, M_layer, Y_layer, R_raster, G_raster, B_raster]
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
        C_layer = parameters[0].valueAsText
        M_layer = parameters[1].valueAsText
        Y_layer = parameters[2].valueAsText
        R_raster = parameters[3].valueAsText
        G_raster = parameters[4].valueAsText
        B_raster = parameters[5].valueAsText

        C, C_a = get_raster_array(C_layer)
        M, M_a = get_raster_array(M_layer)
        Y, Y_a = get_raster_array(Y_layer)
        R_a, G_a, B_a = dippy.cmy_to_rgb(C_a, M_a, Y_a)
        save_raster_array(R_raster, R_a, C)
        save_raster_array(G_raster, G_a, M)
        save_raster_array(B_raster, B_a, Y)
        return

class ConvertRGBToHSI(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert RGB to HSI"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        R_layer = arcpy.Parameter(
            name='R_layer',
            displayName='R Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        G_layer = arcpy.Parameter(
            name='G_layer',
            displayName='G Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        B_layer = arcpy.Parameter(
            name='B_layer',
            displayName='B Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        H_raster = arcpy.Parameter(
            name='H_raster',
            displayName='Output H Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        S_raster = arcpy.Parameter(
            name='S_raster',
            displayName='Output S Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        I_raster = arcpy.Parameter(
            name='I_raster',
            displayName='Output I Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [R_layer, G_layer, B_layer, H_raster, S_raster, I_raster]
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
        R_layer = parameters[0].valueAsText
        G_layer = parameters[1].valueAsText
        B_layer = parameters[2].valueAsText
        H_raster = parameters[3].valueAsText
        S_raster = parameters[4].valueAsText
        I_raster = parameters[5].valueAsText

        R, R_a = get_raster_array(R_layer)
        G, G_a = get_raster_array(G_layer)
        B, B_a = get_raster_array(B_layer)
        H_a, S_a, I_a = dippy.rgb_to_hsi(R_a, G_a, B_a)
        save_raster_array(H_raster, H_a, R)
        save_raster_array(S_raster, S_a, G)
        save_raster_array(I_raster, I_a, B)
        return

class ConvertHSIToRGB(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert HSI to RGB"
        self.description = self.label
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        H_layer = arcpy.Parameter(
            name='H_layer',
            displayName='H Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        S_layer = arcpy.Parameter(
            name='S_layer',
            displayName='S Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        I_layer = arcpy.Parameter(
            name='I_layer',
            displayName='I Layer',
            direction='Input',
            datatype='GPRasterLayer',
            parameterType='Required',
        )
        R_raster = arcpy.Parameter(
            name='R_raster',
            displayName='Output R Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        G_raster = arcpy.Parameter(
            name='G_raster',
            displayName='Output G Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        B_raster = arcpy.Parameter(
            name='B_raster',
            displayName='Output B Raster',
            direction='Output',
            datatype='DERasterDataset',
            parameterType='Required',
        )
        params = [H_layer, S_layer, I_layer, R_raster, G_raster, B_raster]
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
        H_layer = parameters[0].valueAsText
        S_layer = parameters[1].valueAsText
        I_layer = parameters[2].valueAsText
        R_raster = parameters[3].valueAsText
        G_raster = parameters[4].valueAsText
        B_raster = parameters[5].valueAsText

        H, H_a = get_raster_array(H_layer)
        S, S_a = get_raster_array(S_layer)
        I, I_a = get_raster_array(I_layer)
        R_a, G_a, B_a = dippy.hsi_to_rgb(H_a, S_a, I_a)
        save_raster_array(R_raster, R_a, H)
        save_raster_array(G_raster, G_a, S)
        save_raster_array(B_raster, B_a, I)
        return
