# Name: isarithmic.py
# Description: Interpolates a surface from points using kriging.
# Then, generates a polygon feature class to be used in isarithmic mapping. 
# Requirements: Spatial Analyst Extension
# Import system modules

import arcpy
from arcpy import env
from arcpy.sa import *

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
    
try:
    # The local variables to be changed to adjust to your dataset
    ##########################################################################################################################
    
    # Set environment settings
    folder = "E:/GitHub/spatialInterpolation/"
    outWorkspace = folder + "weather.gdb"
    arcpy.env.workspace = outWorkspace
    outCoordinateSystem = arcpy.SpatialReference("WGS 1984 UTM Zone 15N")
    arcpy.env.outputCoordinateSystem = outCoordinateSystem
    arcpy.env.overwriteOutput = True

    # Set input csv file variables
    input_table = "yearly"
    inFeatures = "precipitation"
    x_coords = "Longitude"
    y_coords = "Latitude"
    input_table_srid = arcpy.SpatialReference(4326)

    # output file name to be used for interpolation and the final output shapefile that contains the isarithmic map
    outputFileName = "prec2018"
    # Set kriging (interpolation) variable
    field = "2018_PREC"
    # Set input Clipping feature class for clipping the interpolated raster
    inClipFeatures = "iowa_counties"

    # The number of classes for equal interval classification to generate classified values of polygons for isarithmic mapping:
    numOfClasses = 8
    
    ############################################################################################################################
    
    # Make the XY event layer...
    arcpy.MakeXYEventLayer_management(input_table, x_coords, y_coords, inFeatures, input_table_srid)

    # Print the total rows
    print(arcpy.GetCount_management(inFeatures))
    print("observations successfully read from table to add xy event table with coordinates.")

    inFeaturesDesc = arcpy.Describe(inFeatures)
    projectedInFeatures = input_table + "_points"
    # Determine if the input has a defined coordinate system, can't project it if it does not
    if inFeaturesDesc.spatialReference.Name == "Unknown":
        print('skipped this feature class due to undefined coordinate system: ' + inFeaturesDesc)
    else: 
        # run project tool to change the projection to UTM Zone 15N
        arcpy.Project_management(inFeatures, projectedInFeatures, outCoordinateSystem)
        print("Input point feature layer was successfully reprojected to: " + projectedInFeatures)

    projectedInFeaturesDesc = arcpy.Describe(projectedInFeatures)
    width = abs(projectedInFeaturesDesc.extent.XMax - projectedInFeaturesDesc.extent.XMin)
    height = abs(projectedInFeaturesDesc.extent.YMax - projectedInFeaturesDesc.extent.YMin)
    dim = -1
    if width > height:
        dim = height
    else:
        dim = width
    # Set the cell size
    # It is the shorter of the width or the height of the extent of the input point features,
    # in the input spatial reference, divided by 2000 (as compared to arcpy default 250).
    arcpy.env.cellSize = 1.0 * dim / 2000

    # INTERPOLATION
    majorRange = 2.6
    partialSill = 542
    nugget = 0
    # Set complex variables
    kModelOrdinary = KrigingModelOrdinary("CIRCULAR", majorRange = majorRange, partialSill = partialSill, nugget = nugget)

    # Execute Kriging - no need to save it
    outKriging = Kriging(projectedInFeatures, field, kModelOrdinary)
    #outKriging.save(outputFileName)
    print("Kriging was successfully performed on " + field + " attribute of the table: " + input_table + " Output:", outKriging)

    #clip the kriging raster using Iowa counties feature class
    #First, reproject iowa counties to UTM Zone 15N (output coordinate system)
    projectedClipFeatures = inClipFeatures + "_utm"
    clipFeaturesDesc = arcpy.Describe(inClipFeatures)

    # Determine if the input has a defined coordinate system, can't project it if it does not
    if clipFeaturesDesc.spatialReference.Name == "Unknown":
        print('skipped this feature class due to undefined coordinate system: ' + inClipFeatures)
    else: 
        # run project tool to change the projection to UTM Zone 15N
        arcpy.Project_management(inClipFeatures, projectedClipFeatures, outCoordinateSystem)
        print("Input clip feature class " + inClipFeatures + " was successfully reprojected")

    clippedRaster = outputFileName + "c"
    arcpy.Clip_management(outKriging, "#", out_raster = clippedRaster, in_template_dataset = projectedClipFeatures, clipping_geometry = "ClippingGeometry")
    print("Precipitation raster " + clippedRaster + " was clipped successfully.")
    
    #delete the initial raster after clipping
    arcpy.Delete_management(outKriging)
    print("Original kriging raster was deleted successfully.")

    #Create isorithmic map from the raster
    #First, Convert floating raster to integer raster using reclassify.
    #Get the values of all cells in the raster and create breaks using equal interval
    minValue = 99999999;
    maxValue = -1;

    rstArray = arcpy.RasterToNumPyArray(clippedRaster) # Change rasterFile to numpy array
    rows, cols = rstArray.shape                     # Return the rows, columns
    for rowNum in xrange(rows):                     # Loop through the rows
        for colNum in xrange(cols):                 # Loop through the row's columns
            value = rstArray.item(rowNum, colNum)
            # avoid 0 (nodata) values for the surroundings of the raster
            if value > maxValue and value > 0: 
                maxValue = value
            if value < minValue and value > 0:
                minValue = value

    interval = 1.0* (maxValue - minValue) / numOfClasses
    remapValuesArray = [] # each element consist of start value, end value and new value for constructing RemapRange
    breakValue = minValue

    classRangeLookUp = {}
    for i in range(0, numOfClasses):
        classCode = i+1
        startVal = breakValue
        endVal = breakValue + interval
        remap = [startVal, endVal, classCode]
        classRangeLookUp.update( {classCode : str(round(startVal, 3)) + " - " + str(round(endVal, 3))} )
        remapValuesArray.append(remap)
        breakValue += interval

    remapRange = RemapRange(remapValuesArray)
    print("Breaks:")
    print(remapRange)
    #Spatial Analyst > Reclass > Reclassify the raster into a number of class breaks based on cell (precipitation) values
    outReclassify = Reclassify(clippedRaster, "Value", remapRange, "NODATA")
    print("Floating point raster was successfully reclassified into an integer raster")

    #Second, Conversion Tools > From Raster > Raster to Polygon
    arcpy.RasterToPolygon_conversion(outReclassify, outputFileName, "SIMPLIFY", "VALUE")
    print("Isarithmic polygon file: " + outputFileName + " was successfully created.")

    #add new column to store the range of values for each class in the classification
    newField = "range"
    arcpy.AddField_management(outputFileName, newField, "TEXT", 40)

    #assign values to the new field
    #existing field for classification
    existingField = "gridcode"
    cursor = arcpy.UpdateCursor(outputFileName)
    for row in cursor:
        row.setValue(newField, classRangeLookUp.get(row.getValue(existingField)))
        cursor.updateRow(row)
    print("New field with class range values were successfully added to the output polygon.")

    # USE THE CODE BELOW FOR PREPARING THE DATA FOR WEB MAPPING ONLY############################################################
    # the code below is for exporting the output feature class (isarithmic polygon) to a shapefile to be used in web mapping
    # shapefile will be reprojected to wgs84 and then converted to topojson using mapshaper.org
##    arcpy.FeatureClassToShapefile_conversion(outputFileName, folder)
##    print("UTM isarithmic shapefile was created successfully.")
##    shapeFileName = folder + outputFileName + ".shp"
##    projectedShapeFile = folder + outputFileName + "_wgs84.shp"
##    arcpy.Project_management(shapeFileName, projectedShapeFile, arcpy.SpatialReference(4326))
##    print("UTM isarithmic shapefile was reprojected to wgs84 for web mapping.")
##    arcpy.Delete_management(shapeFileName)
##    print("UTM isarithmic shapefile was deleted successfully.")
    #############################################################################################################################
    
except Exception:
    e = sys.exc_info()[1]
    print(e.args[0])
