coding_guide = 0 #avoids some sort of coding interpretation bugs
# Prepared for open source release August 2022

log_path = r'C:\Users\Jinyan\Documents\ArcGIS\Projects\CampusUlleval2\REMEDY_GIS_RiskTool_JZ\log'
lyr_path = r'C:\Users\Jinyan\Documents\ArcGIS\Projects\CampusUlleval2\REMEDY_GIS_RiskTool_JZ\lyr'
# ASRElib_path = r'C:\Users\Jinyan\Desktop\NGI_GIBV\ASRETimoshenko\ASRE_Timo2\x64\Release\ASRE_Timo2.dll'

import arcpy
import sys
import os
import traceback
import logging.handlers
sys.path.append(log_path)
from osgeo import gdal
from osgeo import ogr
import importlib
import Utils
import Utils_arcpy
import BegrensSkade
import BegrensSkadeLib
import numpy as np
from numpy import linalg as la
importlib.reload(Utils)
importlib.reload(Utils_arcpy)
importlib.reload(BegrensSkade)
importlib.reload(BegrensSkadeLib)

from datetime import datetime

CALCULATION_RANGE = 380

##############  SETUP LOGGERS ##############################
maxLoggerFileSize = 2 * 1024 * 1024
logger = logging.getLogger("BegrensSkade_SSIpreprocessing")
if not len(logger.handlers):
    logFile = log_path + "//BegrensSkadeII_ArcGISPro_SSIpreprocessing.log"
    hdlr = logging.handlers.RotatingFileHandler(logFile, "a", maxLoggerFileSize, 20)
    formatter = logging.Formatter("%(asctime)s %(levelname)s Thread %(thread)d %(message)s ")
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)
############################################################


##############  READ PARAMETERS ############################
building_polys_fl = arcpy.GetParameter(0)
excavation_polys_fl = arcpy.GetParameter(1)
output_ws = arcpy.GetParameterAsText(2)
coord_syst = arcpy.GetParameterAsText(3)
BUILDING_PIXEL_SIZE = arcpy.GetParameter(5)
feature_name = 'SSIpreprocessing'
output_spatial_ref = arcpy.SpatialReference()
output_spatial_ref.loadFromString(coord_syst)
output_proj = output_spatial_ref.PCSCode

##############  GET INPUT PROJECTIONS ####################
building_spatial_ref = arcpy.Describe(building_polys_fl).spatialReference
excavation_spatial_ref = arcpy.Describe(excavation_polys_fl).spatialReference
###  GET EXCAVATION ANS BUILDINGS ON SAME PROJECTION #####
excavation_polys_matched = False
if  excavation_spatial_ref != building_spatial_ref:
    arcpy.AddMessage("Matching input projections before clip..")
    excavation_polys_matched = output_ws + os.sep + "exc_match.shp"
    arcpy.Project_management(excavation_polys_fl, excavation_polys_matched, building_spatial_ref)
    excavation_polys_fl = excavation_polys_matched
################ GET EXCAVATION INFO #####################
excavation_outline_as_json = Utils_arcpy.getConstructionAsJson(excavation_polys_fl)
buildingsClipExtent = Utils_arcpy.getBuildingsClipExtentFromConstruction(excavation_outline_as_json, CALCULATION_RANGE, building_spatial_ref, logger)
################ EXTRACTING BUILDINGS ##################
buildings_clip = output_ws + os.sep + "buildings_clip.shp"
logger.debug("TIME - Starting extraction of buildings")
Utils_arcpy.extractBuildingsFromFL(building_polys_fl, buildingsClipExtent, buildings_clip, logger)
logger.info("TIME - Done extraction of buildings.")
######  PROJECT BUILDING AND EXCAVATION TO OUTPUT  ########
building_polys_projected = False
excavation_polys_projected = False
buildings_clip_projected = False

if building_spatial_ref != output_spatial_ref:

    arcpy.AddMessage("Projecting bulidings polygon..")
    buildings_clip_projected = output_ws + os.sep + "buil_proj.shp"
    arcpy.Project_management(buildings_clip, buildings_clip_projected, output_proj)
    building_polys_fl = buildings_clip_projected

    arcpy.AddMessage("Projecting excavation polygon..")
    excavation_polys_projected = output_ws + os.sep + "exc_proj.shp"
    arcpy.Project_management(excavation_polys_fl, excavation_polys_projected, output_proj)
    excavation_polys_fl = excavation_polys_projected

    excavation_outline_as_json = Utils_arcpy.getConstructionAsJson(excavation_polys_fl)

else:
    building_polys_fl = buildings_clip

############  RUN BEGRENS SKADE CORE FUNCTIONS   ##############
arcpy.AddMessage("Running mainBegrensSkade_Preprocessing...")
logger.debug("Running mainBegrensSkade_Preprocessing...")

parameter_log = open(output_ws + os.sep + "GIBV_Excavation_params_" + feature_name + ".txt", "w")
parameter_log.write("PARAMETER LOG FOR JZ GIBV Excavation \n")
parameter_log.write(feature_name+ "\n")
parameter_log.write("Time of run: " + str(datetime.today()) + "\n")
parameter_log.write("-----------------------------\n")
parameter_log.write("Output coordinate system code: " + str(output_proj) + "\n")

#The buildings are first converted from vector polygons to pixels, using the polygon to raster conversion tool.
logger.debug("Opening shapefile: {}".format(building_polys_fl))
shp_driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = shp_driver.Open(building_polys_fl, 1)
building_polygons = dataSource.GetLayer()
foundation = None
structure = None
status = None


### TODO: Convert polygons to rasters here. Now I just read the raster field to 
building_raster = output_ws + "\\" + "temp_raster"

# bid = 1
# arcpy.AddField_management(building_polygons, "building_ID", "SHORT", field_is_nullable=True, field_is_required=False)
# for bld in building_polygons:
#     bld.setValue('building_ID', bid)
#     bid += 1

# arcpy.PolygonToRaster_conversion( building_polygons, 'phi', building_raster, "CELL_CENTER","", BUILDING_PIXEL_SIZE)
arcpy.conversion.PolygonToRaster(building_polys_fl, 'buildingID', building_raster, "CELL_CENTER", "", BUILDING_PIXEL_SIZE)

#Important that all data is projected into the raster coord syst
poly_proj = Utils_arcpy.getProjCodeFromFC(building_polys_fl)
building_raster_proj = Utils_arcpy.getProjCodeFromFC(building_raster)

#arcpy.AddMessage("poly proj: " + str(poly_proj))
#arcpy.AddMessage("building_raster proj: " + str(building_raster_proj))

desc_building_ras = arcpy.Describe(building_raster)
cellsize_building_raster = desc_building_ras.meanCellWidth
extent_building = desc_building_ras.extent

#arcpy.AddMessage("cellsize_building_raster: " + str(cellsize_building_raster))
#arcpy.AddMessage("extent_building_raster: " + str(extent_building))

xmin_buil = extent_building.XMin
ymin_buil = extent_building.YMin
xmax_buil = extent_building.XMax
ymax_buil = extent_building.YMax

xlen_buil = xmax_buil - xmin_buil
ylen_buil = ymax_buil - ymin_buil

#Convert the raster to a numpy array for further analysis of second area moment.
builRasterArray = arcpy.RasterToNumPyArray(building_raster, nodata_to_value=0)

arcpy.AddMessage(len(builRasterArray[0]))
arcpy.AddMessage(builRasterArray[0])
arcpy.AddMessage("builRasterArray: "+str(builRasterArray))

nrows_buil, ncols_buil = builRasterArray.shape
arcpy.AddMessage("nrows, ncols: " + str(nrows_buil) + ", " + str(ncols_buil))

dx_buil = xlen_buil / ncols_buil
dy_buil = ylen_buil / nrows_buil

arcpy.AddMessage("xmin, xlen   building ras: " + str(round(xmin_buil, 0)) + ", " + str(round(xlen_buil, 0)))
arcpy.AddMessage("ymin, ylen   building ras: " + str(round(ymin_buil, 0)) + ", " + str(round(ylen_buil, 0)))
arcpy.AddMessage("ncols, nrows building ras: " + str(ncols_buil) + "; " + str(nrows_buil))
arcpy.AddMessage("dx, dy       building ras: " + str(round(dx_buil, 1)) + ", " + str(round(dy_buil, 1)))

#Assembling the data structure "building_pixel_dict" that contains the coordinates at of each pixel for a particular building
building_pixel_dict = {}
for col in range(0, ncols_buil - 1):
    for row in range(0, nrows_buil - 1):
        
        auto_BID = int(builRasterArray.item(row, col))  # Building ID is the value of the cell
        # arcpy.AddMessage("auto_BID: " + str(auto_BID))

        #bid_test2 = builRasterArray.bygning_ID
        #arcpy.AddMessage("test bid: " + str(bid_test2))

        x = xmin_buil + dx_buil * (col)  # Coordinates of the pixel
        y = ymin_buil + dy_buil * (nrows_buil - row)
        
        if auto_BID == 0:
            #Outside a building
            continue
        elif auto_BID in building_pixel_dict.keys():
            # arcpy.AddMessage("auto_BID: " + str(builRasterArray.item(0, 0)))
            building_pixel_dict[auto_BID].append([x, y])
        else:
            building_pixel_dict[auto_BID] = [[x, y]]


count = 0

#Create a polyline feature for visualization of the principal axes (direction and magnitude)
arcpy.CreateFeatureclass_management(output_ws, "temp_lines1", "POLYLINE")
dir_lines1 = output_ws + "\\" + "temp_lines1"
#https://pro.arcgis.com/en/pro-app/latest/tool-reference/data-management/add-field.htm
arcpy.management.AddField(dir_lines1, "buildingID","LONG")
arcpy.management.AddField(dir_lines1, "axisID","LONG")
arcpy.management.AddField(dir_lines1, "bfoot","DOUBLE")
field_names = ['SHAPE@', 'buildingID', 'axisID', 'bfoot'] #longer axis has axisID=1, short axis has axisID=2
dir_lines_cursor1 = arcpy.da.InsertCursor(dir_lines1, field_names)

#Create a polyline feature for visualization of the principal axes (direction and magnitude)
arcpy.CreateFeatureclass_management(output_ws, "temp_lines2", "POLYLINE")
dir_lines2 = output_ws + "\\" + "temp_lines2"
arcpy.management.AddField(dir_lines2, "buildingID","LONG")
arcpy.management.AddField(dir_lines2, "axisID","LONG")
arcpy.management.AddField(dir_lines2, "bfoot","DOUBLE")#https://pro.arcgis.com/en/pro-app/latest/tool-reference/data-management/add-field.htm
field_names = ['SHAPE@', 'buildingID', 'axisID', 'bfoot'] #longer axis has axisID=1, short axis has axisID=2
dir_lines_cursor2 = arcpy.da.InsertCursor(dir_lines2, field_names)

# theInsertCur = arcpy.InsertCursor(dir_lines)

arcpy.AddMessage("building_pixel_dict length: " + str(len(building_pixel_dict.items())))

for bid, pixels in building_pixel_dict.items():
    #Area moment components:
    Ix = 0
    Iy = 0

    for pixel in pixels:
        x = pixel[0]
        y = pixel[1]
        Ix += x
        Iy += y

    centeroid_x = Ix / len(pixels)
    centeroid_y = Iy / len(pixels)

    #Second area moment components:
    Ixx = 0
    Iyy = 0
    Ixy = 0 #symmetric matrix, Ixy = Iyx

    for pixel in pixels:
        x = pixel[0]
        y = pixel[1]
        cx = x - centeroid_x
        cy = y - centeroid_y
        Ixx += cx*cx
        Iyy += cy*cy
        Ixy += cx*cy

    #arcpy.AddMessage("BOID, cent_x, cent_y, Ixx, Iyy, Ixy: " + str(bid) +", " + str(centeroid_x) + ", " + str(centeroid_y) + ", "+ str(Ixx) +", "+ str(Iyy) + ", " + str(Ixy))

    I = np.array([[Ixx, Ixy], [Ixy, Iyy]])
    eigenvalues, eigenvectors = la.eig(I)

    dir_x = 0
    dir_y = 0
    dirlen = 0
    if eigenvalues[0]> eigenvalues[1] and eigenvalues[1] > 0.001:
        dir_x1, dir_y1 = eigenvectors[:, 0] #unit vector
        dir_x2, dir_y2 = eigenvectors[:, 1] #unit vector
        dir1 = eigenvectors[:, 0]
        dir2 = eigenvectors[:, 1]
        dirlen = eigenvalues[0]/eigenvalues[1]
    elif eigenvalues[1]> eigenvalues[0] and eigenvalues[0] > 0.001:
        dir_x1, dir_y1 = eigenvectors[:, 1] #unit vector
        dir_x2, dir_y2 = eigenvectors[:, 0] #unit vector
        dir1 = eigenvectors[:, 1]
        dir2 = eigenvectors[:, 0]
        dirlen = eigenvalues[1]/eigenvalues[0]
    else:
        dir_x1 = 1
        dir_y1 = 0
        dir_x2 = 0
        dir_y2 = 1
        dirlen = 0

    #dir1 = [centeroid_x + dir_x*dirlen, centeroid_y + dir_y*dirlen]
    #dir2 = [centeroid_x - dir_x*dirlen, centeroid_y - dir_y*dirlen]

    #loop to find "length" of building along first dominant direction
    step = 0.5
    dist1 = step
    dist2 = step
    cont = True
    security = 0
    #axis 1
    step_rate = 0.5
    try:
        while cont and security < 10000:
            cont = False
            for pixel in pixels:

                x = pixel[0]
                y = pixel[1]

                line1x = centeroid_x + dir_x1*dist1
                line1y = centeroid_y + dir_y1*dist1
                line2x = centeroid_x - dir_x1*dist2
                line2y = centeroid_y - dir_y1*dist2

                if abs(line1x - x) < step_rate*step and abs(line1y - y) < step_rate*step:

                    dist1 += step
                    cont = True
                    arcpy.AddMessage("HELLOOO1")
                    break

                if abs(line2x - x) < step_rate*step and abs(line2y - y) < step_rate*step:
                    dist2 += step
                    cont = True
                    arcpy.AddMessage("HELLOOO2")
                    break

            security += 1

    except:
        dist1 = 0
        dist2 = 0


    axis1Length = dist1+dist2

    #arcpy.AddMessage("Dominanit build len (" + str(bid) + "): "+ str(dominant_buil_len))

    if centeroid_x < 10000 or centeroid_y < 1000000 or centeroid_x > 1000000 or centeroid_y > 10000000:
        arcpy.AddError("Skipping: " + str(dir1[0]) + ", "+ str(dir1[1]) + ", "+ str(dir2[0]) + ", "+ str(dir2[1]))
        continue
    

    arcpy.AddMessage("Lines directions: " + str(dir1[0]) + ", "+ str(dir1[1]) + ", "+ str(dir2[0]) + ", "+ str(dir2[1]))
        # polyline = arcpy.Polyline(arcpy.Array([arcpy.Point(dir1[0], dir1[1]), arcpy.Point(dir2[0], dir2[1])]))
    polyline1 = arcpy.Polyline(arcpy.Array([arcpy.Point(line1x, line1y), arcpy.Point(line2x, line2y)]))

    #loop to find "length" of building along first dominant direction
    step = 0.5
    dist1 = step
    dist2 = step
    cont = True
    security = 0
    step_rate = 0.5
    try:
        while cont and security < 10000:
            cont = False
            for pixel in pixels:

                x = pixel[0]
                y = pixel[1]

                line1x = centeroid_x + dir_x2*dist1
                line1y = centeroid_y + dir_y2*dist1
                line2x = centeroid_x - dir_x2*dist2
                line2y = centeroid_y - dir_y2*dist2

                if abs(line1x - x) < step_rate*step and abs(line1y - y) < step_rate*step:

                    dist1 += step
                    cont = True
                    arcpy.AddMessage("HELLOOO1")
                    break

                if abs(line2x - x) < step_rate*step and abs(line2y - y) < step_rate*step:
                    dist2 += step
                    cont = True
                    arcpy.AddMessage("HELLOOO2")
                    break

            security += 1

    except:
        dist1 = 0
        dist2 = 0

    axis2Length = dist1+dist2

    # building_geom_dict[bid] = [dir1, dir2, dirlen, dominant_buil_len]

    #arcpy.AddMessage("Dominanit build len (" + str(bid) + "): "+ str(dominant_buil_len))
    

    arcpy.AddMessage("Lines directions: " + str(dir1[0]) + ", "+ str(dir1[1]) + ", "+ str(dir2[0]) + ", "+ str(dir2[1]))
    polyline2 = arcpy.Polyline(arcpy.Array([arcpy.Point(line1x, line1y), arcpy.Point(line2x, line2y)]))
    try:
        arcpy.AddMessage("Lines directions: " + str(dir1[0]) + ", "+ str(dir1[1]) + ", "+ str(dir2[0]) + ", "+ str(dir2[1]))
        dir_lines_cursor1.insertRow((polyline1, bid, 1, axis2Length))
        dir_lines_cursor2.insertRow((polyline2, bid, 2, axis1Length))
    except:
        arcpy.AddError("Failed at: " + str(dir1[0]) + ", "+ str(dir1[1])+ ", " + str(dir2[0]) +", "+ str(dir2[1]))


del dir_lines_cursor1
del dir_lines_cursor2


#TODO: Add estimating building elevation from DEM model
#TODO: Add automatic estimation of building stiffness and other properties (perhaps add brails to estimate opening ratio)