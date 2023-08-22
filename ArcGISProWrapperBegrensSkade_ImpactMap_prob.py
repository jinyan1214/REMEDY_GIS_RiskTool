coding_guide = 0 #avoids some sort of coding interpretation bugs
# Prepared for open source release August 2022

log_path = r'C:\Users\Jinyan\Documents\ArcGIS\Projects\CampusUlleval2\REMEDY_GIS_RiskTool_JZ\log'
lyr_path = r'C:\Users\Jinyan\Documents\ArcGIS\Projects\CampusUlleval2\REMEDY_GIS_RiskTool_JZ\lyr'

import arcpy
import sys
import os
import traceback
import logging.handlers
import pathlib
sys.path.append(log_path)
import importlib
import Utils
import Utils_arcpy
import BegrensSkade
import BegrensSkadeLib
importlib.reload(Utils)
importlib.reload(Utils_arcpy)
importlib.reload(BegrensSkade)
importlib.reload(BegrensSkadeLib)

CALCULATION_RANGE = 380

##############  SETUP LOGGERS ##############################
maxLoggerFileSize = 2 * 1024 * 1024
logger = logging.getLogger('BegrensSkade_IMPACTMAP')
if not len(logger.handlers):
    logFile = log_path + '//BegrensSkadeII_ArcGISPro_IMPACTMAP.log'
    hdlr = logging.handlers.RotatingFileHandler(logFile, 'a', maxLoggerFileSize, 20)
    formatter = logging.Formatter('%(asctime)s %(levelname)s Thread %(thread)d %(message)s ')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)
############################################################

##############  READ PARAMETERS ############################
excavation_polys_fl = arcpy.GetParameter(0)
output_folder = arcpy.GetParameterAsText(1)
output_name = arcpy.GetParameterAsText(2)
output_resolution = arcpy.GetParameterAsText(3)
coord_syst = arcpy.GetParameterAsText(4)

logger.info("------------------------------NEW RUN: " + str(output_name)+"-------------------------------")

sr = arcpy.SpatialReference()
sr.loadFromString(coord_syst)
output_proj = sr.PCSCode

bShortterm = arcpy.GetParameter(5)
if bShortterm:
    excavation_depth = arcpy.GetParameter(6)
    short_term_curve = arcpy.GetParameterAsText(7)
else:
    excavation_depth = None
    short_term_curve = None

dtb_raster = arcpy.GetParameter(8)
pw_reduction_curve = arcpy.GetParameterAsText(9)
porewp_red = arcpy.GetParameter(10)
dry_crust_thk = arcpy.GetParameter(11)
dep_groundwater = arcpy.GetParameter(12)
density_sat = arcpy.GetParameter(13)
OCR = arcpy.GetParameter(14)
janbu_ref_stress = arcpy.GetParameter(15)
janbu_const = arcpy.GetParameter(16)
janbu_m = arcpy.GetParameter(17)
consolidation_time = arcpy.GetParameter(18)
dvmax_stats = arcpy.GetParameter(19)
dvmax_mean = dvmax_stats.getTrueValue(0,0)
dvmax_cv = dvmax_stats.getTrueValue(0,1)
dvmax_range = dvmax_stats.getTrueValue(0,2)
dvmax_nugget_ratio = dvmax_stats.getTrueValue(0,3)
eta_stats = arcpy.GetParameter(20)
eta_mean = eta_stats.getTrueValue(0,0)
eta_cv = eta_stats.getTrueValue(0,1)
eta_range = eta_stats.getTrueValue(0,2)
eta_nugget_ratio = eta_stats.getTrueValue(0,3)
outQuantile = arcpy.GetParameter(20)

#bContours = arcpy.GetParameter(19)
#if bContours:
#    contour_interval = arcpy.GetParameter(20)

##############  SET PROJECTION ###########################
exc_proj = Utils_arcpy.getProjCodeFromFC(excavation_polys_fl)

######  PROJECT BUILDING AND EXCAVATION TO OUTPUT  ########
excavation_polys_projected = False
if exc_proj != output_proj:
    arcpy.AddMessage("Projecting excavation polygon..")
    excavation_polys_projected = output_folder + os.sep + "exc_proj.shp"
    arcpy.Project_management(excavation_polys_fl, excavation_polys_projected, output_proj)
    excavation_polys_fl = excavation_polys_projected

################ GET EXCAVATION INFO #####################
excavation_outline_as_json = Utils_arcpy.getConstructionAsJson(excavation_polys_fl)
excavation_desc = arcpy.Describe(excavation_polys_fl)

############### HANDELING OF INPUT RASTER ################
# If necessary, projects raster to the working projection
raster_desc = arcpy.Describe(dtb_raster)
dtb_raster_proj = raster_desc.SpatialReference.PCSCode
dtb_proj_raster = False
if (str(dtb_raster_proj) != str(output_proj)):
    logger.info("START raster projection")
    arcpy.AddMessage("Projecting raster...")
    dtb_proj_raster = "proj_raster"
    if os.path.exists(dtb_proj_raster):
        os.remove(dtb_proj_raster)
    arcpy.ProjectRaster_management(dtb_raster, dtb_proj_raster, output_proj)
    dtb_raster = dtb_proj_raster
    logger.info("DONE raster projection")

#Checks if raster area and requested resolution is too demanding and if clipping is necessary
for row in arcpy.da.SearchCursor(excavation_polys_fl, ['SHAPE@']):
    extent = row[0].extent
xmin = extent.XMin - CALCULATION_RANGE
xmax = extent.XMax + CALCULATION_RANGE
ymin = extent.YMin - CALCULATION_RANGE
ymax = extent.YMax + CALCULATION_RANGE
extent_str = str(xmin)+ " " + str(ymin) + " " + str(xmax) + " " + str(ymax)
area = abs(ymax-ymin)*abs(xmax-xmin)
dtb_clip_raster = False
if float(output_resolution)/area < 10/820000:
    logger.debug("START raster clipping")
    arcpy.AddWarning("High output resolution and/or large raster - clipping raster!")
    dtb_clip_raster = "clip_raster"
    arcpy.management.Clip(dtb_raster, extent_str, dtb_clip_raster)
    dtb_raster = dtb_clip_raster
    logger.debug("DONE raster clipping")

#Resampels raster to the specified resolution
dtb_raster_resample = "resampl_raster"
res = str(output_resolution) + " " + str(output_resolution)
logger.debug("START raster resampling")
arcpy.Resample_management(dtb_raster, dtb_raster_resample, res, "NEAREST")
dtb_raster = dtb_raster_resample
logger.debug("DONE raster resampling")
n_cols = arcpy.management.GetRasterProperties(dtb_raster, "COLUMNCOUNT")
n_rows = arcpy.management.GetRasterProperties(dtb_raster, "ROWCOUNT")
logger.info("Dtb raster cols and rows after resampling: " + str(n_cols) + ", " + str(n_rows) )

#Create a tif file from the raster. Necessary for input to GDAL.
raster_desc = arcpy.Describe(dtb_raster)
if raster_desc.extension != ".tif":
    logger.info("START raster to TIFF conversion")
    arcpy.AddMessage("Converting raster...")
    dtb_raster_tiff = output_folder + os.sep + raster_desc.name + ".tif"
    #Delete existing rasters with the same name
    if os.path.exists(dtb_raster_tiff):
        os.remove(dtb_raster_tiff)
    arcpy.RasterToOtherFormat_conversion(raster_desc.name, output_folder,"TIFF")
    dtb_raster_str = dtb_raster_tiff
    logger.info("DONE raster to TIFF conversion")
else:
    dtb_raster_str = output_folder + os.sep + raster_desc.name + ".tif"

############  RUN BEGRENS SKADE CORE FUNCTIONS   ##############
arcpy.AddMessage("Running mainBegrensSkade_ImpactMap...")
try:
    logger.info("Before run impactMap")
    outputFiles = BegrensSkade.mainBegrensSkade_ImpactMap(
        logger,
        excavation_outline_as_json,
        output_folder,
        output_name,
        CALCULATION_RANGE,
        output_proj,
        dtb_raster_str,
        pw_reduction_curve,
        dry_crust_thk,
        dep_groundwater,
        density_sat,
        OCR,
        porewp_red,
        janbu_ref_stress,
        janbu_const,
        janbu_m,
        consolidation_time,
        bShortterm,
        excavation_depth=excavation_depth,
        short_term_curve=short_term_curve,
        dvmax_mean=dvmax_mean,
        dvmax_cv=dvmax_cv,
        dvmax_range = dvmax_range,
        dvmax_nugget_ratio=dvmax_nugget_ratio,
        eta_mean=eta_mean,
        eta_cv=eta_cv,
        eta_range=eta_range,
        eta_nugget_ratio = eta_nugget_ratio,
        outQuantile = outQuantile
        )
except Exception:
    # Print original traceback info
    arcpy.AddError("UNEXPECTED ERROR:\n" + traceback.format_exc())
    arcpy.AddError(sys.exc_info()[1])
    sys.exit()

#################### HANDLE THE RESULT #######################
result_raster = output_folder + os.sep + output_name + ".tif"

if excavation_polys_projected:
    arcpy.Delete_management(excavation_polys_projected)

if dtb_proj_raster:
    try:
        arcpy.Delete_management(output_folder + os.sep + dtb_proj_raster + ".tif")
    except:
        arcpy.AddWarning("Failed to delete temporary raster (projected) ")
if dtb_clip_raster:
    try:
        arcpy.Delete_management(output_folder + os.sep + dtb_clip_raster + ".tif")
    except:
        arcpy.AddWarning("Failed to delete temporary raster (clip) ")
if dtb_raster_resample:
    try:
        arcpy.Delete_management(output_folder + os.sep + dtb_raster_resample + ".tif")
    except:
        arcpy.AddWarning("Failed to delete temporary raster (resample) ")

#Parameter 19: Beregn kotelinjer - bContours
#Parameter 20: Ekvidistanse - contour_interval
#if bContours:
#    if not CheckOutLicense("Spatial", "Spatial Analyst"):
#        if not CheckOutLicense("3D", "3D Analyst"):
#            arcpy.AddError("No Spatial Analyst or 3D Analyst Extension available - terminating")
#            quit()
#        else:
#            extUsed = "3D"
#    else:
#        extUsed = "Spatial"
#    arcpy.AddMessage("...Extension checked out")
#    result_contours = output_folder + os.sep + "contours.shp"
#    arcpy.sa.Contour(result_raster, result_contours, contour_interval, base_contour = 0,z_factor = 1, contour_type = 'CONTOUR')

arcpy.AddMessage("Adding symbology layer to map...")
p = arcpy.mp.ArcGISProject('CURRENT')
pMap = p.activeMap
#newLyr = pMap.addDataFromPath(result_contours)
lyr_impactmap = lyr_path + os.sep + "SETTLEMENT_IMPACT_FIELD.lyrx"
result_lyr = Utils_arcpy.addLayer(pMap, result_raster, lyr_impactmap, output_name)

logger.info("------------------------------DONE-------------------------------\n ")

