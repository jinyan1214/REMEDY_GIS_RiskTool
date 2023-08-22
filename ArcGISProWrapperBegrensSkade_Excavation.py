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
logger = logging.getLogger("BegrensSkade_EXCAVATION")
if not len(logger.handlers):
    logFile = log_path + "//BegrensSkadeII_ArcGISPro_EXCAVATION.log"
    hdlr = logging.handlers.RotatingFileHandler(logFile, "a", maxLoggerFileSize, 20)
    formatter = logging.Formatter("%(asctime)s %(levelname)s Thread %(thread)d %(message)s ")
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)
############################################################


##############  READ PARAMETERS ############################
building_polys_fl = arcpy.GetParameter(0)
excavation_polys_fl = arcpy.GetParameter(1)
output_folder = arcpy.GetParameterAsText(2)
feature_name = arcpy.GetParameterAsText(3)
coord_syst = arcpy.GetParameterAsText(4)

output_spatial_ref = arcpy.SpatialReference()
output_spatial_ref.loadFromString(coord_syst)
output_proj = output_spatial_ref.PCSCode

corner_name = feature_name + "_CORNER"
wall_name = feature_name + "_WALL"
building_name = feature_name + "_BUILDING"

bShortterm = arcpy.GetParameter(5)
if bShortterm:
    excavation_depth = arcpy.GetParameter(6)
    short_term_curve = arcpy.GetParameterAsText(7)
else:
    excavation_depth = None
    short_term_curve = None

bLongterm = arcpy.GetParameter(9)
if bShortterm == False and bLongterm == False:
    arcpy.AddError("Please choose Short term or Long term settlements, or both")
    sys.exit()

if bLongterm:
    dtb_raster = arcpy.GetParameter(9)
    pw_reduction_curve = arcpy.GetParameterAsText(10)
    porewp_red = arcpy.GetParameter(11)
    dry_crust_thk = arcpy.GetParameter(12)
    dep_groundwater = arcpy.GetParameter(13)
    density_sat = arcpy.GetParameter(14)
    OCR = arcpy.GetParameter(15)
    janbu_ref_stress = arcpy.GetParameter(16)
    janbu_const = arcpy.GetParameter(17)
    janbu_m = arcpy.GetParameter(18)
    consolidation_time = arcpy.GetParameter(19)
else:
    dtb_raster = None
    pw_reduction_curve = None
    dry_crust_thk = None
    dep_groundwater = None
    density_sat = None
    OCR = None
    porewp_red = None
    janbu_ref_stress = None
    janbu_const = None
    janbu_m = None
    consolidation_time = None

bVulnerability = arcpy.GetParameter(21)
logger.debug('bVulnerability is : '+str(bVulnerability))
if bVulnerability:
    fields = arcpy.ListFields(building_polys_fl)
    field_map = {}
    for idx, field in enumerate(fields):
        field_map[field.name] = idx - 2
    vuln_idx_count = 0

    try:
        foundation_field = field_map[arcpy.GetParameterAsText(22)]
        vuln_idx_count += 1
    except:
        foundation_field = None
    try:
        structure_field = field_map[arcpy.GetParameterAsText(23)]
        vuln_idx_count += 1
    except:
        structure_field = None
    try:
        status_field = field_map[arcpy.GetParameterAsText(24)]
        vuln_idx_count += 1
    except:
        status_field = None

    if vuln_idx_count > 0:
        arcpy.AddMessage("Vulnerability enabled")
    else:
        arcpy.AddWarning("No valid vulnerability input - disabling vulnerability!")
        bVulnerability = False
else:
    arcpy.AddMessage("Vulnerability disabled")
    foundation_field = None
    structure_field = None
    status_field = None

#TODO: Get the input for probabilsitc settlement curve
if short_term_curve == "Zhao et al. (2023) Probabilistic":
    bProbCurve = True
else:
    bProbCurve = False
logger.debug('short_term_curve is : '+str(short_term_curve))
logger.debug('bProbCurve is : '+str(bProbCurve))
dvmaxOverHeRF = None
dloverdvRV = None
etaRF = None
MCsampleSize = None
outputQuantile = None
if bProbCurve:
    #dvmax here is the short-term of dvmaxOverHe
    dvamxOverList = str(arcpy.GetParameter(32)).split()
    mean = float(dvamxOverList[0])/100.0
    cv = float(dvamxOverList[1])
    dist = dvamxOverList[2]
    range = float(dvamxOverList[3])
    nug = float(dvamxOverList[4])
    dvmaxOverHeRF = BegrensSkadeLib.RandomFieldmodel(
        mean, cv, dist, 'Gaussian', range, nug
        )
    etaList = str(arcpy.GetParameter(33)).split()
    mean = float(etaList[0])
    cv = float(etaList[1])
    dist = etaList[2]
    range = float(etaList[3])
    nug = float(etaList[4])
    etaRF = BegrensSkadeLib.RandomFieldmodel(
        mean, cv, dist, 'Gaussian', range, nug
        )
    dlOverdvList = str(arcpy.GetParameter(34)).split()
    mean = float(dlOverdvList[0])
    cv = float(dlOverdvList[1])
    dist = dlOverdvList[2]
    dloverdvRV = BegrensSkadeLib.RandomVariablemodel(mean, cv, dist)
    #dloverdvRF = BegrensSkadeLib.RandomFieldmodel(dvmaxNominal, dvmaxStd, dvmaxDist, dvmaxVRM, dvmaxRange, dvmaxNug)
    #etaRF = BegrensSkadeLib.RandomFieldmodel(dvmaxNominal, dvmaxStd, dvmaxDist, dvmaxVRM, dvmaxRange, dvmaxNug)
    MCsampleSize = arcpy.GetParameter(35)
    # MCsampleSize = int(MCsampleSize)
    outputQuantile = arcpy.GetParameter(36)
    logger.debug('etaRF.mean is: '+str(etaRF.mean))
#TODO: Get the input for deterministic SSI
dvmax = None
eta = None
dlmax = None
Es = None
nis = None
buildingFeatureFieldListDeterm = None
bSSI = arcpy.GetParameter(31)
logger.debug('bSSI is : '+str(bSSI))
if bSSI:
    if not bShortterm:
        arcpy.AddError('Soil-Structure Interaction analysis only supports short-term ground movements.')
    if short_term_curve != 'Zhao et al. (2022) Deterministic':
        logger.debug('short_term_curve is: '+str(short_term_curve))
        arcpy.AddError("Timoshenko beam SSI (Deterministic) only supports short-term ground movement of Zhao et al. (2022) Deterministic")
    # try:
    # Get paremeters for ground movement
    gf_params = str(arcpy.GetParameter(8)).split()
    logger.debug('gf_params is : '+str(gf_params))
    dvmaxOverHe = float(gf_params[0])
    logger.debug('dvmaxOverHe is : '+str(dvmaxOverHe))
    eta = float(gf_params[1])
    dlmaxOverdvmax = float(gf_params[2])
    He = arcpy.GetParameter(6)
    dvmax = dvmaxOverHe/100.0*He
    # Get field name for Eb, Es, phi_int, dfoot, bfoot
    fields = arcpy.ListFields(building_polys_fl)
    logger.debug('fields is : '+str(fields))
    field_map = {}
    for idx, field in enumerate(fields):
        field_map[field.name] = idx - 2
    logger.debug('field_map is : '+str(field_map))
    SSI_idx_count = 0
    Eb_field = field_map[arcpy.GetParameterAsText(37)]
    SSI_idx_count += 1
    phi_int_field = field_map[arcpy.GetParameterAsText(38)]
    SSI_idx_count += 1
    dfoot_field = field_map[arcpy.GetParameterAsText(39)]
    SSI_idx_count += 1
    EoverG_field = field_map[arcpy.GetParameterAsText(40)]
    SSI_idx_count += 1  
    qz_field = field_map[arcpy.GetParameterAsText(41)]
    SSI_idx_count += 1 
    Es = arcpy.GetParameter(42)
    SSI_idx_count += 1
    nis = arcpy.GetParameter(43)
    SSI_idx_count += 1 
    buildingFeatureFieldListDeterm = [Eb_field, phi_int_field, dfoot_field, EoverG_field, qz_field]  
    # except:
    #     arcpy.AddWarning("No valid deterministic SSI input - disabling SSI!")

#TODO: Get the input for probabilistic SSI
EsRF = None
buildingFeatureFieldListProb = None
bProbSSI = arcpy.GetParameter(44)
logger.debug('bProbSSI is : '+str(bProbSSI))
if bProbSSI:
    if not bProbCurve:
        arcpy.AddError('Probabilistic Soil-Structure Interaction analysis requires probabilistic short-term ground movement') 
    try:
        # Get paremeters for Es randomrfield
        EsRFList = str(arcpy.GetParameter(51)).split()
        mean = float(EsRFList[0])
        cv = float(EsRFList[1])
        dist = EsRFList[2]
        range = float(EsRFList[3])
        nug = float(EsRFList[4])
        EsRF = BegrensSkadeLib.RandomFieldmodel(mean, cv, dist, 'Gaussian', range, nug)
        
        # Get parameters for Ebmean, Ebstd, phi_int, dfoot, bfoot
        fields = arcpy.ListFields(building_polys_fl)
        field_map = {}
        for idx, field in enumerate(fields):
            field_map[field.name] = idx - 2
        SSI_idx_count = 0
        Ebmean_field = field_map[arcpy.GetParameterAsText(45)]
        SSI_idx_count += 1
        Ebcv_field = field_map[arcpy.GetParameterAsText(46)]
        SSI_idx_count += 1
        phi_int_field = field_map[arcpy.GetParameterAsText(38)]
        SSI_idx_count += 1
        dfoot_field = field_map[arcpy.GetParameterAsText(39)]
        SSI_idx_count += 1
        EoverGmean_field = field_map[arcpy.GetParameterAsText(47)]
        SSI_idx_count += 1
        EoverGcv_field = field_map[arcpy.GetParameterAsText(48)]
        SSI_idx_count += 1
        qzmean_field = field_map[arcpy.GetParameterAsText(49)]
        SSI_idx_count += 1
        qzcv_field = field_map[arcpy.GetParameterAsText(50)]
        SSI_idx_count += 1
        nis = arcpy.GetParameter(43)
        SSI_idx_count += 1 
        buildingFeatureFieldListProb = [Ebmean_field, Ebcv_field,
                                          EoverGmean_field, EoverGcv_field,
                                          qzmean_field, qzcv_field,
                                          phi_int_field, dfoot_field
                                          ]
    except:
        arcpy.AddWarning("No valid probabilistic SSI input - disabling SSI!")
##############  GET INPUT PROJECTIONS ####################
building_spatial_ref = arcpy.Describe(building_polys_fl).spatialReference
excavation_spatial_ref = arcpy.Describe(excavation_polys_fl).spatialReference

###  GET EXCAVATION ANS BUILDINGS ON SAME PROJECTION #####
excavation_polys_matched = False
if  excavation_spatial_ref != building_spatial_ref:
    arcpy.AddMessage("Matching input projections before clip..")
    excavation_polys_matched = output_folder + os.sep + "exc_match.shp"
    arcpy.Project_management(excavation_polys_fl, excavation_polys_matched, building_spatial_ref)
    excavation_polys_fl = excavation_polys_matched

################ GET EXCAVATION INFO #####################
excavation_outline_as_json = Utils_arcpy.getConstructionAsJson(excavation_polys_fl)
buildingsClipExtent = Utils_arcpy.getBuildingsClipExtentFromConstruction(excavation_outline_as_json, CALCULATION_RANGE, building_spatial_ref, logger)

################ EXTRACTING BUILDINGS ##################
buildings_clip = output_folder + os.sep + "buildings_clip.shp"
logger.debug("TIME - Starting extraction of buildings")
Utils_arcpy.extractBuildingsFromFL(building_polys_fl, buildingsClipExtent, buildings_clip, logger)
logger.info("TIME - Done extraction of buildings.")

#arcpy.SelectLayerByAttribute_management(building_polys_fl, "CLEAR_SELECTION")

######  PROJECT BUILDING AND EXCAVATION TO OUTPUT  ########
building_polys_projected = False
excavation_polys_projected = False
buildings_clip_projected = False

if building_spatial_ref != output_spatial_ref:

    arcpy.AddMessage("Projecting bulidings polygon..")
    buildings_clip_projected = output_folder + os.sep + "buil_proj.shp"
    arcpy.Project_management(buildings_clip, buildings_clip_projected, output_proj)
    building_polys_fl = buildings_clip_projected

    arcpy.AddMessage("Projecting excavation polygon..")
    excavation_polys_projected = output_folder + os.sep + "exc_proj.shp"
    arcpy.Project_management(excavation_polys_fl, excavation_polys_projected, output_proj)
    excavation_polys_fl = excavation_polys_projected

    excavation_outline_as_json = Utils_arcpy.getConstructionAsJson(excavation_polys_fl)

else:
    building_polys_fl = buildings_clip



############### HANDELING OF INPUT RASTER ################

dtb_proj_raster = False
if bLongterm:
    # If necessary, projects raster to the working projection
    raster_desc = arcpy.Describe(dtb_raster)
    dtb_raster_proj = raster_desc.SpatialReference.PCSCode
    if (str(dtb_raster_proj) != str(output_proj)):
        logger.info("START raster projection")
        arcpy.AddMessage("Projecting bedrock raster...")
        dtb_proj_raster = "temp_raster"
        if os.path.exists(dtb_proj_raster):
            os.remove(dtb_proj_raster)
        arcpy.ProjectRaster_management(dtb_raster, dtb_proj_raster, output_proj)
        dtb_raster = dtb_proj_raster
    logger.info("DONE raster projection")

    # Create a tif file from the raster. Necessary for input to GDAL.
    raster_desc = arcpy.Describe(dtb_raster)
    if raster_desc.extension != ".tif":
        logger.info("START raster to TIFF conversion")
        arcpy.AddMessage("Converting bedrock raster...")
        dtb_raster_tiff = output_folder + os.sep + raster_desc.name + ".tif"
        #Delete existing rasters with the same name
        if os.path.exists(dtb_raster_tiff):
            os.remove(dtb_raster_tiff)
        arcpy.RasterToOtherFormat_conversion(raster_desc.name, output_folder, "TIFF")
        dtb_raster_str = dtb_raster_tiff
    logger.info("DONE raster to TIFF conversion")
else:
    dtb_raster_str = None

############  RUN BEGRENS SKADE CORE FUNCTIONS   ##############
arcpy.AddMessage("Running mainBegrensSkade_Excavation...")

try:
    outputFiles = BegrensSkade.mainBegrensSkade_Excavation(
        logger,
        building_polys_fl,
        excavation_outline_as_json,
        output_folder,
        feature_name,
        output_proj,
        bShortterm,
        excavation_depth=excavation_depth,
        short_term_curve=short_term_curve,
        bLongterm=bLongterm,
        dtb_raster=dtb_raster_str,
        pw_reduction_curve=pw_reduction_curve,
        dry_crust_thk=dry_crust_thk,
        dep_groundwater=dep_groundwater,
        density_sat=density_sat,
        OCR=OCR,
        porewp_red=porewp_red,
        janbu_ref_stress=janbu_ref_stress,
        janbu_const=janbu_const,
        janbu_m=janbu_m,
        consolidation_time=consolidation_time,
        bVulnerability=bVulnerability,
        fieldNameFoundation=foundation_field,
        fieldNameStructure=structure_field,
        fieldNameStatus=status_field,
        bProbCurve = bProbCurve,
        dvmaxRF = dvmaxOverHeRF,
        dloverdvRV = dloverdvRV,
        etaRF = etaRF,
        MCsampleSize = MCsampleSize,
        outputQuantile = outputQuantile,
        bSSI=bSSI,
        buildingFeatureFieldListDeterm = buildingFeatureFieldListDeterm,
        Es = Es,
        nis = nis,
        bProbSSI=bProbSSI,
        EsRF=EsRF,
        buildingFeatureFieldListProb = buildingFeatureFieldListProb,
        )
except Exception:
    # Print original traceback info
    arcpy.AddError("UNEXPECTED ERROR:\n" + traceback.format_exc())
    arcpy.AddError(sys.exc_info()[1])
    sys.exit()

if dtb_proj_raster:
    try:
        arcpy.Delete_management(output_folder + os.sep + dtb_proj_raster + ".tif")
    except:
        arcpy.AddWarning("Failed to delete temporary bedrock raster!")

arcpy.Delete_management(buildings_clip)
if excavation_polys_matched:
    arcpy.Delete_management(excavation_polys_matched)
if building_polys_projected:
    arcpy.Delete_management(building_polys_projected)
if buildings_clip_projected:
    arcpy.Delete_management(buildings_clip_projected)
if excavation_polys_projected:
    arcpy.Delete_management(excavation_polys_projected)

############################### HANDLE THE RESULT ############################################
buildings_Shapefile_result = outputFiles[0]
walls_Shapefile_result = outputFiles[1]
corners_Shapefile_result = outputFiles[2]

arcpy.AddMessage("Adding symbology layer to map...")
p = arcpy.mp.ArcGISProject("CURRENT")
pMap = p.activeMap

if bVulnerability:
    addRiskAngle = Utils.setBooleanParameter(arcpy.GetParameter(24))
    addRiskSettl = Utils.setBooleanParameter(arcpy.GetParameter(25))
addImpactAngle = Utils.setBooleanParameter(arcpy.GetParameter(26))
addImpactSettl = Utils.setBooleanParameter(arcpy.GetParameter(27))
addWalls = Utils.setBooleanParameter(arcpy.GetParameter(28))
addCorners = Utils.setBooleanParameter(arcpy.GetParameter(29))

lyr_corners = lyr_path + os.sep + "CORNER_SV.lyrx"
lyr_walls = lyr_path + os.sep + "WALL_ANGLE.lyrx"
lyr_building_sv_max = lyr_path + os.sep + "BUILDING_TOTAL_SV_MAX_mm.lyrx"
lyr_building_a_max = lyr_path + os.sep + "BUILDING_TOTAL_ANGLE_MAX.lyrx"
lyr_building_risk_sv = lyr_path + os.sep + "BUILDING_RISK_SV_gdal.lyrx"
lyr_building_risk_a = lyr_path + os.sep + "BUILDING_RISK_ANGLE_gdal.lyrx"
lyr_group = lyr_path + os.sep + "GIBV_RUN_.lyrx"

lyr_group = pMap.addLayer(arcpy.mp.LayerFile(lyr_group), "TOP")[0]
lyr_group.name = feature_name

# if addCorners:
#     Utils_arcpy.addLayerToGroup(pMap, corners_Shapefile_result, lyr_corners, lyr_group)
# if addWalls:
#     Utils_arcpy.addLayerToGroup(pMap, walls_Shapefile_result, lyr_walls, lyr_group)
# if bVulnerability:
#     if addRiskAngle:
#         Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_risk_a,lyr_group)
#     if addRiskSettl:
#         Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_risk_sv, lyr_group)
# if addImpactAngle:
#     Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_a_max, lyr_group)
# if addImpactSettl:
#     Utils_arcpy.addLayerToGroup(pMap, buildings_Shapefile_result, lyr_building_sv_max, lyr_group)

#arcpy.SelectLayerByAttribute_management(buildings_Shapefile_result, "CLEAR_SELECTION")


logger.info("------------------------------DONE-------------------------------")



