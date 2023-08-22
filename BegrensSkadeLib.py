import numpy as np
from osgeo import gdal
from osgeo import ogr
import logging.handlers
# import Logger
import Utils
import Utils_arcpy
import scipy.stats
from ctypes import *
import _ctypes
import multiprocessing
import arcpy
import sys
from numpy import linalg as la
import os
ASRE_path = r"C:\Users\Jinyan\Desktop\NGI_GIBV\ASRETimoshenko\ASRE_Timo2\x64\Release\ASRE_Timo2.dll"

'''
This library contains classes and related functions that are specific to the BegrensSkade GIBV programs,
empirical curves for short and long term settlements, 
geotechnical calculators for settlements in time and space,
and vulnerability index functions.
'''

# For debugging or quality control of stepwise, detailed geotechnical calulations,
# the user has the possibility to enable output of a csv log containing such details:
write_csv = False

################ CLASSES FOR CORNER, WALLS AND BUILDINGS ##################

class Corner:
    def __init__(self, cid, x, y):
        self.cid = cid  # corner id to keep track on the line direction
        self.x = x
        self.y = y
        self.dtb = None
        self.sv_short = None
        self.sv_long = None
        self.sh_short = None
        self.near_dist = None
        self.near_angle = None
        self.porewp_red = None
        self.sv_shortSample = None

class Wall:
    def __init__(self, wid, corner_start, corner_end, slope_angle=None, hor_tensile_strain=None, principal_tensile_strain=None):
        self.wid = wid
        self.corner_start = corner_start
        self.corner_end = corner_end
        self.slope_angle = slope_angle
        self.hor_tensile_strain = hor_tensile_strain
        self.principal_tensile_strain = principal_tensile_strain

class RandomFieldmodel:
    #A model to describe stationary random field model
    def __init__(self, mean, cv, distribution, variogramModel, range, nuggetOverVar):
        self.distribution = distribution # string, e.g., 'LogNormal'
        self.mean = mean # mean value
        self.cv = cv # coefficient of variation
        self.vgm = variogramModel #string, e.g., 'Gaussian'
        self.range = range
        self.nuggetOverVar = nuggetOverVar

class RandomVariablemodel:
    #A model to describe stationary random field model
    def __init__(self, mean, cv, distribution):
        self.distribution = distribution # string, e.g., 'LogNormal'
        self.mean = mean # mean value
        self.cv = cv # coefficient of variation


class Building:
    def __init__(self, bid, corners, area=None, circumf=None, foundation=None, structure=None, status=None, logger=None,
                 buildingFeatureListDeterm = None, Es= None, nis = None, buildingFeatureListProb = None):
        self.bid = bid
        self.corners = corners  # array[corner]
        self.walls = []  # array[wall]
        self.adjacent = None

        self.area = area
        self.circumf = circumf
        self.foundation = foundation
        self.structure = structure
        self.status = status

        self.max_sh_sho = None
        self.max_sv_sho = None
        self.max_sv_tot = None
        self.max_angle = None
        self.max_strain = None
        self.max_pstrai = None

        self.vulnerability = None
        self.risk_totset = None
        self.risk_angle = None
        self.logger = logger

        self.maximumStrain = None
        self.damageState = None
        self.logger.debug('In building buildingFeatureListDeterm is not None: '+str(buildingFeatureListDeterm is not None))
        self.logger.debug('In building buildingFeatureListDeterm is: '+str(buildingFeatureListDeterm))
        if buildingFeatureListDeterm is not None:
            self.Eb = buildingFeatureListDeterm[0] * 1.0e9
            self.phi_int = buildingFeatureListDeterm[1]
            self.dfoot = buildingFeatureListDeterm[2]
            self.EoverG = buildingFeatureListDeterm[3]
            self.qz_foot = buildingFeatureListDeterm[4] *1.0e3
            self.bfoot1 = None
            self.bfoot2 = None
            self.Es = Es * 1.0e6
            self.nis = nis
        if buildingFeatureListProb is not None:
            self.Ebmean = buildingFeatureListProb[0]
            self.Ebcv = buildingFeatureListProb[1]
            self.EoverGmean = buildingFeatureListProb[2]
            self.EoverGcv = buildingFeatureListProb[3]
            self.qz_foot_mean = buildingFeatureListProb[4]
            self.qz_foot_cv = buildingFeatureListProb[5]
            self.phi_int = buildingFeatureListProb[6]
            self.dfoot = buildingFeatureListProb[7]
            self.bfoot1 = None
            self.bfoot2 = None
        self.axis1meshX = []#principal axis 1 mesh
        self.axis1meshY = []#principal axis 1 mesh
        self.axis2meshX = []#principal axis 2 mesh
        self.axis2meshY = []#principal axis 1 mesh
        # Greenfield disp of axis 1, V: vertical; L: Longitudinal along axis; T: Translational from axis
        self.dispV1 = np.array([]); self.dispL1 = np.array([]); self.dispT1 = np.array([])
        # Greenfield disp of axis 2, V: vertical; L: Longitudinal along axis; T: Translational from axis
        self.dispV2 = np.array([]); self.dispL2 = np.array([]); self.dispT2 = np.array([])


    def filter_duplicates(self):

        security = 0
        i = 0
        while i < len(self.corners):
            x_prev = self.corners[i - 1].x
            y_prev = self.corners[i - 1].y
            x = self.corners[i].x
            y = self.corners[i].y
            dist = np.sqrt((y - y_prev) ** 2 + (x - x_prev) ** 2)
            if dist < 0.1:
                # Removing duplicate
                del self.corners[i - 1]
            else:
                i += 1

            if security > 10000:
                self.logger.error("While overrun in filter duplicates")
                raise Exception("While overrun in filter duplicates")
            security += 1

    def filter_straights(self, WALL_CORNER_ANGLE_THRESHOLD):
        # remove points that are along a straigth wall
        security = 0
        # i starts at 0, that means, first prev point has index -1, in this way first and last point are compared.
        i = 0
        while i < (len(self.corners) - 1):
            oid = self.corners[i].cid
            x_prev = self.corners[i - 1].x
            y_prev = self.corners[i - 1].y
            x = self.corners[i].x
            y = self.corners[i].y
            x_next = self.corners[i + 1].x
            y_next = self.corners[i + 1].y
            try:
                prev_dir, prev_quadr = Utils.get_angle(x_prev, y_prev, x, y)
                next_dir, next_quadr = Utils.get_angle(x, y, x_next, y_next)
            except Exception("Angle calculation failed, skipping building " + str(self.bid)):
                break

            prev_dir_deg = 180 * prev_dir / np.pi
            next_dir_deg = 180 * next_dir / np.pi

            if abs(prev_dir_deg - next_dir_deg) < WALL_CORNER_ANGLE_THRESHOLD:  # degree
                # Removing straight pnt
                del self.corners[i]
            else:
                i += 1

            if security > 10000:
                self.logger.error("While overrun in filter straights")
                raise Exception("While overrun in filter straights")
            security += 1

    def create_walls(self):
        cid = 0

        for i in range(0, len(self.corners)):

            x_prev = self.corners[i - 1].x
            y_prev = self.corners[i - 1].y
            x = self.corners[i].x
            y = self.corners[i].y

            xy_dist = np.sqrt((x - x_prev) ** 2 + (y - y_prev) ** 2)

            sv_short_prev = self.corners[i - 1].sv_short
            sv_short = self.corners[i].sv_short

            sh_short_prev = self.corners[i - 1].sh_short
            sh_short = self.corners[i].sh_short

            sv_long_prev = self.corners[i - 1].sv_long
            sv_long = self.corners[i].sv_long

            sv_short_term_diff = abs(sv_short_prev - sv_short)
            sh_short_term_diff = abs(sh_short_prev - sh_short)
            sv_total_diff = abs(
                (sv_short_prev + sv_long_prev) - (sv_short + sv_long))

            hor_tensile_strain = sh_short_term_diff / xy_dist

            if xy_dist < 1:
                slope_angle = 0
            else:
                slope_angle = sv_total_diff / xy_dist

            try:
                thetamax = 0.5 * np.arctan(slope_angle / hor_tensile_strain)
            except:
                thetamax = None

            try:
                principal_tensile_strain = hor_tensile_strain * (np.cos(thetamax)) ** 2 + slope_angle * np.sin(thetamax) * np.cos(thetamax)
            except:
                principal_tensile_strain = None

            self.walls.append(Wall(cid,
                                   self.corners[i-1], self.corners[i], slope_angle, hor_tensile_strain, principal_tensile_strain))

            cid += 1
    def equivalentBeamMesh(self, elem_size, logger):
        #TODO 1: find the two principal axis of the building footprint
        #TODO 2: mesh each 
        cornerX = []
        cornerY = []
        principal1 = []
        self.bfoot1 = 1
        self.bfoot2 = 1
        principal2 = []
        # For axis 1
        b_corner1 = self.corners[0] #TODO: change to end1 of axis 1
        b_corner2 = self.corners[1] #TODO: change to end2 of axis 1
        logger.debug("Building {bid} axis 1 corner1 X is {meshx}, Y is {meshy}".
                     format(bid=self.bid, meshx = b_corner1.x, meshy = b_corner1.y))
        logger.debug("Building {bid} axis 1 corner2 X is {meshx}, Y is {meshy}".
                     format(bid=self.bid, meshx = b_corner2.x, meshy = b_corner2.y))
        CONSTR_RESAMPLE_LEN = elem_size
        self.axis1meshX = []
        self.axis1meshY = []
        wall_len = np.sqrt((b_corner1.x-b_corner2.x)**2 + (b_corner1.y-b_corner2.y)**2)
        if wall_len > CONSTR_RESAMPLE_LEN:
            n_segments = int(np.ceil(wall_len/CONSTR_RESAMPLE_LEN))
            seg_len = wall_len / n_segments
            for seg in range(0, n_segments):
                self.axis1meshX.append(b_corner1.x + seg * seg_len * \
                    (b_corner2.x - b_corner1.x) / wall_len)
                self.axis1meshY.append(b_corner1.y + seg * seg_len * \
                    (b_corner2.y - b_corner1.y) / wall_len)
        else:
            self.axis1meshX.append(b_corner1.x)
            self.axis1meshX.append(b_corner2.x)
            self.axis1meshY.append(b_corner1.y)
            self.axis1meshY.append(b_corner2.y)
            logger.error("Building {bid} wall length less than resample lent".format(bid=self.bid))
        logger.debug("Building {bid} axis 1 mesh X is {meshx}".format(bid=self.bid, meshx = str(self.axis1meshX)))
        logger.debug("Building {bid} axis 1 mesh Y is {meshy}".format(bid=self.bid, meshy = str(self.axis1meshY)))
        # For axis 2
        b_corner1 = self.corners[0] #TODO: change to end1 of axis 2
        b_corner2 = self.corners[3] #TODO: change to end2 of axis 2
        logger.debug("Building {bid} axis 2 corner1 X is {meshx}, Y is {meshy}".
                     format(bid=self.bid, meshx = b_corner1.x, meshy = b_corner1.y))
        logger.debug("Building {bid} axis 2 corner2 X is {meshx}, Y is {meshy}".
                     format(bid=self.bid, meshx = b_corner2.x, meshy = b_corner2.y))
        CONSTR_RESAMPLE_LEN = elem_size
        wall_len = np.sqrt((b_corner1.x-b_corner2.x)**2 + (b_corner1.y-b_corner2.y)**2)
        if wall_len > CONSTR_RESAMPLE_LEN:
            n_segments = int(np.ceil(wall_len/CONSTR_RESAMPLE_LEN))
            seg_len = wall_len / n_segments
            for seg in range(0, n_segments):
                self.axis2meshX.append(b_corner1.x + seg * seg_len * \
                    (b_corner2.x - b_corner1.x) / wall_len)
                self.axis2meshY.append(b_corner1.y + seg * seg_len * \
                    (b_corner2.y - b_corner1.y) / wall_len)
        else:
            self.axis2meshX.append(b_corner1.x)
            self.axis2meshX.append(b_corner2.x)
            self.axis2meshY.append(b_corner1.y)
            self.axis2meshY.append(b_corner2.y)
            logger.error("Building {bid} wall length less than resample lent".format(bid=self.bid))
        logger.debug("Building {bid} axis 2 mesh X is {meshx}".format(bid=self.bid, meshx = str(self.axis2meshX)))
        logger.debug("Building {bid} axis 2 mesh Y is {meshy}".format(bid=self.bid, meshy = str(self.axis2meshY)))
    def greenFieldDisp_zhao2022_determine(self, dvmaxOverHe, eta, excavation_depth, dlmaxOverdvmax, construction_area_corners, logger):
    #TODO 1: calculate dispX, dispY, dispZ for axis 1 and axis 2
    #TODO 2: calculate dispV, dispL, dispT for axis 1 and axis 2
        dispV1 = np.zeros(len(self.axis1meshX))
        dispL1 = np.zeros(len(self.axis1meshX))
        dispT1 = np.zeros(len(self.axis1meshX))
        dispV2 = np.zeros(len(self.axis2meshX))
        dispL2 = np.zeros(len(self.axis2meshX))
        dispT2 = np.zeros(len(self.axis2meshX))
        # For axis 1
        near_dist_corner1 = np.sqrt(near_analysis_sqr(self.axis1meshX[ 0], self.axis1meshY[ 0], construction_area_corners))
        near_dist_corner2 = np.sqrt(near_analysis_sqr(self.axis1meshX[-1], self.axis1meshY[-1], construction_area_corners))
        #axis1_angle is the angel from north to building axis using the axis end closer to excavation as origin, unit is degree
        if near_dist_corner2 < near_dist_corner1:
            # axis1_angle = Utils.getAngleFromDir(self.axis1meshX[0], self.axis1meshY[0],
            #                                         self.axis1meshX[-1], self.axis1meshY[-1])
            vecAxis1 = np.array([self.axis1meshX[0] - self.axis1meshX[-1], self.axis1meshY[0] - self.axis1meshY[-1]])
        else:
            # axis1_angle = Utils.getAngleFromDir(self.axis1meshX[-1], self.axis1meshY[-1],
            #                                         self.axis1meshX[0], self.axis1meshY[0])
            vecAxis1 = np.array([self.axis1meshX[-1] - self.axis1meshX[0], self.axis1meshY[-1] - self.axis1meshY[0]])
        #rotate vecAxis1 clockwise by 90 degree: https://limnu.com/sketch-easy-90-degree-rotate-vectors/#:~:text=Normally%20rotating%20vectors%20involves%20matrix,swap%20X%20and%20Y%20values.
        vecAxis1_perpend = np.array([vecAxis1[1], -vecAxis1[0]])
        for i in range(len(self.axis1meshX)):
            # near_dist, near_angle = near_analysis(
            #     self.axis1meshX[i], self.axis1meshY[i], construction_area_corners)
            near_point_ind, near_dist_sqr = near_analysis_sqrNcorner(
                self.axis1meshX[i], self.axis1meshY[i], construction_area_corners)
            near_dist = np.sqrt(near_dist_sqr)
            #near_angle is the angle from east to near_point using building corner as origin
            dispV1[i]=get_sv_short_Zhao2022(near_dist, dvmaxOverHe, eta, excavation_depth)
            dispH = get_sh_short_Zhao2022(near_dist, dlmaxOverdvmax*dvmaxOverHe, eta, excavation_depth)
            vecdH = np.array([construction_area_corners[near_point_ind].x-self.axis1meshX[i], 
                             construction_area_corners[near_point_ind].y-self.axis1meshY[i]])
            dispL1[i] = dispH/near_dist * np.dot(vecdH, vecAxis1)/np.linalg.norm(vecAxis1)
            dispT1[i] = dispH/near_dist * np.dot(vecdH, vecAxis1_perpend)/np.linalg.norm(vecAxis1_perpend)
        
        # For axis 2
        near_dist_corner1 = np.sqrt(near_analysis_sqr(self.axis2meshX[ 0], self.axis2meshY[ 0], construction_area_corners))
        near_dist_corner2 = np.sqrt(near_analysis_sqr(self.axis2meshX[-1], self.axis2meshY[-1], construction_area_corners))
        #axis1_angle is the angel from north to building axis using the axis end closer to excavation as origin, unit is degree
        if near_dist_corner2 < near_dist_corner1:
            # axis1_angle = Utils.getAngleFromDir(self.axis1meshX[0], self.axis1meshY[0],
            #                                         self.axis1meshX[-1], self.axis1meshY[-1])
            vecAxis2 = np.array([self.axis2meshX[0] - self.axis2meshX[-1], self.axis2meshY[0] - self.axis2meshY[-1]])
        else:
            # axis1_angle = Utils.getAngleFromDir(self.axis1meshX[-1], self.axis1meshY[-1],
            #                                         self.axis1meshX[0], self.axis1meshY[0])
            vecAxis2 = np.array([self.axis2meshX[-1] - self.axis2meshX[0], self.axis2meshY[-1] - self.axis2meshY[0]])
        #rotate vecAxis2 clockwise by 90 degree: https://limnu.com/sketch-easy-90-degree-rotate-vectors/#:~:text=Normally%20rotating%20vectors%20involves%20matrix,swap%20X%20and%20Y%20values.
        vecAxis2_perpend = np.array([vecAxis2[1], -vecAxis2[0]])
        for i in range(len(self.axis2meshX)):
            # near_dist, near_angle = near_analysis(
            #     self.axis1meshX[i], self.axis1meshY[i], construction_area_corners)
            near_point_ind, near_dist_sqr = near_analysis_sqrNcorner(
                self.axis2meshX[i], self.axis2meshY[i], construction_area_corners)
            near_dist = np.sqrt(near_dist_sqr)
            #near_angle is the angle from east to near_point using building corner as origin
            dispV2[i]=get_sv_short_Zhao2022(near_dist, dvmaxOverHe, eta, excavation_depth)
            dispH = get_sh_short_Zhao2022(near_dist, dlmaxOverdvmax*dvmaxOverHe, eta, excavation_depth)
            vecdH = np.array([construction_area_corners[near_point_ind].x-self.axis2meshX[i], 
                             construction_area_corners[near_point_ind].y-self.axis2meshY[i]])
            dispL2[i] = dispH/near_dist * np.dot(vecdH,vecAxis2)/np.linalg.norm(vecAxis2)
            dispT2[i] = dispH/near_dist * np.dot(vecdH,vecAxis2_perpend)/np.linalg.norm(vecAxis2_perpend)
        logger.debug('dispL1 is: '+str(dispL1))
        return [dispL1, dispT1, dispV1], [dispL2, dispT2, dispV2]
    
    def calculate_damageState_ASRE_timo(self, axis1_gf, axis2_gf, logger):
        #TODO 1: load ASRE with CDLL
        logger.debug('ASRE_path is?: '+str(ASRE_path))
        ASRE_clib = CDLL(r"C:\Users\Jinyan\Desktop\NGI_GIBV\ASRETimoshenko\ASRE_Timo2\x64\Release\ASRE_Timo2.dll")
        logger.debug('ASRE_clib is: '+str(ASRE_clib))
        ASRE_clib.run.argtypes = [c_int, 
                          np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.float64), 
                              np.ctypeslib.ndpointer(dtype=np.float64), 
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double
                            #   POINTER(c_double)
                            #   np.ctypeslib.ndpointer(dtype=np.float64)
                              ]
        ASRE_clib.run.restype = POINTER(c_double * 3)
        logger.debug('ASRE_clib is 2: '+str(ASRE_clib))
        #TODO 2: calcualte principal strain
        #For axis 1:
        nnode = len(self.axis1meshX)
        dispV = axis1_gf[2]; dispL = axis1_gf[0]; dispT = axis1_gf[1];
        ni_foot = self.EoverG/2-1
        mu_int = np.tan(self.phi_int/180.0*np.pi)
        logger.debug('ASRE_clib is 3: '+str(ASRE_clib))
        logger.debug('type(dispV) is 3: '+str(type(dispV)))
        logger.debug('type(nnode) is 3: '+str(type(nnode)))
        meshX = np.array(self.axis1meshX)
        meshY = np.array(self.axis1meshY)
        meshZ = np.zeros(nnode)
        logger.debug("JINYAN before ASRE_clib, {0}; {0};".format(
            type(self.Eb), type(dispV)
        ))

        ###########Test values#####################
        # nnode = 31
        # meshX = np.linspace(-15, 15, nnode)
        # meshY = np.zeros(nnode)
        # meshZ = np.zeros(nnode)
        # dispV = np.ones(nnode)/1000
        # dispL = np.ones(nnode)/1000
        # dispT = np.zeros(nnode)
        # Eb = 70000000000 #70GPa
        # dfoot = 1.5 #1.5m
        # EoverG = 2.6
        # bfoot1 = 10 #10m
        # EsNominal = 25000000 #25MPa
        # nis = 0.25
        # ni_foot = 0.3 #nu = E/2G - 1 https://en.wikipedia.org/wiki/Elastic_modulus
        # mu_int = np.tan(30/180*np.pi)
        # qz_foot = 40.5*1000*10*30/30 #40.5KPa on a 10*30 foot
        ###########################################
        result = 0
        try:
            result = ASRE_clib.run(nnode, meshX, meshY, meshZ, dispV, dispL, dispT, self.Eb, self.EoverG, 
                       self.Es, self.nis, self.dfoot, self.bfoot1, ni_foot, mu_int, self.qz_foot)
            # result = ASRE_clib.run(nnode, meshX, meshY, meshZ,
            #                     dispV, dispL, dispT, self.Eb, self.EoverG, 
            #             self.Es, self.nis, self.dfoot, self.bfoot1, ni_foot, mu_int, self.qz_foot)
            resultList = [i for i in result.contents]
            resultArray1 = np.array(resultList)
        except:
            libHandle = ASRE_clib._handle
            del ASRE_clib
            _ctypes.FreeLibrary(libHandle)
            logger.error('ASRE_clib.run() failed at axis 1 of building bid: '+str(self.bid))
        if result:
            # libHandle = ASRE_clib._handle
            # del ASRE_clib
            # _ctypes.FreeLibrary(libHandle)
            logger.info('ASRE_clib.run() completed at axis 1 of building bid: '+str(self.bid)+
                       'resultArray is: ' + str(resultArray1))

        #For axis 2:
        nnode = len(self.axis2meshX)
        dispV = axis2_gf[2]; dispL = axis2_gf[0]; dispT = axis2_gf[1];
        ni_foot = self.EoverG/2-1
        mu_int = np.tan(self.phi_int/180.0*np.pi)
        meshX = np.array(self.axis2meshX)
        meshY = np.array(self.axis2meshY)
        meshZ = np.zeros(nnode)
        result = 0
        try:
            result = ASRE_clib.run(nnode, meshX, meshY, meshZ, dispV, dispL, dispT, self.Eb, self.EoverG, 
                       self.Es, self.nis, self.dfoot, self.bfoot2, ni_foot, mu_int, self.qz_foot)
            # result = ASRE_clib.run(nnode, meshX, meshY, meshZ,
            #                     dispV, dispL, dispT, self.Eb, self.EoverG, 
            #             self.Es, self.nis, self.dfoot, self.bfoot1, ni_foot, mu_int, self.qz_foot)
            resultList = [i for i in result.contents]
            resultArray2 = np.array(resultList)
        except:
            libHandle = ASRE_clib._handle
            del ASRE_clib
            _ctypes.FreeLibrary(libHandle)
            logger.error('ASRE_clib.run() failed at axis 2 of building bid: '+str(self.bid))
        if result:
            libHandle = ASRE_clib._handle
            del ASRE_clib
            _ctypes.FreeLibrary(libHandle)
            logger.info('ASRE_clib.run() completed at axis 2 of building bid: '+str(self.bid)+
                       'resultArray is: ' + str(resultArray2))
        # #TODO 3: find damage state with limiting strain method
        maximumStrain = max(resultArray1.max(), resultArray2.max())
        if maximumStrain < 0.05/100.0:
            DS = 0
        elif maximumStrain < 0.075/100.0:
            DS = 1
        elif maximumStrain < 0.15/100.0:
            DS = 2
        elif maximumStrain < 0.3/100.0:
            DS = 3
        else:
            DS = 4
        return maximumStrain, DS
    
    def calculate_damageState_ASRE_timoMC(self, logger):
        #For each sample
            #TODO 1: load ASRE with CDLL
            #TODO 2: calcualte principal strain
        #TODO 3: calculate the statistics
        #TODO 4: 
        # Test for multiple processing
        import joblib, multiprocessing
        multiprocessing.set_executable(os.path.join(sys.exec_prefix, 'python.exe'))
        num_cores = multiprocessing.cpu_count()
        logger.debug("multiprocessing test" + str(os.path.join(sys.exec_prefix, 'python.exe')) + str(multiprocessing.cpu_count()))
        input = range(10)
        jobs = []
        for i in input:
            jobs.append((logger, i))
        with multiprocessing.Pool(processes=2) as pool: # Create the pool object 
            res = pool.starmap(runParallel_SSI, jobs)    
        # res = joblib.Parallel(n_jobs=num_cores)(joblib.delayed(BegrensSkadeLib.runParallel_SSI)(logger, i) for i in input)
        logger.debug("multiprocessing results" + str(res))
        return 1, 2

def runParallel_SSI(logger, i):
    log_path = r'C:\Users\Jinyan\Documents\ArcGIS\Projects\CampusUlleval2\REMEDY_GIS_RiskTool_JZ\log\subProcess'
    curr_proc = multiprocessing.current_process()
    with open(log_path + os.sep + "GIBV_Excavation_subprocess_" +str(curr_proc.name) + ".txt", "a") as debug_log:
        debug_log.write('current process:' +str(curr_proc.name) + str(curr_proc._identity) + 'i: ' + str(i) + "\n")
    return i

####################### SHORT TERM SETTLEMENT CURVES #######################

def get_sv_short_a(near_dist, excavation_depth):
    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 3

    if x < 0.3:
        return [(0.15 + ((0.5-0.15)/0.3)*x)*excavation_depth*0.01, W]

    elif x >= 0.3 and x < 1:
        return [(0.5 + ((0.1-0.5)/(1-0.3))*(x-0.3))*excavation_depth*0.01 , W]

    elif x >= 1 and x < 1.5:
        return [(0.1 + ((-0.1)/(1.5-1))*(x-1)) *excavation_depth*0.01, W]

    else:
        return [0.0, W]

def get_sv_short_b(near_dist, excavation_depth):

    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 3

    if x < 0.5:
        return [(0.3 + ((1-0.3)/0.5)*x)*excavation_depth*0.01, W]

    elif x >= 0.5 and x < 2:
        return [(1 + ((0.2-1)/(2-0.5))*(x-0.5))*excavation_depth*0.01 , W]

    elif x >= 2 and x < 3:
        return [(0.2 + ((-0.2)/(3-2))*(x-2))*excavation_depth*0.01,W]

    else:
        return [0.0, W]

def get_sv_short_c(near_dist, excavation_depth):

    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 5

    if x < 0.7:
        return [(1 + ((2-1)/0.7)*x)*excavation_depth*0.01, W]

    elif x >= 0.7 and x < 2.5:
        return [(2 + ((0.5-2)/(2.5-0.7))*(x-0.7))*excavation_depth*0.01 , W]

    elif x >= 2.5 and x < 4:
        return [(0.5 + ((-0.5)/(4-2.5))*(x-2.5)) *excavation_depth*0.01, W]

    else:
        return [0.0, W]

def get_sv_short_d(near_dist, excavation_depth):

    x = near_dist / excavation_depth #normalized distance from byggegrop
    W = 5

    if x < 1:
        return [(1.5 + ((3 - 1.5) / 1) * x)*excavation_depth*0.01, W]

    elif x >= 1 and x < 3:
        return [(3 + ((0.75 - 3) / (3 - 1)) * (x - 1))*excavation_depth*0.01, W]

    elif x >= 3 and x < 5:
        return [(0.75 + ((-0.75) / (5 - 3)) * (x - 3))*excavation_depth*0.01, W]

    else:
        return [0.0, W]


################## PECK SETTLEMENT CURVE FOR SOIL TUNNELS ##################

def get_sv_short_Peck(near_dist, tunnel_depth, tunnel_diameter, volume_loss, trough_width):
    D = tunnel_diameter
    i = trough_width * tunnel_depth
    s0 = volume_loss/100*((np.pi*D**2)/4)/(np.sqrt(2*np.pi)*i) 
    s = s0 * np.exp(-(near_dist**2)/(2*i**2)) 
    return s 


################## JANBU LONG TERM SETTLEMENT CALCULATION #################

def get_sv_long_janbu(dtb, dry_crust_thk, dep_groundwater, density_sat, OCR, porewp_red, p_ref, janbu_const, janbu_m, consolidation_time):
    density_water = 10  # kN/m3
    permeability = 1e-9  # (m/s)
    adj = False
    dep_to_clay = max(dry_crust_thk, dep_groundwater)
    clay_thk = dtb - dep_to_clay
    sv_acc = 0  # accumulated horizontal settlement in clay layer

    janbu_M_avgs = []  #depth average of janbu_M to be used in terzagi time equation
    METHOD = "LYR"
    #METHOD = "AVG"

    for clay_dep_i in range(int(round(dep_to_clay, 0)), int(round(dtb, 0))):

        clay_dep_i += 0.5
        pz0 = density_sat * clay_dep_i
        uz0 = density_water * (clay_dep_i - dep_groundwater)
        pz0_eff = pz0 - uz0
        pz_pre = OCR * pz0_eff

        if pz_pre < p_ref:
            #This can happen right below the surface where pz0_eff is low, implies negative stiffness and does not make sense
            raise Exception("Negative stiffness in Janbu model. Try to decrease the reference pressure or increase the OCR")

        #porewp_red = porewp_red * longterm_porewr

        uz0_tot = density_water * (dtb - dep_groundwater)
        if porewp_red > uz0_tot:
            # arcpy.AddMessage("porewp_red_atdist: " + str(porewp_red_atdist))
            # arcpy.AddMessage("density_water * (dtb- dep_groundwater): " + str(density_water * (dtb - dep_groundwater)))
            adj = True
            porewp_red = density_water * (dtb - dep_groundwater)

        clay_dep_lyr = (clay_dep_i - dep_to_clay)

        du = (porewp_red/clay_thk)*clay_dep_lyr
        janbu_M = janbu_const * janbu_m * pz_pre
        janbu_M_avgs.append(janbu_M)

        if pz0_eff + du < pz_pre:
            # implicitly multiplied by 1 m clay thickness (iteration step)
            sv = du / janbu_M
            # regime = "flat"

        else:
            # implicitly multiplied by 1 m clay thickness (iteration step)
            sv = (pz_pre - pz0_eff)/janbu_M + (1/janbu_m) * \
                np.log((pz0_eff + du - p_ref)/(pz_pre - p_ref))
            # regime = "linear"

        # Calculation of time-dependent consolidation. Default is maxiumn consolidation (1000 years - gives 99.9 % consolidation for 50 m clay)
        if METHOD == "LYR":
            t = 60 * 60 * 24 * 365 * consolidation_time  # sec
            c = janbu_M * permeability / density_water
            T = c * t / (clay_thk) ** 2
            sv_inf = sv
            sv = sv * U_draintop_b(T, 10)

        sv_acc += sv




    #Calculation of time-dependent consolidation. Default is maxiumn consolidation (1000 years - gives 99.9 % consolidation for 50 m clay)
    if METHOD == "AVG" and sv_acc > 0:
        janbu_M_avg = sum(janbu_M_avgs)/len(janbu_M_avgs)
        t = 60 * 60 * 24 * 365 * consolidation_time  # sec
        c = janbu_M_avg * permeability / density_water
        T = c * t / (clay_thk) ** 2
        sv_acc = sv_acc * U_draintop_b(T, 10)

    return sv_acc, adj


################## LONG TERM POREWATER REDUCTION CURVES ###################

def get_longterm_porewr_min(near_dist):
    x = near_dist
    # water_density = 10  # kN/m3
    if x < 220:

        # old polynomial regresseion, not optimal
        # return (0.000002*x**2 - 0.0044*x + 0.7059) #* water_density * byggegrop_depth #AOL may 6th 2019: Uncertain about this formula.

        J = 25
        K = 0.62
        L = 0.6
        M = 0.96

        return M * K * np.exp((-(x ** L) / J)) + (1 - M) * K * (-x / J)

    else:
        return 0

def get_longterm_porewr_max(near_dist):
    x = near_dist
    # water_density = 10  # kN/m3
    if x < 360:

        # old polynomial regresseion, not optimal
        # return (0.000002*x**2 - 0.0044*x + 0.7059) #* water_density * byggegrop_depth #AOL may 6th 2019: Uncertain about this formula.

        J = 26
        K = 1
        L = 0.625
        M = 0.985

        return M * K * np.exp((-(x ** L) / J)) + (1 - M) * K * (-x / J)

    else:
        return 0

def get_longterm_porewr_mean(near_dist):
    x = near_dist
    # water_density = 10  # kN/m3
    if x < 340:

        # old polynomial regresseion, not optimal
        # return (0.000002*x**2 - 0.0044*x + 0.7059) #* water_density * byggegrop_depth #AOL may 6th 2019: Uncertain about this formula.

        J = 30
        K = 0.8
        L = 0.68
        M = 0.985

        return M * K * np.exp((-(x ** L) / J)) + (1 - M) * K * (-x / J)

    else:
        return 0


#################### TERZAGI CONSOLIDATION TIME CURVES ####################

def U_drainboth(T, m_max):

    sum = 0

    for m in range(0, m_max):
        M = 0.5 * np.pi * (2*m + 1)

        sum += (2/(M**2))*np.exp(-T*M**2)

    return 1 - sum

def U_draintop_a(T, m_max): #Trekant med spiss oppad

    sum = 0

    for m in range(0, m_max):

        M = 0.5 * np.pi * (2*m - 1)

        sum += 2*((-1)**(m+1))*np.exp(-T*M**2)/M**3

    return 1 - sum

def U_draintop_b(T, m_max): #Trekant med spiss nedad

    return 2*U_drainboth(T, m_max) - U_draintop_a(T, m_max)


######################## VULNERABILITY INDEX FUNCTIONS ####################

def get_buil_len_cvi(buil_len):
    if buil_len < 10:
        return 0
    elif buil_len < 15:
        return 5
    elif buil_len < 30:
        return 20
    else:
        return 50

def get_buil_shape_cvi(squareness):
    if squareness < 0.35:
        return 50
    elif squareness < 0.5:
        return 20
    elif squareness < 0.75:
        return 5
    else:
        return 0

def get_buil_impact_totset_cvi(tot_set):
    if tot_set < 0.010:
        return 1
    elif tot_set < 0.050:
        return 2
    elif tot_set < 0.075:
        return 3
    else:
        return 4

def get_buil_impact_angle_cvi(angle):
    if angle < (1/500):
        return 1
    elif angle < (1/200):
        return 2
    elif angle < (1/50):
        return 3
    else:
        return 4

def get_buil_vuln_cvi(vuln):
    if vuln < 0.25:
        return 1
    elif vuln < 0.5:
        return 2
    elif vuln < 0.75:
        return 3
    else:
        return 4

def get_risk_cvi(vuln_cvi, impact_cvi):
    if vuln_cvi * impact_cvi < 3:
        return 1
    elif vuln_cvi * impact_cvi < 5 and vuln_cvi < 4 and impact_cvi < 4:
        return 2
    elif vuln_cvi * impact_cvi < 8:
        return 3
    elif vuln_cvi * impact_cvi < 12:
        return 4
    else:
        return 5


################ GIBV SPECIFIC SHAPEFILE I/O FUNCTIONALITY #################

def get_construction_corners(shapeFile, CONSTR_RESAMPLE_LEN, logger):
    construction_area_corner_pnts = []
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile, 0)
    layer = dataSource.GetLayer()
    cid = 0

    for feature in layer:
        # gjson = feature.ExportToJson()
        # pj = json.loads(gjson)
        # geom = feature.geometry()

        cid = 1

        for part in feature.geometry():
            # prev_X = part[0].X
            # prev_Y = part[0].Y
            first = True
            for p in part.GetPoints():
                X = p[0]
                Y = p[1]
                if first:
                    prev_X = X
                    prev_Y = Y
                    first = False
                    continue

                # print("X: {}  Y:{}".format(x,y))

                wall_len = np.sqrt((X-prev_X)**2 + (Y-prev_Y)**2)
                if wall_len > CONSTR_RESAMPLE_LEN:
                    n_segments = int(np.ceil(wall_len/CONSTR_RESAMPLE_LEN))
                    seg_len = wall_len / n_segments
                    for seg in range(0, n_segments):
                        X_sampl = prev_X + seg * seg_len * \
                            (X - prev_X) / wall_len
                        Y_sampl = prev_Y + seg * seg_len * \
                            (Y - prev_Y) / wall_len

                        construction_area_corner_pnts.append(
                            Corner(cid, X_sampl, Y_sampl))
                else:
                    construction_area_corner_pnts.append(Corner(cid, X, Y))
                cid += 1
                prev_X = X
                prev_Y = Y
    return construction_area_corner_pnts

def get_buildings(features, fieldNameFoundation=None, fieldNameStructure=None, fieldNameStatus=None,
                   buildingFeatureFieldNameListDeterm = None, buildingFeatureFieldNameListProb = None, Es = None, nis = None,logger=None) -> []:
    input_buildings = []
    bid = 1
    numBuildings = 0
    numSkippedBuildings = 0

    if (not logger is None):
        logger.debug("Opening shapefile: {}".format(features))

    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = shp_driver.Open(features, 1)
    layer = dataSource.GetLayer()
    foundation = None
    structure = None
    status = None

    for buildingFeature in layer:
        numBuildings +=1
        cid = 0

        g = buildingFeature.geometry().GetGeometryRef(0)

        if (g is None):
            # logger.debug("Geometry is None, jumping to next")
            continue
        if (g.GetPoints() is None):
            # logger.debug("g.GetPoints() is None, jumping to next")
            continue

        corners = []

        for pnt in g.GetPoints():
            corners.append(Corner(cid, pnt[0], pnt[1]))
            cid += 1

        if (not fieldNameFoundation is None):
            foundation = buildingFeature.GetField(fieldNameFoundation)
        if (not fieldNameStructure is None):
            structure = buildingFeature.GetField(fieldNameStructure)
        if (not fieldNameStatus is None):
            status = buildingFeature.GetField(fieldNameStatus)
        
        buildingFeatureListDeterm = None
        if buildingFeatureFieldNameListDeterm is not None:
            buildingFeatureListDeterm = []
            for i in buildingFeatureFieldNameListDeterm:
                buildingFeatureListDeterm.append(
                    buildingFeature.GetField(i)
                )
        buildingFeatureListProb = None
        if buildingFeatureFieldNameListProb is not None:
            buildingFeatureListProb = []
            for i in buildingFeatureFieldNameListProb:
                buildingFeatureListProb.append(
                    buildingFeature.GetField(i)
                )
        input_buildings.append(
            Building(bid, corners, g.GetArea(), g.Length(),
                     foundation=foundation, structure=structure, status=status, 
                     buildingFeatureListDeterm=buildingFeatureListDeterm,
                     Es=Es, nis=nis,
                     buildingFeatureListProb=buildingFeatureListProb,logger=logger))
        #logger.debug(f"Building {bid}, num corners: {len(corners)}, area: {g.GetArea()}, length: {g.Length()} - get_buildings")
        bid += 1
    return input_buildings

def get_buildings_with_dtb(features, dtb_filename, fieldNameFoundation=None, fieldNameStructure=None, fieldNameStatus=None, logger=None) -> []:
    input_buildings = []
    bid = 1
    numBuildings = len(input_buildings)
    numSkippedBuildings = 0

    if (not logger is None):
        logger.debug("Opening shapefile: {}".format(features))

    src_ds = gdal.Open(str(dtb_filename))
    gt = src_ds.GetGeoTransform()
    rb = src_ds.GetRasterBand(1)
    no_data_value = rb.GetNoDataValue()
    logger.debug(f"No data value for raster = {no_data_value}")

    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = shp_driver.Open(features, 1)
    layer = dataSource.GetLayer()
    foundation = None
    structure = None
    status = None

    for buildingFeature in layer:

        cid = 0

        g = buildingFeature.geometry().GetGeometryRef(0)

        if (g is None):
            # logger.debug("Geometry is None, jumping to next")
            continue
        if (g.GetPoints() is None):
            # logger.debug("g.GetPoints() is None, jumping to next")
            continue

        corners = []
        foundNullVal = False


        for pnt in g.GetPoints():
            if pnt:
                px = int((pnt[0] - gt[0]) / gt[1])  # x pixel
                py = int((pnt[1] - gt[3]) / gt[5])  # y pixel
                dtb = rb.ReadAsArray(px, py, 1, 1)                
                if dtb is None or dtb == no_data_value:
                    numSkippedBuildings += 1
                    foundNullVal = True
                    #if one of the corners has no dbt value, skip to next building
                    #logger.debug("Skipping buildingNumber of skipped buildings: {}".format(numSkippedBuildings))
                    break
                corners.append(Corner(cid, pnt[0], pnt[1]))              
                cid += 1

        if foundNullVal:
            continue
        if (not fieldNameFoundation is None):
            foundation = buildingFeature.GetField(fieldNameFoundation)
        if (not fieldNameStructure is None):
            structure = buildingFeature.GetField(fieldNameStructure)
        if (not fieldNameStatus is None):
            status = buildingFeature.GetField(fieldNameStatus)

        input_buildings.append(
            Building(bid, corners, g.GetArea(), g.Length(),
                     foundation=foundation, structure=structure, status=status, logger=logger))
        bid += 1
    logger.debug(f"Number of skipped buildings: {numSkippedBuildings} of {numBuildings} - get_buildings_with_dtb")
    return input_buildings

def get_construction_corners_from_ArcGIS_json(jsonBody, CONSTR_RESAMPLE_LEN, logger):
    construction_area_corner_pnts = []

    # logger.debug("json: {}".format(jsonBody))
    if "features" in jsonBody:
        logger.debug("FINNES")
        jsonBody = jsonBody["features"]
        if "geometry" in jsonBody[0]:
            logger.debug("features FINNES")
            jsonBody = jsonBody[0]["geometry"]
        else:
            logger.debug("features FINNES IKKE")
    else:
        logger.debug("Cant fint fatures in json body, getting rings directly")

    # logger.debug("json etter avskrelling: {}".format(jsonBody))
    # features = jsonBody["features"]
    # geometry = features[0]["geometry"]
    rings = jsonBody["rings"]
    # logger.debug("rings: {}".format(rings))
    for ring in rings:
        # gjson = feature.ExportToJson()
        # pj = json.loads(gjson)
        # geom = feature.geometry()
        logger.debug("Getting corners in excavation")

        cid = 1
        first = True
        for p in ring:
            # logger.debug("point: {}".format(p))
            # for p in part.GetPoints():
            X = p[0]
            Y = p[1]
            if first:
                prev_X = X
                prev_Y = Y
                first = False
                continue

            # print("X: {}  Y:{}".format(x,y))

            wall_len = np.sqrt((X-prev_X)**2 + (Y-prev_Y)**2)
            if wall_len > CONSTR_RESAMPLE_LEN:
                n_segments = int(np.ceil(wall_len/CONSTR_RESAMPLE_LEN))
                seg_len = wall_len / n_segments
                for seg in range(0, n_segments):
                    X_sampl = prev_X + seg * seg_len * \
                        (X - prev_X) / wall_len
                    Y_sampl = prev_Y + seg * seg_len * \
                        (Y - prev_Y) / wall_len

                    construction_area_corner_pnts.append(
                        Corner(cid, X_sampl, Y_sampl))
            else:
                construction_area_corner_pnts.append(Corner(cid, X, Y))
            cid += 1
            prev_X = X
            prev_Y = Y
    # logger.debug("construction_area_corner_pnts: {}".format(
    #    construction_area_corner_pnts))
    return construction_area_corner_pnts

def createBuildingCornersDict(shapeFile, fieldNames, bLongterm, short_term_curve, byggegrop_depth, hor_vert_ratio, dry_crust_thk, dep_groundwater,
                              density_sat, OCR, porewp_red, janbu_ref_stress, janbu_m, pw_reduction_curve, logger) -> {}:

    bid_prev = "start"
    pnt_array = []
    result_building_corner_dict = {}

    first_time = True

    # Loop through the newly created point feature that contains the building corners.
    # Sh and sv and sv_long are calculated for every points.
    # Store corner information in dictionary building_corners_dict.
    # The key is the building_id (normally object ID), the value is a table containing the corner points for that particular building.
    # The purpose of this is to have a faster accessible data structure for the corner points, that can be edited/fitered for duplicates and straigt-wall points.
    count_all = 0
    count_adj = 0

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapeFile, 0)
    layer = dataSource.GetLayer()

    for pnt in layer:
        oid = pnt.GetField(0)
        bid = pnt.GetField(1)
        x = pnt.GetField(2)
        y = pnt.GetField(3)
        near_dist = pnt.GetField(4)

        if bLongterm:
            dtb = pnt.GetField(6)  # 5??

        # arcpy.AddMessage(str(oid) + ", " + str(bid) + ", " + str(x) + ", " + str(y) + ", " + str(near_dist) + ", " + str(angle))


        if short_term_curve in ["Wall not to bedrock - regression", "0,5 % av byggegropdybde", "Spunt installert til berg med høy sikkerhet", "Norm_setning_0.5"]:
            sv_short, W = get_sv_short_a(near_dist, byggegrop_depth)
        elif short_term_curve in ["Wall not to bedrock - discrete", "1 % av byggegropdybde", "Spunt installert til berg med lav sikkerhet", "Norm_setning_1"]:
            sv_short, W = get_sv_short_b(near_dist, byggegrop_depth)
        elif short_term_curve in ["Tie-back anchors - regression", "2 % av byggegropdybde", "Svevespunt høy sikkerhet", "Norm_setning_2"]:
            sv_short, W = get_sv_short_c(near_dist, byggegrop_depth)
        elif short_term_curve in ["Tie-back anchors - regression", "3 % av byggegropdybde", "Svevespunt lav sikkerhet", "Norm_setning_3"]:
            sv_short, W = get_sv_short_d(near_dist, byggegrop_depth)
        else:
            raise Exception("Not a valid regression curve: " + str(short_term_curve) )

        norm_dist = near_dist / byggegrop_depth  # normalized distance from byggegrop

        sh_short = -hor_vert_ratio * (1 + 2*norm_dist/float(W))*sv_short

        if bLongterm:
            sv_long, red_adj = get_sv_long_janbu(
                 dtb, dry_crust_thk, dep_groundwater, density_sat, OCR, porewp_red, janbu_ref_stress, janbu_m)
            if red_adj:
                count_adj += 1
        else:
            sv_long = 0

        if bid_prev == "start":
            bid_prev = bid

        if bid != bid_prev:
            # New building. Store line data.
            # arcpy.AddMessage("Adding new building: " + str(bid_prev) + ", " + str(len(pnt_array)))
            result_building_corner_dict[bid_prev] = pnt_array
            pnt_array = []
            # arcpy.AddMessage([oid, x, y, near_dist, sv_short, sh_short, sv_long])
            pnt_array.append(
                [oid, x, y, near_dist, sv_short, sh_short, sv_long])

        else:
            pnt_array.append(
                [oid, x, y, near_dist, sv_short, sh_short, sv_long])

        bid_prev = bid
        count_all += 1

    if bLongterm:
        logger.info("Percentage adjustment of porewp red. due to shallow bedrock : " +
                    str(round(100*(count_adj/count_all), 1)))

    # write the last building
    result_building_corner_dict[bid_prev] = pnt_array
    return result_building_corner_dict

def near_analysis_sqr(xref, yref, construction_area_corners):

    near_dist_sqr = 999999

    for corner in construction_area_corners:
        x = corner.x
        y = corner.y

        dist_sqr = (xref-x)**2+(yref-y)**2

        if dist_sqr < near_dist_sqr:
            near_dist_sqr = dist_sqr

    return near_dist_sqr

def near_analysis_sqrNcorner(xref, yref, construction_area_corners):

    near_dist_sqr = 999999
    near_dist_corner_ind = 0

    for ind, corner in enumerate(construction_area_corners):
        x = corner.x
        y = corner.y

        dist_sqr = (xref-x)**2+(yref-y)**2

        if dist_sqr < near_dist_sqr:
            near_dist_sqr = dist_sqr
            near_dist_corner_ind = ind

    return near_dist_corner_ind, near_dist_sqr

def gaussianAutoCorr(h, c, nug, r):
    if h!=0:
        gamma_h = c*(1-np.exp(-h/r))+nug
    else:
        gamma_h = nug
    gamma_inf = c+nug
    rho = 1-gamma_h/gamma_inf
    return rho

#https://scipy-cookbook.readthedocs.io/items/CorrelatedRandomSamples.html
def generateCorrGauSample(sample_size, eig_vec, eig_val):
    c = np.dot(eig_vec, np.diag(np.sqrt(eig_val)))
    x = scipy.stats.norm.rvs(size = (eig_val.size, sample_size))
    return np.dot(c, x)
# dvmax is the ratio between max dv and He, He is in the unit of meter. returned dv is millimeter
def get_sv_short_Zhao2022(dist2wall, dvmax, eta, He):
    sig = 0.46
    dist2wall = dist2wall/He
    dv = 1.14/(dist2wall/eta + 0.39) * 1/sig/np.sqrt(2*np.pi)*np.exp(
        -(np.log(dist2wall/eta+0.39)-0.095)**2/(2*sig**2))*dvmax*He
    return dv
# dvmax is the ratio between max dv and He, He is in the unit of meter. returned dv is millimeter
def get_sh_short_Zhao2022(dist2wall, dlmax, eta, He):
    sig = 0.44
    dist2wall = dist2wall/He
    return 2.14/(dist2wall/eta + 0.82) * 1/sig/np.sqrt(2*np.pi) * np.exp(
        -(np.log(dist2wall/eta+0.82)-0.80)**2/(2*sig**2))*dlmax*He

# mesh of the equivalent beam, vertical disp, longitudinal disp, transverse disp, beam stiffness, beam shear, soil stiffness
def SSIequivalentBeam_Zhao2022(meshX, meshY, dispV, dispL, dispT, Eb, Gb, Es, phi_int, logger):
    #TODO:
    #Hard coded variables, need to be inputed using excavation wrapper
    dfoot = 10
    bfoot = 5
    ni_foot = 0.3

    logger.debug('ASRE_clib: before loading')
    # ASRE_clib = CDLL(r"C:\Users\Jinyan\Desktop\NGI_GIBV\REMEDY_GIS_RiskTool_JZ\ASRETimoshenko.dll")
    ASRE_clib = CDLL(r"C:\Users\Jinyan\Desktop\NGI_GIBV\ASRETimoshenko\ASRETimoshenko\x64\Release\ASRETimoshenko.dll")
    # ASRE_clib = CDLL(r"C:\Users\Jinyan\Desktop\ASRE3DcppDLL\out\build\x64-Release\bin\ASRE3Dlib.dll")
    # ASRE_clib.solve.restype = np.ctypeslib.ndpointer(dtype = np.float64)
    # ASRE_clib.solve.argtypes = [POINTER(c_char),
    #                             np.ctypeslib.ndpointer(dtype=np.float64),
    #                             np.ctypeslib.ndpointer(dtype=np.float64),
    #                             c_int]
    logger.debug('ASRE_clib: '+str(ASRE_clib))
    # ASRE_clib.run.restype = c_int
    # ASRE_clib.run.restype = c_int*3
    # ASRE_clib.run.restype = np.ctypeslib.ndpointer(dtype=c_int, shape=(3,))
    ASRE_clib.run.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.float64), 
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              np.ctypeslib.ndpointer(dtype=np.float64),
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              c_double,
                              POINTER(c_double)
                            #   np.ctypeslib.ndpointer(dtype=np.float64)
                              ]
    # meshX = np.array(meshX)
    # meshY = np.array(meshY)
    result = (c_double*3)() # the output of ASRE_clib is an array with a single element. Use an array in case more output is needed. 
    logger.debug('ASRE_clib result before: '+str(np.array(result))+' ' + str(result[0]))
    status = ASRE_clib.run(meshX, meshY, dispV, dispL, dispT, Eb, Gb, Es, dfoot, bfoot, ni_foot, phi_int, result)
    logger.debug('ASRE_clib result later: '+str(np.array(result)) +' ' + str(result[0]))
    # array_pointer = cast(status, POINTER(ArrayType))
    logger.debug('ASRE_clib: '+str(status))
    logger.debug('ASRE_clib: debug 1')
    libHandle = ASRE_clib._handle
    logger.debug('ASRE_clib: debug 2')
    del ASRE_clib
    logger.debug('ASRE_clib: debug 3')
    _ctypes.FreeLibrary(libHandle)
    logger.debug('ASRE_clib: debug 4')
    return np.array(result)



def near_analysis(xref, yref, construction_area_corners):

    near_dist = 999999
    near_angle = 0  # unit vector

    for corner in construction_area_corners:
        x = corner.x
        y = corner.y

        dist = np.sqrt((xref-x)**2+(yref-y)**2)
        angle = Utils.getAngleFromDir(x, y, xref, yref)

        if dist < near_dist:
            near_dist = dist
            near_angle = angle

    return near_dist, near_angle

def writeBuildingsToShape(buildingsFN, buildings, projection, filterValue, logger):
    fields = [["bid", "long"], ["circumf", "float"], ["foundation", "float"], ["structure", "float"], ["status", "float"],
              ["max_sv_shr", "float"], ["max_sh_shr", "float"], [
                  "max_sv_tot", "float"], ["max_angle", "float"],
              ["max_strain", "float"], ["max_p_stra", "float"], ["vulnerab", "float"], ["risk_tots", "float"], ["risk_angle", "float"], ["filter_val", "string"],
              ["maxStrain", "float"], ["DS", "float"]]

    try:
        numBuildings = 0
        Utils.createShapefile(buildingsFN, "polygon", projection, fields)
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDataSource = shpDriver.Open(buildingsFN, 1)
        outLayer = outDataSource.GetLayer()
        featureDefn = outLayer.GetLayerDefn()

        for building in buildings:
            corners = []
            for corner in building.corners:
                c = [corner.x, corner.y]
                corners.append(c)
            geom = Utils.createPolygon(corners)
            # logger.debug(
            #     "writeBuildingsToShape, polygon created:{}".format(geom))
            outFeature = ogr.Feature(featureDefn)

            # logger.debug(
            #     "writeBuildingsToShape, feature created:{}".format(outFeature))
            outFeature.SetGeometry(geom)
            # logger.debug(
            #     "writeBuildingsToShape, feature created and geometry sat:{}".format(outFeature))
            Utils.addValueToField(
                outFeature, "bid", building.bid, outLayer, logger)
            Utils.addValueToField(outFeature, "circumf",
                                  building.circumf, outLayer, logger)
            Utils.addValueToField(
                outFeature, "foundation", building.foundation, outLayer, logger)
            Utils.addValueToField(outFeature, "structure",
                                  building.structure, outLayer, logger)
            Utils.addValueToField(outFeature, "status",
                                  building.status, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_sv_shr", building.max_sv_short, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_sh_shr", building.max_sh_short, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_sv_tot", building.max_sv_total, outLayer, logger)
            Utils.addValueToField(outFeature, "max_angle",
                                  building.max_angle, outLayer, logger)
            Utils.addValueToField(
                outFeature, "max_strain", building.max_strain, outLayer, logger)
            Utils.addValueToField(outFeature, "max_p_stra",
                                  building.max_principal_strain, outLayer, logger)
            Utils.addValueToField(
                outFeature, "vulnerab", building.vulnerability, outLayer, logger)
            Utils.addValueToField(
                outFeature, "risk_tots", building.risk_totset, outLayer, logger)
            Utils.addValueToField(
                outFeature, "risk_angle", building.risk_angle, outLayer, logger)
            Utils.addValueToField(outFeature, "filter_val",
                                  filterValue, outLayer, logger)
            Utils.addValueToField(outFeature, "maxStrain",
                                  building.maximumStrain, outLayer, logger)
            Utils.addValueToField(outFeature, "DS",
                                  building.damageState, outLayer, logger)

            # logger.debug(
            #     "writeBuildingsToShape, feature created and all attributes added:{}".format(outFeature))
            outLayer.CreateFeature(outFeature)
            # logger.debug(
            #     "writeBuildingsToShape, feature created and added to outLayer:{}".format(outFeature))
            # outLayer = None
            # logger.debug("Building nr: {0}: geom {1}".format(numBuildings, geom.ExportToWkt()))
            outFeature = None
            numBuildings += 1
        # del outFeature
        # Save and close the data source
        outDataSource = None
        logger.debug(
            "Number of buildings written to shape: {}".format(numBuildings))

    except Exception as e:
        logger.error("Error writing buildings: {}".format(type(e)))
        logger.error("Error writing buildings: {}".format(e.args))
        logger.error("Error writing buildings: {}".format(e))
        raise e



def writeWallsToShape(wallsFN, buildings, working_proj, filterValue, logger):

    fields = [["wid", "long"], ["slope_ang", "float"], [
        "h_te_stra", "float"], ["p_te_stra", "float"], ["filter_val", "string"]]

    try:
        Utils.createShapefile(wallsFN, "line", working_proj, fields)
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDataSource = shpDriver.Open(wallsFN, 1)
        outLayer = outDataSource.GetLayer()
        featureDefn = outLayer.GetLayerDefn()
        numWalls = 0
        for building in buildings:
            for wall in building.walls:

                coord = [[wall.corner_start.x, wall.corner_start.y],
                         [wall.corner_end.x, wall.corner_end.y]]

                # logger.debug(
                #     "writeWallsToShape, coord created:{}".format(coord))
                geom = Utils.createLine(coord, logger)
                # logger.debug(
                #     "writeWallsToShape, geom created:{}".format(geom))

                outFeature = ogr.Feature(featureDefn)
                outFeature.SetGeometry(geom)
                Utils.addValueToField(
                    outFeature, "wid", wall.wid, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "filter_val", filterValue, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "slope_ang", wall.slope_angle, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "h_te_stra", wall.hor_tensile_strain, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "p_te_stra", wall.principal_tensile_strain, outLayer, logger)
                outLayer.CreateFeature(outFeature)
                numWalls += 1
                outFeature = None
            # del outFeature
        # Save and close the data source
        outDataSource = None
        logger.debug("Number of walls written to shape: {}".format(numWalls))

    except Exception as e:
        logger.error(type(e))
        logger.error(e.args)
        logger.error(e)
        raise e

def writeCornersToShape(cornersFN, buildings, projection, filterValue, logger):

    #logger.debug(f"Writing corners to shape, Number of buildings:{len(buildings)}, shapefile: {cornersFN}")

    fields = [["cid", "long"], ["dtb", "float"], ["sv_short", "float"], ["sv_long", "float"], ["sv_tot", "float"], [
        "sh_short", "float"], ["near_dist", "float"], ["near_angle", "float"], ["filter_val", "string"], ["porewp_red", "float"]]

    try:
        Utils.createShapefile(cornersFN, "point", projection, fields)
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDataSource = shpDriver.Open(cornersFN, 1)
        outLayer = outDataSource.GetLayer()
        featureDefn = outLayer.GetLayerDefn()

        numCorners = 0
        for building in buildings:
            for corner in building.corners:

                geom = Utils.createPoint(
                    corner.x, corner.y)
                outFeature = ogr.Feature(featureDefn)
                outFeature.SetGeometry(geom)
                Utils.addValueToField(
                    outFeature, "cid", corner.cid, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "dtb", corner.dtb, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "filter_val", filterValue, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sv_short", corner.sv_short, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sv_long", corner.sv_long, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sh_short", corner.sh_short, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "near_dist", corner.near_dist, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "near_angle", corner.near_angle, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "sv_tot", corner.sv_tot, outLayer, logger)
                Utils.addValueToField(
                    outFeature, "porewp_red", corner.porewp_red, outLayer, logger)
                outLayer.CreateFeature(outFeature)
                outFeature = None
                numCorners += 1
            # del outFeature
        # Save and close the data source
        # del outDataSource
        outDataSource = None
        logger.debug(
            "Number of corners written to shape: {}".format(numCorners))
    except Exception as e:
        logger.error(type(e))
        logger.error(e.args)
        logger.error(e)
        raise e

def appendZValuesFromRaster(corner_list, rasterFN, logger=None):
    # if (not logger is None):
    # logger.debug("Logger sent, corner_list: {}".format(corner_list))
    # logger.debug("Logger sent, rasterFN: {}".format(rasterFN))

    src_ds = gdal.Open(str(rasterFN))
    gt = src_ds.GetGeoTransform()
    rb = src_ds.GetRasterBand(1)

    result_corner_list = []

    # numWithoutBr = 0
    for corner_i in corner_list:

        # Convert from map to pixel coordinates.
        # Only works for geotransforms with no rotation.
        px = int((corner_i.x - gt[0]) / gt[1])  # x pixel
        py = int((corner_i.y - gt[3]) / gt[5])  # y pixel

        dtb = rb.ReadAsArray(px, py, 1, 1)

        if (dtb is None):
            # if (not logger is None):
            #     logger.debug("Returning None")
            return None
            # numWithoutBr += 1
            # continue

        if (dtb[0][0] < -50):
            # if (not logger is None):
                # logger.debug("DTB: {}".format(dtb[0][0]))
            return None
        corner_i.dtb = dtb[0][0]
        result_corner_list.append(corner_i)

        # if (numWithoutBr > 0):
        #      logger.debug("Skipping corner in point with no detph to bedrock for {} corners".format(numWithoutBr))
    return result_corner_list


def sampleSv_short_bProbCurve(buildings, dvmaxOverHeRF, etaRF, He, construction_area_corners, n_sample, logger):
    adjacent_buidingID = []
    #collect cornerX and cornerY
    CALCULATION_RANGE=380
    # cornerX = []
    # cornerY = []
    near_dist_corner_ind = []
    near_dist = []
    for bid, bld in enumerate(buildings):
        near_dist_corner_ind_bld = []
        near_dist_bld = []
        # cornerX_bld = []
        # cornerY_bld = []
        adjacent = True
        for corner in bld.corners:
            near_dist_corner_ind_i, near_dist_sqr_i = near_analysis_sqrNcorner(
                corner.x, corner.y, construction_area_corners)
            if near_dist_sqr_i > CALCULATION_RANGE**2:
                corner.sv_shortSample = 0
                bld.adjacent = False
                adjacent = False    
                break   
            # cornerX_bld.append(corner.x)
            # cornerY_bld.append(corner.y)
            near_dist_corner_ind_bld.append(near_dist_corner_ind_i)
            near_dist_bld.append(near_dist_sqr_i)
        if adjacent:
            bld.adjacent = True
            adjacent_buidingID.append(bid)
            # cornerX = cornerX + cornerX_bld
            # cornerY = cornerY + cornerY_bld
            near_dist_corner_ind = near_dist_corner_ind + near_dist_corner_ind_bld
            near_dist = near_dist + near_dist_bld
        # near_dist_corner_ind.append(near_dist_corner_ind_bld)
        # near_dist.append(near_dist_bld)
    tot_corner = len(near_dist_corner_ind)
    log_dvmax_cov = np.zeros([tot_corner, tot_corner])
    log_eta_cov = np.zeros([tot_corner, tot_corner])

    log_dvmax_var = np.log(1+dvmaxOverHeRF.cv**2)
    log_eta_var = np.log(1+etaRF.cv**2)
    log_dvmax_mean = np.log(dvmaxOverHeRF.mean)-log_dvmax_var/2
    log_eta_mean = np.log(etaRF.mean)-log_eta_var/2

    log_dvmax_nug = dvmaxOverHeRF.nuggetOverVar*log_dvmax_var
    log_eta_nug = etaRF.nuggetOverVar*log_eta_var

    for corner_i in range(tot_corner):
        for corner_j in range(tot_corner):
            if corner_i == corner_j:
                log_dvmax_cov[corner_i, corner_j] = log_dvmax_var
                log_eta_cov[corner_i, corner_j] = log_eta_var
            else:
                construction_corner_i = construction_area_corners[near_dist_corner_ind[corner_i]]
                construction_corner_j = construction_area_corners[near_dist_corner_ind[corner_j]]
                dist_ij = np.sqrt((construction_corner_i.x - construction_corner_j.x)**2 + 
                                    (construction_corner_i.y - construction_corner_j.y)**2)
                log_dvmax_cov[corner_i, corner_j] = gaussianAutoCorr(dist_ij, (log_dvmax_var-log_dvmax_nug), log_dvmax_nug, dvmaxOverHeRF.range)*log_dvmax_var + log_dvmax_nug
                log_eta_cov[corner_i, corner_j] = gaussianAutoCorr(dist_ij, (log_eta_var-log_eta_nug), log_eta_nug, etaRF.range)*log_eta_var + log_eta_nug
    logger.debug("Covariance matrix computed")
    ######## Generate samples for log_dvmax and log_eta
    log_dvmax_eig_val, log_dvmax_eig_vec = scipy.linalg.eigh(log_dvmax_cov)
    log_dvmax_eig_val = np.flip(log_dvmax_eig_val)
    log_dvmax_eig_vec = np.flip(log_dvmax_eig_vec, axis = 1)

    log_eta_eig_val, log_eta_eig_vec = scipy.linalg.eigh(log_eta_cov)
    log_eta_eig_val = np.flip(log_eta_eig_val)
    log_eta_eig_vec = np.flip(log_eta_eig_vec, axis = 1)

    log_eta_cov = None # Release the memory used to store log_eta_cov
    log_dvmax_cov = None # Release the memory used to store log_eta_cov
    
    # Find 99% explained variance
    eig_val_accu = 0.
    eig_val_sum = log_dvmax_eig_val.sum()
    for i in range(0, log_dvmax_eig_val.size):
        eig_val_accu += log_dvmax_eig_val[i]
        if eig_val_accu/eig_val_sum > 0.99:
            break
    dvmax_num_eig = i+1

    eig_val_accu = 0.
    eig_val_sum = log_eta_eig_val.sum()
    for i in range(0, log_eta_eig_val.size):
        eig_val_accu += log_eta_eig_val[i]
        if eig_val_accu/eig_val_sum > 0.99:
            break
    eta_num_eig = i+1

    log_dvmax_samples = generateCorrGauSample(n_sample, log_dvmax_eig_vec[:, 0:dvmax_num_eig], log_dvmax_eig_val[0:dvmax_num_eig]) + log_dvmax_mean
    log_eta_samples = generateCorrGauSample(n_sample, log_eta_eig_vec[:, 0:eta_num_eig], log_eta_eig_val[0:eta_num_eig]) + log_eta_mean
    # log_dvmax_samples is corener_txt by n_sample matrix
    dvmax_samples = np.exp(log_dvmax_samples)
    logger.debug("median of dvmaxOverHe at corner 1: "+str(np.mean(dvmax_samples[0,:])))
    eta_samples = np.exp(log_eta_samples)
    logger.debug("median of eta at corner 1: "+str(np.mean(eta_samples[0,:])))
    logger.debug("Samples generated")
    # Release some useless memory
    log_dvmax_eig_val = None
    log_dvmax_eig_vec = None
    log_eta_eig_val = None
    log_eta_eig_vec = None
    log_dvmax_samples = None
    log_eta_samples = None
    ################# Calculate dv
    short_dv_results = np.zeros([tot_corner, n_sample], dtype=np.float32)
    for ind in range(0, tot_corner):
        dist = near_dist[ind]
        if dist > CALCULATION_RANGE**2:
            # px += 1
            # new_progress = int(100 * px / tot_px)
            # if new_progress > progress:
            #     progress = new_progress
            #     logger.info("Progress: " + str(progress) + " %. Outside.")
            continue
        for i_sample in range(0, n_sample):
            short_dv_results[ind, i_sample] = get_sv_short_Zhao2022(dist, dvmax_samples[ind, i_sample], eta_samples[ind, i_sample], He)
        # px += 1
        # new_progress = int(100*px/tot_px)
        # if  new_progress > progress:
        #     progress = new_progress
        #     logger.info("Progress: " + str(progress) + " %")
    logger.debug("short_dv_results: "+str(short_dv_results))
    ################# Save short_dv_results to corner.sv_shortSample
    ind = 0
    for bld_id in adjacent_buidingID:
        bld = buildings[bld_id]
        for corner in bld.corners:
            corner.sv_shortSample = short_dv_results[ind,:]
            logger.debug("Corner {0} shor_dv_sample is {1}: ".format(ind, str(corner.sv_shortSample)))
            ind = ind + 1
    logger.debug("Short term adjacent_buidingID: "+str(adjacent_buidingID))
    logger.info("Short term settlement sampled at all building corners")



def getPrincipleAxis(building_polygon):
    BUILDING_PIXEL_SIZE = 1
    building_raster = "temp_building_raster"
    arcpy.PolyganoToRaster_conversion(building_polygon, 'building_ID', building_raster, "CELL_CENTER","", BUILDING_PIXEL_SIZE)
    poly_proj = Utils_arcpy.getProjCodeFromFC(building_polygon)
    building_raster_proj = Utils_arcpy.getProjCodeFromFC(building_raster)
    desc_building_ras = arcpy.Describe(building_raster)
    cellsize_building_raster = desc_building_ras.meanCellWidth
    extent_building = desc_building_ras.extent
    xmin_buil = extent_building.XMin
    ymin_buil = extent_building.YMin
    xmax_buil = extent_building.XMax
    ymax_buil = extent_building.YMax

    xlen_buil = xmax_buil - xmin_buil
    ylen_buil = ymax_buil - ymin_buil
    
    builRasterArray = arcpy.RasterToNumPyArray(building_raster, nodata_to_value=0)
    nrows_buil, ncols_buil = builRasterArray.shape
    dx_buil = xlen_buil / ncols_buil
    dy_buil = ylen_buil / nrows_buil

    building_pixel_coor = []
    #Assembling the data structure "building_pixel_dict" that contains the coordinates at of each pixel for a particular building
    for col in range(0, ncols_buil - 1):
        for row in range(0, nrows_buil - 1):

            auto_BID = int(builRasterArray.item(row, col))  # Building ID is the value of the cell

            #bid_test2 = builRasterArray.bygning_ID
            #arcpy.AddMessage("test bid: " + str(bid_test2))

            x = xmin_buil + dx_buil * (col)  # Coordinates of the pixel
            y = ymin_buil + dy_buil * (nrows_buil - row)

            if auto_BID == 0:
                #Outside a building
                continue
            # elif auto_BID in building_pixel_dict.keys():
            #     building_pixel_dict[auto_BID].append([x, y])
            else:
                building_pixel_coor.append([x, y])
    count = 0
    # arcpy.CreateFeatureclass_management(output_ws, "temp_dir_lines", "POLYLINE")
    # dir_lines = output_ws + "\\" + "temp_dir_lines"
    # field_names = ['SHAPE@']
    # dir_lines_cursor = arcpy.da.InsertCursor(dir_lines, field_names)

    Ix = 0
    Iy = 0
    for pixel in building_pixel_coor:
        x = pixel[0]
        y = pixel[1]
        Ix += x
        Iy += y

    centeroid_x = Ix / len(building_pixel_coor)
    centeroid_y = Iy / len(building_pixel_coor)

    #Second area moment components:
    Ixx = 0
    Iyy = 0
    Ixy = 0 #symmetric matrix, Ixy = Iyx

    for pixel in building_pixel_coor:
        x = pixel[0]
        y = pixel[1]
        cx = x - centeroid_x
        cy = y - centeroid_y
        Ixx += cx*cx
        Iyy += cy*cy
        Ixy += cx*cy

    I = np.array([[Ixx, Ixy], [Ixy, Iyy]])
    eigenvalues, eigenvectors = la.eig(I)

    dir_x = 0
    dir_y = 0
    dirlen = 0
    if eigenvalues[0]> eigenvalues[1] and eigenvalues[1] > 0.001:
        dir_x, dir_y = eigenvectors[:, 0] #unit vector
        dir1 = eigenvectors[:, 0]
        dir2 = eigenvectors[:, 1]
        dirlen = eigenvalues[0]/eigenvalues[1]
    elif eigenvalues[1]> eigenvalues[0] and eigenvalues[0] > 0.001:
        dir_x, dir_y = eigenvectors[:, 1] #unit vector
        dir1 = eigenvectors[:, 1]
        dir2 = eigenvectors[:, 0]
        dirlen = eigenvalues[1]/eigenvalues[0]
    else:
        dir_x = 1
        dir_y = 1
        dirlen = 0

    #loop to find "length" of building along dominant direction
    step = 0.5
    dist1 = step
    dist2 = step
    cont = True
    security = 0

    try:
        while cont and security < 10000:
            cont = False
            for pixel in building_pixel_coor:

                x = pixel[0]
                y = pixel[1]

                line1x = centeroid_x + dir_x*dist1
                line1y = centeroid_y + dir_y*dist1
                line2x = centeroid_x + dir_x*dist2
                line2y = centeroid_y + dir_y*dist2

                if abs(line1x - x) < 1.5*step and abs(line1y - y) < 1.5*step:

                    dist1 += step
                    cont = True
                    #arcpy.AddMessage("HELLOOO1")
                    break

                if abs(line2x - x) < 1.5*step and abs(line2y - y) < 1.5*step:
                    dist2 += step
                    cont = True
                    #arcpy.AddMessage("HELLOOO2")
                    break

            security += 1

    except:
        dist1 = 0
        dist2 = 0

    axis1end1 = [line1x, line1y]
    axis1end2 = [line2x, line2y]
    dominant_buil_len_axis1 = dist1 + dist2

    step = 0.5
    dist1 = step
    dist2 = step
    cont = True
    security = 0
    dir_x = dir2[0]
    dir_y = dir2[1]

    try:
        while cont and security < 10000:
            cont = False
            for pixel in building_pixel_coor:

                x = pixel[0]
                y = pixel[1]

                line1x = centeroid_x + dir_x*dist1
                line1y = centeroid_y + dir_y*dist1
                line2x = centeroid_x + dir_x*dist2
                line2y = centeroid_y + dir_y*dist2

                if abs(line1x - x) < 1.5*step and abs(line1y - y) < 1.5*step:

                    dist1 += step
                    cont = True
                    #arcpy.AddMessage("HELLOOO1")
                    break

                if abs(line2x - x) < 1.5*step and abs(line2y - y) < 1.5*step:
                    dist2 += step
                    cont = True
                    #arcpy.AddMessage("HELLOOO2")
                    break

            security += 1

    except:
        dist1 = 0
        dist2 = 0

    axis2end1 = [line1x, line1y]
    axis2end2 = [line2x, line2y]
    dominant_buil_len_axis2 = dist1 + dist2




def meshBuildings(buildings, beamElemSize, logger):
    """This function call equivalentBeamMesh function to all
    bld in the buildings. The axis(1or2)meash(XorY) fields
    in buildings are initialized.

    Args:
        buildings: A list of building object.
        beamElemSize: The approximate size of beam elements.
        logger: A logger to write to

    Returns:
        Void.
    """
    for bld in buildings:
        bld.equivalentBeamMesh(beamElemSize, logger)
# def sample_gfDisp_short_bProbSSI(buildings, dvmaxRF, etaRF, excavation_depth, construction_area_corners, MCsampleSize, logger):
#     #TODO: loop and register near buildings
#     nodeX1 = []
#     nodeY1 = []
#     near_dist_corner_ind1 = []
#     near_dist1 = []
#     nodeX2 = []
#     nodeY2 = []
#     near_dist_corner_ind2 = []
#     near_dist2 = []
#     for bld in buildings:
#         bld.equivalentBeamMesh(1, logger)
#         nodeX1.append(bld.axis1meshX)
#         nodeX2.append(bld.axis2meshX)
#         nodeY1.append(bld.axis1meshY)
#         nodeY2.append(bld.axis2meshY)

def sample_gfDispNEs_short_bProbSSI_singleAxis(buildings, axis, dvmaxOverHeRF, etaRF, dloverdvRV, EsRF, construction_area_corners, n_sample, logger):
    adjacent_buidingID = []
    CALCULATION_RANGE=380
    near_dist_corner_ind = []
    near_dist = []
    for bid, bld in enumerate(buildings):
        if axis==1:
            nnode_bld = len(bld.axis1meshX)
            axisMeshX_bld = bld.axis1meshX
            axisMeshY_bld = bld.axis1meshY
        else:
            nnode_bld = len(bld.axis2meshX)
            axisMeshX_bld = bld.axis2meshX
            axisMeshY_bld = bld.axis2meshY
        near_dist_corner_ind_bld = []
        near_dist_bld = []
        adjacent = True
        for ind in range(nnode_bld):
            near_dist_corner_ind_i, near_dist_sqr_i = near_analysis_sqrNcorner(
                axisMeshX_bld[ind], axisMeshY_bld[ind], construction_area_corners)
            if near_dist_sqr_i > CALCULATION_RANGE**2:
                bld.adjacent = False
                adjacent = False    
                break   
            # axis1X_bld.append(axisMeshX[ind])
            # axis1Y_bld.append(axisMeshY[ind])
            near_dist_corner_ind_bld.append(near_dist_corner_ind_i)
            near_dist_bld.append(near_dist_sqr_i)
        if adjacent:
            bld.adjacent = True
            adjacent_buidingID.append(bid)
            # axis1X = axis1X + axis1X_bld
            # axis1Y = axis1Y + axis1Y_bld
            near_dist_corner_ind = near_dist_corner_ind + near_dist_corner_ind_bld
            near_dist = near_dist + near_dist_bld
    tot_nodes = len(near_dist)
    log_dvmax_cov = np.zeros([tot_nodes, tot_nodes])
    log_eta_cov = np.zeros([tot_nodes, tot_nodes])
    log_Es_cov = np.zeros([tot_nodes, tot_nodes])

    log_dvmax_var = np.log(1+dvmaxOverHeRF.cv**2)
    log_eta_var = np.log(1+etaRF.cv**2)
    log_Es_var = np.log(1+EsRF.cv**2)
    log_dvmax_mean = np.log(dvmaxOverHeRF.mean)-log_dvmax_var/2
    log_eta_mean = np.log(etaRF.mean)-log_eta_var/2
    log_Es_mean = np.log(EsRF.mean)-log_Es_var/2
    

    log_dvmax_nug = dvmaxOverHeRF.nuggetOverVar*log_dvmax_var
    log_eta_nug = etaRF.nuggetOverVar*log_eta_var
    log_Es_nug = EsRF.nuggetOverVar*log_Es_var

    for corner_i in range(tot_nodes):
        for corner_j in range(tot_nodes):
            if corner_i == corner_j:
                log_dvmax_cov[corner_i, corner_j] = log_dvmax_var
                log_eta_cov[corner_i, corner_j] = log_eta_var
                log_Es_cov[corner_i, corner_j] = log_Es_var
            else:
                construction_corner_i = construction_area_corners[near_dist_corner_ind[corner_i]]
                construction_corner_j = construction_area_corners[near_dist_corner_ind[corner_j]]
                dist_ij = np.sqrt((construction_corner_i.x - construction_corner_j.x)**2 + 
                                    (construction_corner_i.y - construction_corner_j.y)**2)
                log_dvmax_cov[corner_i, corner_j] = gaussianAutoCorr(dist_ij, (log_dvmax_var-log_dvmax_nug), log_dvmax_nug, dvmaxOverHeRF.range)*log_dvmax_var + log_dvmax_nug
                log_eta_cov[corner_i, corner_j] = gaussianAutoCorr(dist_ij, (log_eta_var-log_eta_nug), log_eta_nug, etaRF.range)*log_eta_var + log_eta_nug
                log_eta_cov[corner_i, corner_j] = gaussianAutoCorr(dist_ij, (log_Es_var-log_Es_nug), log_Es_nug, EsRF.range)*log_Es_var + log_Es_nug
    logger.debug("Covariance matrix for axis 1 computed")
    ######## Generate samples for log_dvmax and log_eta
    log_dvmax_eig_val, log_dvmax_eig_vec = scipy.linalg.eigh(log_dvmax_cov)
    log_dvmax_eig_val = np.flip(log_dvmax_eig_val)
    log_dvmax_eig_vec = np.flip(log_dvmax_eig_vec, axis = 1)

    log_eta_eig_val, log_eta_eig_vec = scipy.linalg.eigh(log_eta_cov)
    log_eta_eig_val = np.flip(log_eta_eig_val)
    log_eta_eig_vec = np.flip(log_eta_eig_vec, axis = 1)

    log_Es_eig_val, log_Es_eig_vec = scipy.linalg.eigh(log_Es_cov)
    log_Es_eig_val = np.flip(log_Es_eig_val)
    log_Es_eig_vec = np.flip(log_Es_eig_vec, axis = 1)

    log_eta_cov = None # Release the memory used to store log_eta_cov
    log_dvmax_cov = None # Release the memory used to store log_dvmax_cov
    log_Es_cov = None # Release the memory used to store log_Es_cov
    
    # Find 99% explained variance
    eig_val_accu = 0.
    eig_val_sum = log_dvmax_eig_val.sum()
    for i in range(0, log_dvmax_eig_val.size):
        eig_val_accu += log_dvmax_eig_val[i]
        if eig_val_accu/eig_val_sum > 0.99:
            break
    dvmax_num_eig = i+1

    eig_val_accu = 0.
    eig_val_sum = log_eta_eig_val.sum()
    for i in range(0, log_eta_eig_val.size):
        eig_val_accu += log_eta_eig_val[i]
        if eig_val_accu/eig_val_sum > 0.99:
            break
    eta_num_eig = i+1

    eig_val_accu = 0.
    eig_val_sum = log_Es_eig_val.sum()
    for i in range(0, log_Es_eig_val.size):
        eig_val_accu += log_Es_eig_val[i]
        if eig_val_accu/eig_val_sum > 0.99:
            break
    Es_num_eig = i+1

    log_dvmax_samples = generateCorrGauSample(n_sample, log_dvmax_eig_vec[:, 0:dvmax_num_eig], log_dvmax_eig_val[0:dvmax_num_eig]) + log_dvmax_mean
    log_eta_samples = generateCorrGauSample(n_sample, log_eta_eig_vec[:, 0:eta_num_eig], log_eta_eig_val[0:eta_num_eig]) + log_eta_mean
    log_Es_samples = generateCorrGauSample(n_sample, log_Es_eig_vec[:, 0:Es_num_eig], log_Es_eig_val[0:Es_num_eig]) + log_Es_mean
    # log_dvmax_samples is corener_txt by n_sample matrix
    dvmax_samples = np.exp(log_dvmax_samples)
    logger.debug("median of dvmaxOverHe at bld 1, node 1 is {0}, at bld1, node 2 is {1}: ".format(str(np.mean(dvmax_samples[0,:])), str(np.mean(dvmax_samples[1,:]))))
    eta_samples = np.exp(log_eta_samples)
    logger.debug("median of eta at bld 1, node 1 is {0}, at bld1, node 2 is {1}: ".format(str(np.mean(eta_samples[0,:])), str(np.mean(eta_samples[1,:]))))
    Es_samples = np.exp(log_Es_samples)
    logger.debug("median of Es at bld 1, node 1 is {0}, at bld1, node 2 is {1}: ".format(str(np.mean(Es_samples[0,:])), str(np.mean(Es_samples[1,:]))))
    logger.debug("Samples generated")
    # Release some useless memory
    log_dvmax_eig_val = None
    log_dvmax_eig_vec = None
    log_eta_eig_val = None
    log_eta_eig_vec = None
    log_Es_eig_val = None
    log_Es_eig_vec = None
    log_dvmax_samples = None
    log_eta_samples = None
    log_Es_samples = None
    ################# Save short_dv_results to corner.sv_shortSample
    count = 0
    for bld_ind in adjacent_buidingID:
        bld = buildings[bld_ind]
        if axis==1:
            nnode_bld = len(bld.axis1meshX)
        else:
            nnode_bld = len(bld.axis2meshX)
        bld_dvmaxOverHe_Sample = dvmax_samples[count:count+nnode_bld,:]
        bld_etaSample = eta_samples[count:count+nnode_bld,:]
        bld_EsSample = Es_samples[count:count+nnode_bld,:]
        bld_dlmaxOverdvmaxSample = np.random.normal(dloverdvRV.mean, dloverdvRV.mean * dloverdvRV.cv, n_sample)
        if axis==1:
            bld.dvmaxOverHe_Sample1 = bld_dvmaxOverHe_Sample
            bld.eta_Sample1 = bld_etaSample
            bld.Es_Sample1 = bld_EsSample
            bld.dlmaxOverdvmaxSample1 = bld_dlmaxOverdvmaxSample
        else:
            bld.dvmaxOverHe_Sample2 = bld_dvmaxOverHe_Sample
            bld.eta_Sample2 = bld_etaSample
            bld.Es_Sample2 = bld_EsSample
            bld.dlmaxOverdvmaxSample2 = bld_dlmaxOverdvmaxSample
        count = count+nnode_bld
    logger.info("Short term settlement sampled at all building corners")


def greenFieldDisp_zhao2022_determine(axis1meshX, axis1meshY, axis2meshX, axis2meshY, dvmaxOverHe, eta, excavation_depth, dlmaxOverdvmax, construction_area_corners, logger):
    #TODO 1: calculate dispX, dispY, dispZ for axis 1 and axis 2
    #TODO 2: calculate dispV, dispL, dispT for axis 1 and axis 2
        dispV1 = np.zeros(len(axis1meshX))
        dispL1 = np.zeros(len(axis1meshX))
        dispT1 = np.zeros(len(axis1meshX))
        dispV2 = np.zeros(len(axis2meshX))
        dispL2 = np.zeros(len(axis2meshX))
        dispT2 = np.zeros(len(axis2meshX))
        # For axis 1
        near_dist_corner1 = np.sqrt(near_analysis_sqr(axis1meshX[ 0], axis1meshY[ 0], construction_area_corners))
        near_dist_corner2 = np.sqrt(near_analysis_sqr(axis1meshX[-1], axis1meshY[-1], construction_area_corners))
        if near_dist_corner2 < near_dist_corner1:
            vecAxis1 = np.array([axis1meshX[0] - axis1meshX[-1], axis1meshY[0] - axis1meshY[-1]])
        else:
            # axis1_angle = Utils.getAngleFromDir(self.axis1meshX[-1], self.axis1meshY[-1],
            #                                         self.axis1meshX[0], self.axis1meshY[0])
            vecAxis1 = np.array([axis1meshX[-1] - axis1meshX[0], axis1meshY[-1] - axis1meshY[0]])
        #rotate vecAxis1 clockwise by 90 degree: https://limnu.com/sketch-easy-90-degree-rotate-vectors/#:~:text=Normally%20rotating%20vectors%20involves%20matrix,swap%20X%20and%20Y%20values.
        vecAxis1_perpend = np.array([vecAxis1[1], -vecAxis1[0]])
        for i in range(len(axis1meshX)):
            # near_dist, near_angle = near_analysis(
            #     self.axis1meshX[i], self.axis1meshY[i], construction_area_corners)
            near_point_ind, near_dist_sqr = near_analysis_sqrNcorner(
                axis1meshX[i], axis1meshY[i], construction_area_corners)
            near_dist = np.sqrt(near_dist_sqr)
            #near_angle is the angle from east to near_point using building corner as origin
            dispV1[i]=get_sv_short_Zhao2022(near_dist, dvmaxOverHe, eta, excavation_depth)
            dispH = get_sh_short_Zhao2022(near_dist, dlmaxOverdvmax*dvmaxOverHe, eta, excavation_depth)
            vecdH = np.array([construction_area_corners[near_point_ind].x-axis1meshX[i], 
                             construction_area_corners[near_point_ind].y-axis1meshY[i]])
            dispL1[i] = dispH/near_dist * np.dot(vecdH, vecAxis1)/np.linalg.norm(vecAxis1)
            dispT1[i] = dispH/near_dist * np.dot(vecdH, vecAxis1_perpend)/np.linalg.norm(vecAxis1_perpend)
        
        # For axis 2
        near_dist_corner1 = np.sqrt(near_analysis_sqr(axis2meshX[ 0], axis2meshY[ 0], construction_area_corners))
        near_dist_corner2 = np.sqrt(near_analysis_sqr(axis2meshX[-1], axis2meshY[-1], construction_area_corners))
        #axis1_angle is the angel from north to building axis using the axis end closer to excavation as origin, unit is degree
        if near_dist_corner2 < near_dist_corner1:
            vecAxis2 = np.array([axis2meshX[0] - axis2meshX[-1], axis2meshY[0] - axis2meshY[-1]])
        else:
            vecAxis2 = np.array([axis2meshX[-1] - axis2meshX[0], axis2meshY[-1] - axis2meshY[0]])
        #rotate vecAxis2 clockwise by 90 degree: https://limnu.com/sketch-easy-90-degree-rotate-vectors/#:~:text=Normally%20rotating%20vectors%20involves%20matrix,swap%20X%20and%20Y%20values.
        vecAxis2_perpend = np.array([vecAxis2[1], -vecAxis2[0]])
        for i in range(len(axis2meshX)):
            near_point_ind, near_dist_sqr = near_analysis_sqrNcorner(
                axis2meshX[i], axis2meshY[i], construction_area_corners)
            near_dist = np.sqrt(near_dist_sqr)
            #near_angle is the angle from east to near_point using building corner as origin
            dispV2[i]=get_sv_short_Zhao2022(near_dist, dvmaxOverHe, eta, excavation_depth)
            dispH = get_sh_short_Zhao2022(near_dist, dlmaxOverdvmax*dvmaxOverHe, eta, excavation_depth)
            vecdH = np.array([construction_area_corners[near_point_ind].x-axis2meshX[i], 
                             construction_area_corners[near_point_ind].y-axis2meshY[i]])
            dispL2[i] = dispH/near_dist * np.dot(vecdH,vecAxis2)/np.linalg.norm(vecAxis2)
            dispT2[i] = dispH/near_dist * np.dot(vecdH,vecAxis2_perpend)/np.linalg.norm(vecAxis2_perpend)
        logger.debug('dispL1 is: '+str(dispL1))
        return [dispL1, dispT1, dispV1], [dispL2, dispT2, dispV2]