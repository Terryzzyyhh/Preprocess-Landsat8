import os
from osgeo import gdal
import numpy as np
import shutil
import tarfile
import argparse
import glob
from Py6S import *
import re

# 逐波段辐射定标
def RadiometricCalibration(data2,ImgRasterData,BandId):
    # LandSat8 TM辐射定标参数
    #global data2, ImgRasterData
    parameter_OLI = np.zeros((9, 2))

    # 计算辐射亮度参数
    parameter_OLI[0, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_1.+', data2)[0]).split("=")[1])
    parameter_OLI[1, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_2.+', data2)).split("=")[1])
    parameter_OLI[2, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_3.+', data2)).split("=")[1])
    parameter_OLI[3, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_4.+', data2)).split("=")[1])
    parameter_OLI[4, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_5.+', data2)).split("=")[1])
    parameter_OLI[5, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_6.+', data2)).split("=")[1])
    parameter_OLI[6, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_7.+', data2)).split("=")[1])
    parameter_OLI[7, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_8.+', data2)).split("=")[1])
    parameter_OLI[8, 0] = float(''.join(re.findall('RADIANCE_MULT_BAND_9.+', data2)).split("=")[1])

    parameter_OLI[0, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_1.+', data2)[0]).split("=")[1])
    parameter_OLI[1, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_2.+', data2)).split("=")[1])
    parameter_OLI[2, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_3.+', data2)).split("=")[1])
    parameter_OLI[3, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_4.+', data2)).split("=")[1])
    parameter_OLI[4, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_5.+', data2)).split("=")[1])
    parameter_OLI[5, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_6.+', data2)).split("=")[1])
    parameter_OLI[6, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_7.+', data2)).split("=")[1])
    parameter_OLI[7, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_8.+', data2)).split("=")[1])
    parameter_OLI[8, 1] = float(''.join(re.findall('RADIANCE_ADD_BAND_9.+', data2)).split("=")[1])

    Gain = parameter_OLI[int(BandId) - 1, 0]
    Bias = parameter_OLI[int(BandId) - 1, 1]

    RaCal = np.where(ImgRasterData > 0, Gain * ImgRasterData + Bias, 0)
    #print("最大最小值：",np.max(RaCal),np.min(RaCal))
    return (RaCal)

# 6s大气校正
def AtmosphericCorrection(data2,BandId):
    global data
    # 6S模型
    s = SixS()

    s.geometry = Geometry.User()
    s.geometry.solar_z = 90 - float(''.join(re.findall('SUN_ELEVATION.+', data2)).split("=")[1])
    s.geometry.solar_a = float(''.join(re.findall('SUN_AZIMUTH.+', data2)).split("=")[1])
    s.geometry.view_z = 0
    s.geometry.view_a = 0

    # 日期
    Dateparm = ''.join(re.findall('DATE_ACQUIRED.+', data2)).split("=")
    Date = Dateparm[1].split('-')

    s.geometry.month = int(Date[1])
    s.geometry.day = int(Date[2])

    # 中心经纬度
    point1lat = float(''.join(re.findall('CORNER_UL_LAT_PRODUCT.+', data2)).split("=")[1])
    point1lon = float(''.join(re.findall('CORNER_UL_LON_PRODUCT.+', data2)).split("=")[1])
    point2lat = float(''.join(re.findall('CORNER_UR_LAT_PRODUCT.+', data2)).split("=")[1])
    point2lon = float(''.join(re.findall('CORNER_UR_LON_PRODUCT.+', data2)).split("=")[1])
    point3lat = float(''.join(re.findall('CORNER_LL_LAT_PRODUCT.+', data2)).split("=")[1])
    point3lon = float(''.join(re.findall('CORNER_LL_LON_PRODUCT.+', data2)).split("=")[1])
    point4lat = float(''.join(re.findall('CORNER_LR_LAT_PRODUCT.+', data2)).split("=")[1])
    point4lon = float(''.join(re.findall('CORNER_LR_LON_PRODUCT.+', data2)).split("=")[1])

    sLongitude = (point1lon + point2lon + point3lon + point4lon) / 4
    sLatitude = (point1lat + point2lat + point3lat + point4lat) / 4

    # 大气模式类型
    if sLatitude > -15 and sLatitude <= 15:
        s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)

    if sLatitude > 15 and sLatitude <= 45:
        if s.geometry.month > 4 and s.geometry.month <= 9:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
        else:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)

    if sLatitude > 45 and sLatitude <= 60:
        if s.geometry.month > 4 and s.geometry.month <= 9:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticSummer)
        else:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticWinter)

    # 气溶胶类型大陆
    s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

    # 目标地物？？？？？？
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    # 550nm气溶胶光学厚度,根据日期从MODIS处获取。
    s.visibility=40.0
    s.aot550 = 0.14497

    # 通过研究去区的范围去求DEM高度。
    pointUL = dict()
    pointDR = dict()
    pointUL["lat"] = point1lat
    pointUL["lon"] = point1lon
    pointDR["lat"] = point4lat
    pointDR["lon"] = point2lon
    meanDEM = (MeanDEM(pointUL, pointDR)) * 0.001

    # 研究区海拔、卫星传感器轨道高度
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(meanDEM)
    s.altitudes.set_sensor_satellite_level()

    # 校正波段（根据波段名称）
    if BandId == '1':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B1)

    elif BandId == '2':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B2)

    elif BandId == '3':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B3)

    elif BandId == '4':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B4)

    elif BandId == '5':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B5)

    elif BandId == '6':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B6)

    elif BandId == '7':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B7)

    elif BandId == '8':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B8)

    elif BandId == '9':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B9)

    # 下垫面非均一、朗伯体
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)

    # 运行6s大气模型
    s.run()

    xa = s.outputs.coef_xa
    xb = s.outputs.coef_xb
    xc = s.outputs.coef_xc
    x = s.outputs.values
    return (xa, xb, xc)
def MeanDEM(pointUL, pointDR):
    '''
    计算影像所在区域的平均高程.
    '''
    script_path = os.path.split(os.path.realpath(__file__))[0]
    dem_path = os.path.join(script_path, "GMTED2km.tif")

    try:
        DEMIDataSet = gdal.Open(dem_path)
    except Exception as e:
        pass

    DEMBand = DEMIDataSet.GetRasterBand(1)
    geotransform = DEMIDataSet.GetGeoTransform()
    # DEM分辨率
    pixelWidth = geotransform[1]
    pixelHight = geotransform[5]

    # DEM起始点：左上角，X：经度，Y：纬度
    originX = geotransform[0]
    originY = geotransform[3]

    # 研究区左上角在DEM矩阵中的位置
    yoffset1 = int((originY - pointUL['lat']) / pixelWidth)
    xoffset1 = int((pointUL['lon'] - originX) / (-pixelHight))

    # 研究区右下角在DEM矩阵中的位置
    yoffset2 = int((originY - pointDR['lat']) / pixelWidth)
    xoffset2 = int((pointDR['lon'] - originX) / (-pixelHight))

    # 研究区矩阵行列数
    xx = xoffset2 - xoffset1
    yy = yoffset2 - yoffset1
    # 读取研究区内的数据，并计算高程
    DEMRasterData = DEMBand.ReadAsArray(xoffset1, yoffset1, xx, yy)

    MeanAltitude = np.mean(DEMRasterData)
    return MeanAltitude

def makedir(dirstr):
    if (os.path.exists(dirstr)):
        shutil.rmtree(dirstr)
    os.makedirs(dirstr)


# 解压缩原始文件
def untar(fname, dirs):
    """
    解压tar.gz文件
    :param fname: 压缩文件名
    :param dirs: 解压后的存放路径
    :return: bool
    """
    try:
        t = tarfile.open(fname)
        t.extractall(path = dirs)
        return True
    except Exception as e:
        print(e)
        return False

def parse_arguments(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('--Input_dir', type=str, help='Input dir', default=None)
    parser.add_argument('--Output_dir', type=str, help='Output dir', default=None)

    return parser.parse_args(argv)

def Get_tifinfor(indir):
    inputs=glob.glob(os.path.join(indir,"*tiff"))
    bandsum=7
    dataset=gdal.Open(inputs[0])
    outprj=dataset.GetProjection()
    outtfm=dataset.GetGeoTransform()
    img_w=dataset.RasterXSize
    img_h=dataset.RasterYSize
    #print("波段总数：",bandsum,"\n","x:",img_w,"\n","y:",img_h)
    return bandsum,outtfm,outprj,img_w,img_h

def Stackband(indir,savepath):

    bandsum, outtfm, outprj, img_w, img_h=Get_tifinfor(indir)
    driver = gdal.GetDriverByName("GTiff")
    outdataset = driver.Create(savepath, img_w, img_h, bandsum, gdal.GDT_Int16)
    outdataset.SetProjection(outprj)
    outdataset.SetGeoTransform(outtfm)

    inputs=glob.glob(os.path.join(indir,"*tiff"))
    for i in range(bandsum):
        dataset=gdal.Open(inputs[i])
        banddata=dataset.GetRasterBand(1).ReadAsArray()
        outdataset.GetRasterBand(i+1).WriteArray(banddata)

def Pansharpentif(multitifpath,pantifpath):
    # os.chdir(tiffdir)
    # pantifname = glob.glob("*PAN*.tiff")
    # pantifpath = os.path.join(tiffdir, pantifname[0])
    # multitifname = glob.glob("*AC*.tiff")
    # multitifpath = os.path.join(tiffdir, multitifname[0])

    pssavepath= multitifpath.replace('.tiff', "_PanSharpen.tiff")

    "D:\Anaconda3\python.exe"
    "D:\Anaconda3\Scripts\gdal_pansharpen.py"
    str = (
            'D:\\Anaconda3\python D:\\Anaconda3\Scripts\gdal_pansharpen.py -bitdepth 16 -nodata 0 -of GTiff %s %s %s' % (
        pantifpath, multitifpath, pssavepath))
    os.system(str)


