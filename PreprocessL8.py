
import numpy
from base import *
from multiprocessing import Pool
import glob

def ProcessBlock(tar):
    tardir=os.path.dirname(tar)
    unzipdir = os.path.join(tardir, "unzip")
    resultdir = os.path.join(tardir, "result")
    makedir(unzipdir)
    makedir(resultdir)
    
    # 输入数据路径
    # RootInputPath = parse_arguments(sys.argv[1:]).Input_dir
    # RootOutName = parse_arguments(sys.argv[2:]).Output_dir
    name=os.path.basename(tar)
    print("解压 ",name)
    name=os.path.splitext(name)[0]
    name = os.path.splitext(name)[0]
    RootInputPath=os.path.join(unzipdir,name)
    makedir(RootInputPath)
    untar(tar,RootInputPath)
    print("----解压完成----")
    RootOutName=os.path.join(resultdir,name,"temp")
    makedir(RootOutName)

    # 创建日志文件
    LogFile = open(os.path.join(RootOutName, 'log.txt'), 'w')

    for root, dirs, RSFiles in os.walk(RootInputPath):

        # 判断是否进入最底层
        if len(dirs) == 0:
            # 根据输入输出路径建立生成新文件的路径
            RootInputPathList = RootInputPath.split(os.path.sep)
            RootList = root.split(os.path.sep)
            StartList = len(RootInputPathList)
            EndList = len(RootList)
            outname = RootOutName
            for i in range(StartList, EndList):
                if os.path.exists(os.path.join(outname, RootList[i])) == False:
                    os.makedirs(os.path.join(outname, RootList[i]))
                    outname = os.path.join(outname, RootList[i])
                else:
                    outname = os.path.join(outname, RootList[i])

            MeteData = glob.glob(os.path.join(root, '*MTL.txt'))[0]

            # for MeteData in MeteDatas:
            #     pass

            f = open(MeteData)
            data = f.readlines()
            data2 = ' '.join(data)

            shutil.copyfile(MeteData, os.path.join(outname, os.path.basename(MeteData)))

            if len(os.path.basename(MeteData)) < 10:

                RSbands = glob.glob(os.path.join(root, "B0[1-8].tiff"))
            else:
                RSbands = glob.glob(os.path.join(root, "*B[1-8].TIF"))
            print('影像' + root + '开始大气校正')
            #print(RSbands)
            for tifFile in RSbands:

                BandId = (os.path.basename(tifFile).split('.')[0])[-1]

                # 捕捉打开数据出错异常
                try:
                    IDataSet = gdal.Open(tifFile)
                except Exception as e:
                    print("文件%S打开失败" % tifFile)
                    LogFile.write('\n' + tifFile + '数据打开失败')

                if IDataSet == None:
                    LogFile.write('\n' + tifFile + '数据集读取为空')
                    continue
                else:
                    # 获取行列号
                    cols = IDataSet.RasterXSize
                    rows = IDataSet.RasterYSize
                    ImgBand = IDataSet.GetRasterBand(1)
                    ImgRasterData = ImgBand.ReadAsArray(0, 0, cols, rows)

                    if ImgRasterData is None:
                        LogFile.write('\n' + tifFile + '栅格数据为空')
                        continue
                    else:
                        # 设置输出文件路径
                        outFilename = os.path.join(outname, os.path.basename(tifFile))
                        outFilename=outFilename.replace(".TIF",".tiff")
                        # 如果文件存在就跳过，进行下一波段操作
                        if os.path.isfile(outFilename):
                            print("%s已经完成" % outFilename)
                            continue
                        else:
                            # #辐射校正
                            RaCalRaster = RadiometricCalibration(data2=data2,ImgRasterData=ImgRasterData,BandId=BandId)
                            #print(np.max(RaCalRaster))
                            # 大气校正
                            a, b, c = AtmosphericCorrection(data2=data2,BandId= BandId)
                            y = numpy.where(RaCalRaster != 0, a * RaCalRaster - b, 0)
                            atc = numpy.where(y != 0, (y / (1 + y * c)) * 10000, 0)

                            driver = IDataSet.GetDriver()
                            # 输出栅格数据集
                            outDataset = driver.Create(outFilename, cols, rows, 1, gdal.GDT_UInt16)

                            # 设置投影信息，与原数据一样
                            geoTransform = IDataSet.GetGeoTransform()
                            outDataset.SetGeoTransform(geoTransform)
                            proj = IDataSet.GetProjection()
                            outDataset.SetProjection(proj)

                            outband = outDataset.GetRasterBand(1)
                            outband.SetNoDataValue(0)
                            outband.WriteArray(atc, 0, 0)
                print('第%s波段矫正完成' % BandId)
            # print(root+'计算完成')
            print('\n')
    # 关闭日志文件
    LogFile.close()
    print(name," 波段合成")
    stackname=name+"_"+"B1_B7_Stack.tiff"
    layerstackpath=os.path.join(resultdir,name,stackname)
    Stackband(RootOutName,layerstackpath)
    print("layerstack finish！")

    print(name," 融合")
    pantif=glob.glob(os.path.join(RootInputPath,"*B8.TIF"))[0]
    multitif=layerstackpath
    Pansharpentif(multitifpath=multitif,pantifpath=pantif)



if __name__ == '__main__':
    '''tardir:原始压缩文件路径
    '''
    # ProcessBlock("E:\L8data\LC08_L2SP_120044_20200713_20200912_02_T1",r"E:\L8data\result")

    tardir=r"E:\L8data"
    tarfiles=glob.glob(os.path.join(tardir,"*.tar.gz"))
    #多进程 2核
    L8pool=Pool(2)
    L8pool.map(ProcessBlock,tarfiles)
    L8pool.close()
    L8pool.join()
    
    # for taritem in tarfiles:
    #     ProcessBlock(taritem)
    