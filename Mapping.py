############## Block 1 ##########################

import os
from osgeo import gdal, osr, ogr
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from rasterio import plot
import rasterio

os.environ['PROJ_LIB'] = r'C:\Users\Administrator\AppData\Roaming\Python\Python310\site-packages\pyproj\proj_dir\share\proj'

###### 1 a)
# define the file directory and sentinel images names
working_dir = r'E:\assignment2\A2_data'
os.chdir(working_dir)
sentinel_april_file = "sentinel2_2020_04_17_1045.tiff"
sentinel_may_file = "sentinel2_2020_05_05_1051.tiff"
# open one of the sentinel images and calculate width and height of the raster image
ds_sentinel_may = None
ds_sentinel_may = gdal.Open(sentinel_may_file)
print("sentinel image opened.")

width = ds_sentinel_may.RasterXSize
height = ds_sentinel_may.RasterYSize
print("x size: ", width, " y size: ", height)
print()

###### 1 b)
# use the getProjection() method to get the projection of the raster
# import osr and use SpatialReference to convert the format of projection to wkt
raster_projection = ds_sentinel_may.GetProjection()
spatial_ref = osr.SpatialReference(wkt=raster_projection)
print("EPSG Code:", spatial_ref)
print()

###### 1 c)
# access geotransform information
raster_geotransform = ds_sentinel_may.GetGeoTransform()
# top-left x coordinate
ulx = raster_geotransform[0]
# top-left y coordinate
uly = raster_geotransform[3]
# bottom-right x coordinate
lrx = ulx+raster_geotransform[1]*width
# bottom-left y coordinate
lry = uly+raster_geotransform[5]*height
extent = [ulx, lry, lrx, uly]

print("extent:", extent)
print()

###### 1 d)
# area = (w-e pixel resolution)*(n-s pixel resolution)
pixel_area = raster_geotransform[1]*abs(raster_geotransform[5])
# square meters converted to hectares
print("Pixel area in ha = ", pixel_area / 10000)
print()

###### 1 e)
# use RasterCount() method to calculate the number of bands
nr_bands = ds_sentinel_may.RasterCount
print("There are " + str(nr_bands) + " bands")
# use GetRasterBand() method to get no date value of every band
band = ds_sentinel_may.GetRasterBand(1)
no_data_value = band.GetNoDataValue()
print("No data value", no_data_value)
print()

################# Block 2 ###############################

####### 2 a)
# define the file directory and open the shapefile
aoi_datafile = r'E:\assignment2\A2_data'
os.chdir(aoi_datafile)
ds_aoi = ogr.Open("aoi.shp")
print("AOI shapefile opened")
# use GetLayer() and GetExtent() methods to get the aoi extent
layer = ds_aoi.GetLayer()
extent_AOI = layer.GetExtent()
print("Aoi extent = ", [extent_AOI[0], extent_AOI[1], extent_AOI[2], extent_AOI[3]])

####### 2 b)
# define the file path of aoi shapefile
aoi_path = r'E:\assignment2\A2_data\aoi.shp'
# use GetFeatureCount() to calculate the number of feature
featuresnum_aoi = layer.GetFeatureCount()
print("Feature count = ", featuresnum_aoi)
# use area() method from geopandas library to calculate the value of feature area
for feature in layer:
    feature = gpd.read_file(aoi_path)
    feature = feature.to_crs(3857)
    area = feature.area/10000
print("Area = ", area)

####### 2 c)
# define file paths of shapefile and raster images
input_file = r'E:/assignment2/A2_data/'
output_file = r'E:/assignment2/A2_data/Clip/'
input_shape = r'E:/assignment2/A2_data/aoi.shp'
prefix = '.tiff'

# If the output path does not exist then create a new
if not os.path.exists(output_file):
    os.mkdir(output_file)
# Read raster images one by one for batch cropping
file_all = os.listdir(input_file)
for file_i in file_all:
    if file_i.endswith(prefix):
        file_name = input_file+file_i
        data = gdal.Open(file_name)

        ds = gdal.Warp(output_file+os.path.splitext(file_i)[0]+'_Clip.tiff',
                       data,
                       cropToCutline=True,
                       format='GTiff',
                       cutlineDSName=input_shape,
                       cutlineWhere=None,
                       dstNodata=0)

####################### Block 3 ##########################

####### 3 a)
dataDirectory= r'E:/assignment2/A2_data/Clip/'
os.chdir(dataDirectory)

#april
#define a function:'calculateNDVI(nir_april,red_april)'
#calculate NDVI in April
#NDVI = (NIR - RED) / (NIR+RED)
dsapril=gdal.Open("sentinel2_2020_04_17_1045_Clip.tiff")
def calculateNDVI(nir_april,red_april):
    np.seterr(divide='ignore', invalid='ignore')
    return ((nir_april - red_april) / (nir_april + red_april))

red_april=dsapril.GetRasterBand(1).ReadAsArray().astype(np.float32)
green_april=dsapril.GetRasterBand(2).ReadAsArray().astype(np.float32)
blue_april=dsapril.GetRasterBand(3).ReadAsArray().astype(np.float32)
nir_april=dsapril.GetRasterBand(4).ReadAsArray().astype(np.float32)

#may
#define a function:'calculateNDVI(nir_may,red_may)'
#calculate NDVI in May
#NDVI = (NIR - RED) / (NIR+RED)
dsmay=gdal.Open("sentinel2_2020_05_05_1051_Clip.tiff")
def calculateNDVI(nir_may,red_may):
    np.seterr(divide='ignore', invalid='ignore')
    return ((nir_may - red_may) / (nir_may + red_may))

red_may=dsmay.GetRasterBand(1).ReadAsArray().astype(np.float32)
green_may=dsmay.GetRasterBand(2).ReadAsArray().astype(np.float32)
blue_may=dsmay.GetRasterBand(3).ReadAsArray().astype(np.float32)
nir_may=dsmay.GetRasterBand(4).ReadAsArray().astype(np.float32)

####### 3 b)
#generate NDVI datasets(TIFF) in April
drv_april = gdal.GetDriverByName ( "GTiff" )
dst_dsapril= drv_april.Create ( "output_april.tif", dsapril.RasterXSize,
                          dsapril.RasterYSize, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"] )
dst_dsapril.GetRasterBand(1).WriteArray (calculateNDVI(nir_april,red_april))
dst_dsapril = None


#generate NDVI datasets(TIFF) in May
drv_may = gdal.GetDriverByName ( "GTiff" )
dst_dsmay= drv_may.Create ( "output_may.tif", dsmay.RasterXSize,
                        dsmay.RasterYSize, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"] )
dst_dsmay.GetRasterBand(1).WriteArray (calculateNDVI(nir_may,red_may))
dst_dsmay = None

print()

########################## Block 4 #########################

####### 4 a)
#calculate dNDVI
#define a function:"dNDVI(before, after)"
#dNDVI = NDVI before the fire incident – NDVI after the fire incident
ds_before=gdal.Open("output_april.tif")
ds_after=gdal.Open("output_may.tif")
def dNDVI( before,after ):
    np.seterr(divide='ignore', invalid='ignore')
    return ( before-after )

before = ds_before.GetRasterBand(1).ReadAsArray().astype(np.float32)
after = ds_after.GetRasterBand(1).ReadAsArray().astype(np.float32)

drv_dNDVI = gdal.GetDriverByName ( "GTiff" )
dst_dsdNDVI= drv_dNDVI.Create ( "output_dNDVI.tif", ds_before.RasterXSize,
                                ds_before.RasterYSize, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"] )
dst_dsdNDVI.GetRasterBand(1).WriteArray(dNDVI(before, after))

####### 4 b)
#Reclassification
#dNDVI < 0.1------- new_pixel_value=0--- Unburnt pixels
#0.1<=dNDVI< 0.2 ---new_pixel_value=1--- low severity burns
#0.2<=dNDVI< 0.3 ---new_pixel_value=2--- moderate severity burns
#dNDVI>=0.3-------- new_pixel_value=3--- high severity burns
#use matplotlib method  to reclassify
new_pixel_value = dNDVI(before,after)
new_pixel_value[np.where( new_pixel_value >= 0.3 )] = 3
new_pixel_value[np.where((new_pixel_value >= 0.2) & (new_pixel_value < 0.3))] = 2
new_pixel_value[np.where((new_pixel_value >= 0.1) & (new_pixel_value < 0.2) )] = 1
new_pixel_value[np.where(new_pixel_value < 0.1)] = 0

####### 4 c)
#Save the reclassified dNDVI to a file.
driver = gdal.GetDriverByName ( "GTiff" )
ds_output = driver.Create("dndvi.tif",ds_before.RasterXSize,
                                       ds_before.RasterYSize, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])
ds_output.GetRasterBand(1).WriteArray (new_pixel_value)
#ds_output.SetGeoTransform(ds_sentinel_april_clipped.GetGeoTransform())
#ds_output.SetProjection(ds_sentinel_may_clipped.GetProjection())

####### 4 d)
# Visualize the pixel array of the dNDVI dataset in a suitable plot.
# add title, legend, colormap, and label

import matplotlib.colors as colors
plt.figure(); #ax = plt.subplots()

cmapCHM = colors.ListedColormap(['lightblue','yellow','red','black'])
import matplotlib.patches as mpatches
class1_box = mpatches.Patch(color='lightblue', label='Unburnt pixels')
class2_box = mpatches.Patch(color='yellow', label='Low severity burns')
class3_box = mpatches.Patch(color='red', label='Moderate severity burns')
class4_box = mpatches.Patch(color='black', label='High severity burns')

plt.imshow(new_pixel_value,cmap=cmapCHM)
plt.title('Forest Fire Severity')


plt.legend(handles=[class1_box,class2_box,class3_box,class4_box],
          handlelength=0.7,bbox_to_anchor=(0.7, -0.13),loc='upper left',
           prop={'size': 7.5},borderaxespad=0.)
plt.show()

####### 4 e)
#Estimate area of the forest that experienced moderate and high severity burns in hectare.
#calculate the total of pixels according to the condition(moderate and high severity burns)
sum=0
for i in new_pixel_value:
    for k in i:
        if k ==2:
            sum+=1
        elif k==3:
            sum += 1

#raster resolution is 10m
#the cell size in the raster is 10m
#calculte the area in square meter
area=sum*100

#Conversion square meter to hectare
print("Area in ha = ", area/10000 , "hectare")

############### Block 5 ##################
import rasterio
import numpy as np

#read the datasets
src_april=rasterio.open(r'E:/assignment2/A2_data/Clip/sentinel2_2020_04_17_1045_Clip.tiff')
src_may=rasterio.open(r'E:/assignment2/A2_data/Clip/sentinel2_2020_05_05_1051_Clip.tiff')

#get the band from the datasets
B1, B2, B3, B4 = src_april.read()
b1, b2, b3, b4 = src_may.read()

#get red and near infrared bands
#use numpy to calculate NDVI
np.seterr(divide='ignore', invalid='ignore')
ndvi_april = np.zeros(B1.shape)
ndvi_april = ((B4-B1)/(B1+B4))

np.seterr(divide='ignore', invalid='ignore')
ndvi_may = np.zeros(b1.shape)
ndvi_may = ((b4-b1)/(b4+b1))


#export the result
import matplotlib.pyplot as plt
from rasterio import plot

image_color_april=plt.imshow(ndvi_april)
plt.colorbar()
plt.title('NDVI Before the Fire')
plot.show(ndvi_april)

image_color_may=plt.imshow(ndvi_may)
plt.colorbar()
plt.title('NDVI After the Fire')
plot.show(ndvi_may)