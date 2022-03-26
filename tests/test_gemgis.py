"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke

GemGIS is a Python-based, open-source spatial data processing library.
It is capable of preprocessing spatial data such as vector data
raster data, data obtained from online services and many more data formats.
GemGIS wraps and extends the functionality of packages known to the geo-community
such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE) for more details.

"""

import numpy as np
from owslib import feature
from owslib.feature import wfs100
import pytest
import rasterio
import pandas as pd
from shapely import geometry
import geopandas as gpd
import gempy as gp
import gemgis as gg

gg.download_gemgis_data.download_tutorial_data(filename='test_gemgis.zip',
                                               dirpath='../docs/getting_started/tutorial/data/test_gemgis/')

__all__ = [geometry, feature, wfs100]
# Definition of GeoDataFrames
###########################################################

lines1 = geometry.linestring.LineString([[0.256327195431048, 264.86214748436396],
                                         [10.59346813871597, 276.73370778641777],
                                         [17.134940141888464, 289.089821570188],
                                         [19.150128045807676, 293.313485355882],
                                         [27.79511673965105, 310.571692592952]])

lines2 = geometry.linestring.LineString([[0.1881868620686138, 495.787213546976],
                                         [8.840672956663411, 504.1418419288791],
                                         [41.09257661030145, 546.423052386348],
                                         [71.72834900044251, 604.1435562935771],
                                         [87.60971085508348, 626.5358946968061]])

lines3 = geometry.linestring.LineString([[970.6766251230017, 833.0526164998309],
                                         [959.3724321757514, 800.0232029873156],
                                         [941.2919969974961, 754.801239252839],
                                         [925.3512303815758, 711.3618892206762],
                                         [908.0193245862978, 675.9805395602916]])

gdf_interfaces1_lines = gpd.GeoDataFrame(geometry=[lines1, lines2, lines3], crs='EPSG:4326')
gdf_interfaces1_lines['formation'] = ['Sand1', 'Ton', 'Ton']
gdf_interfaces1_lines['id'] = None

points1 = geometry.point.Point(19.150128045807676, 293.313485355882)

points2 = geometry.point.Point(61.93436666575576, 381.4593263680641)

points3 = geometry.point.Point(109.3578600758187, 480.9455679783049)

points4 = geometry.point.Point(157.812298994796, 615.9994296460927)

points5 = geometry.point.Point(191.3180280345144, 719.0939805375339)

gdf_interfaces1_points = gpd.GeoDataFrame(geometry=[points1, points2, points3, points4, points5], crs='EPSG:4326')
gdf_interfaces1_points['formation'] = 'Ton'
gdf_interfaces1_points['id'] = None

polygons1 = geometry.polygon.Polygon([[0.256327195431048, 264.862147484364],
                                      [10.59346813871597, 276.7337077864178],
                                      [17.13494014188846, 289.089821570188],
                                      [19.15012804580768, 293.313485355882],
                                      [27.79511673965105, 310.571692592952],
                                      [34.41734765644295, 324.1391900810135],
                                      [40.7165429187572, 338.5142767052691],
                                      [49.27698776241503, 352.566327675047],
                                      [55.33390628387104, 364.1148523226231],
                                      [60.98703023722999, 376.3094482791545],
                                      [61.93436666575576, 381.4593263680641],
                                      [74.31225098443318, 404.8981037004269],
                                      [89.49492674488292, 440.4320256929689],
                                      [100.8011746516008, 465.6288067422258],
                                      [109.3578600758187, 480.9455679783049],
                                      [122.121527847126, 511.4998696780527],
                                      [134.7199183717545, 543.1573638168628],
                                      [146.0261662784724, 575.1378936101505],
                                      [154.748128949369, 602.2728885862734],
                                      [157.812298994796, 615.9994296460927],
                                      [170.2538403642964, 655.5737715750863],
                                      [179.944909998626, 686.9082300594188],
                                      [191.3180280345144, 719.0939805375339],
                                      [200.6191918851958, 750.2232183370392],
                                      [210.9563328284808, 774.7739280773408],
                                      [224.8237570742327, 798.4767847239435],
                                      [240.6756130404249, 822.9062405945112],
                                      [255.2122174919194, 846.1648077169024],
                                      [264.903287126249, 861.9935547863074],
                                      [264.903287126249, 861.9935547863074],
                                      [272.7627232387529, 874.2512796291529],
                                      [291.3922107934166, 899.788726360193],
                                      [308.8361361352099, 915.9405090840756],
                                      [323.2790531755591, 927.3449733382452],
                                      [341.462737237453, 937.9069335885562],
                                      [355.0302347255144, 947.2749675684081],
                                      [366.0632917955071, 952.6031383066481],
                                      [384.7495149374586, 957.2890728572154],
                                      [394.4405845717882, 955.9969302393048],
                                      [407.8165849065408, 950.0257745343622],
                                      [423.1907578202995, 941.4603257878103],
                                      [441.2807544710481, 929.1849709176595],
                                      [456.2710238255182, 917.5509910035582],
                                      [474.5534268822465, 900.757833323626],
                                      [490.0591382971738, 879.4374801281008],
                                      [503.949671439713, 861.3474834773522],
                                      [515.0349178336396, 847.9621691518354],
                                      [530.7616307613583, 821.9371336310784],
                                      [545.6212708673304, 798.3555308542095],
                                      [559.8348396643472, 775.7430350407737],
                                      [575.3452301051327, 753.6306550861667],
                                      [592.1384051121127, 724.0573303243491],
                                      [604.0907243277858, 703.706084092257],
                                      [617.9812574703249, 678.1862673885223],
                                      [629.6105410315205, 658.804128119863],
                                      [638.748378903369, 641.7730673689532],
                                      [651.2539298815232, 612.9330651840361],
                                      [660.6219638613752, 594.1969972243322],
                                      [669.3439265322719, 575.4609292646282],
                                      [674.8314717153735, 563.4212086914579],
                                      [685.4957092561546, 542.1882568534299],
                                      [698.7401710897384, 515.053261877307],
                                      [710.6924903054116, 486.9491599377511],
                                      [720.0605242852636, 462.3984501974494],
                                      [727.4903443382495, 442.3702396198349],
                                      [734.1108384779523, 426.3054560058411],
                                      [745.9033766434759, 391.3306062123655],
                                      [759.1478384770597, 360.3191833825108],
                                      [767.5467654934788, 341.5831154228068],
                                      [771.7403495533283, 327.3346871500575],
                                      [779.1760490546743, 305.080086466832],
                                      [785.3137264897498, 279.2372341086196],
                                      [795.967569012817, 245.8899919458189],
                                      [805.9880083763195, 216.8913527944325],
                                      [812.4487214658725, 198.4783204892062],
                                      [818.648370208934, 185.5796796743258],
                                      [826.9853259173669, 173.9276107489045],
                                      [841.3291714050511, 164.9607694960376],
                                      [857.0276417837888, 173.6045750944268],
                                      [867.1028091279114, 188.1570434466119],
                                      [874.1485314711044, 203.323855306371],
                                      [886.7469219957329, 240.1499199168235],
                                      [893.8537063942413, 268.2540218563794],
                                      [900.3144194837944, 285.6979471981728],
                                      [903.701374694373, 295.8908491281679],
                                      [909.6824534636465, 319.6166909183264],
                                      [916.4662022076773, 337.706687569075],
                                      [920.1965028370037, 355.1702158907466],
                                      [929.3876283867835, 383.2547148504243],
                                      [936.1713771308141, 409.0975672086366],
                                      [945.539411110666, 440.7550613474467],
                                      [950.0939225955215, 453.1100392376158],
                                      [955.5535163994733, 474.3507694131226],
                                      [962.6603007979818, 497.6093365355138],
                                      [966.0735779836949, 510.3275149823656],
                                      [972.0889039630482, 529.79176361285],
                                      [971.8062477653798, -0.006899458905866851],
                                      [0.2361374670297192, 0.1546183683329865],
                                      [0.256327195431048, 264.862147484364]])

polygons2 = geometry.polygon.Polygon([[0.256327195431048, 264.862147484364],
                                      [0.1881868620686138, 495.787213546976],
                                      [8.840672956663411, 504.1418419288791],
                                      [41.09257661030145, 546.423052386348],
                                      [71.72834900044251, 604.1435562935771],
                                      [87.60971085508348, 626.5358946968061],
                                      [116.3598841035946, 686.6205264296495],
                                      [140.8016980977082, 751.5687640683381],
                                      [181.6130863080805, 848.1383536684762],
                                      [223.2773388108611, 950.5412472888196],
                                      [238.4673614961476, 988.0127920572999],
                                      [253.1747585693791, 1022.707432912829],
                                      [265.6023564722704, 1045.513138554322],
                                      [278.5237826513765, 1068.771705676713],
                                      [511.6747662706227, 1068.852464590333],
                                      [526.3753184316984, 1045.388234108946],
                                      [560.1099247138658, 990.6172670215257],
                                      [608.8509591448512, 912.3962634589865],
                                      [636.0233035161142, 859.787826958076],
                                      [677.9243082421169, 782.4971293357705],
                                      [700.6304344116447, 739.6185634923892],
                                      [723.285910634351, 693.8358155691311],
                                      [740.0407842579184, 658.5366142184982],
                                      [757.3071124285267, 626.8243574896944],
                                      [772.6673853601612, 603.6205529572973],
                                      [790.2973687137878, 583.524646115289],
                                      [810.139521279569, 572.9321657819203],
                                      [831.5351890703644, 569.0914089904873],
                                      [852.4571920161416, 577.4546649446074],
                                      [869.6801729001976, 594.8650467133475],
                                      [882.8225435370409, 621.7105496080459],
                                      [898.5466471498011, 648.474213176897],
                                      [908.0193245862978, 675.9805395602916],
                                      [925.3512303815758, 711.3618892206762],
                                      [941.2919969974961, 754.801239252839],
                                      [959.3724321757514, 800.0232029873156],
                                      [971.7381074320145, 832.9758676364308],
                                      [972.0889039630482, 529.79176361285],
                                      [966.0735779836949, 510.3275149823656],
                                      [962.6603007979818, 497.6093365355138],
                                      [955.5535163994733, 474.3507694131226],
                                      [950.0939225955215, 453.1100392376158],
                                      [945.539411110666, 440.7550613474467],
                                      [936.1713771308141, 409.0975672086366],
                                      [929.3876283867835, 383.2547148504243],
                                      [920.1965028370037, 355.1702158907466],
                                      [916.4662022076773, 337.706687569075],
                                      [909.6824534636465, 319.6166909183264],
                                      [903.701374694373, 295.8908491281679],
                                      [893.8537063942413, 268.2540218563794],
                                      [886.7469219957329, 240.1499199168235],
                                      [874.1485314711044, 203.323855306371],
                                      [867.1028091279114, 188.1570434466119],
                                      [857.0276417837888, 173.6045750944268],
                                      [841.3291714050511, 164.9607694960376],
                                      [826.9853259173669, 173.9276107489045],
                                      [818.648370208934, 185.5796796743258],
                                      [812.4487214658725, 198.4783204892062],
                                      [805.9880083763195, 216.8913527944325],
                                      [795.967569012817, 245.8899919458189],
                                      [785.3137264897498, 279.2372341086196],
                                      [779.1760490546743, 305.080086466832],
                                      [771.7403495533283, 327.3346871500575],
                                      [767.5467654934788, 341.5831154228068],
                                      [759.1478384770597, 360.3191833825108],
                                      [745.9033766434759, 391.3306062123655],
                                      [734.1108384779523, 426.3054560058411],
                                      [727.4903443382495, 442.3702396198349],
                                      [720.0605242852636, 462.3984501974494],
                                      [710.6924903054116, 486.9491599377511],
                                      [698.7401710897384, 515.053261877307],
                                      [685.4957092561546, 542.1882568534299],
                                      [674.8314717153735, 563.4212086914579],
                                      [669.3439265322719, 575.4609292646282],
                                      [660.6219638613752, 594.1969972243322],
                                      [651.2539298815232, 612.9330651840361],
                                      [638.748378903369, 641.7730673689532],
                                      [629.6105410315205, 658.804128119863],
                                      [617.9812574703249, 678.1862673885223],
                                      [604.0907243277858, 703.706084092257],
                                      [575.3452301051327, 753.6306550861667],
                                      [559.8348396643472, 775.7430350407737],
                                      [545.6212708673304, 798.3555308542095],
                                      [530.7616307613583, 821.9371336310784],
                                      [515.0349178336396, 847.9621691518354],
                                      [503.949671439713, 861.3474834773522],
                                      [490.0591382971738, 879.4374801281008],
                                      [474.5534268822465, 900.757833323626],
                                      [456.2710238255182, 917.5509910035582],
                                      [441.2807544710481, 929.1849709176595],
                                      [423.1907578202995, 941.4603257878103],
                                      [407.8165849065408, 950.0257745343622],
                                      [394.4405845717882, 955.9969302393048],
                                      [384.7495149374586, 957.2890728572154],
                                      [366.0632917955071, 952.6031383066481],
                                      [355.0302347255144, 947.2749675684081],
                                      [341.462737237453, 937.9069335885562],
                                      [323.2790531755591, 927.3449733382452],
                                      [308.8361361352099, 915.9405090840756],
                                      [291.3922107934166, 899.788726360193],
                                      [291.3922107934166, 899.788726360193],
                                      [272.7627232387529, 874.2512796291529],
                                      [264.903287126249, 861.9935547863074],
                                      [255.2122174919194, 846.1648077169024],
                                      [240.6756130404249, 822.9062405945112],
                                      [224.8237570742327, 798.4767847239435],
                                      [210.9563328284808, 774.7739280773408],
                                      [200.6191918851958, 750.2232183370392],
                                      [191.3180280345144, 719.0939805375339],
                                      [179.944909998626, 686.9082300594188],
                                      [170.2538403642964, 655.5737715750863],
                                      [157.812298994796, 615.9994296460927],
                                      [154.748128949369, 602.2728885862734],
                                      [146.0261662784724, 575.1378936101505],
                                      [134.7199183717545, 543.1573638168628],
                                      [122.121527847126, 511.4998696780527],
                                      [109.3578600758187, 480.9455679783049],
                                      [100.8011746516008, 465.6288067422258],
                                      [89.49492674488292, 440.4320256929689],
                                      [74.31225098443318, 404.8981037004269],
                                      [61.93436666575576, 381.4593263680641],
                                      [60.98703023722999, 376.3094482791545],
                                      [55.33390628387104, 364.1148523226231],
                                      [49.27698776241503, 352.566327675047],
                                      [40.7165429187572, 338.5142767052691],
                                      [34.41734765644295, 324.1391900810135],
                                      [27.79511673965105, 310.571692592952],
                                      [19.15012804580768, 293.313485355882],
                                      [17.13494014188846, 289.089821570188],
                                      [10.59346813871597, 276.7337077864178],
                                      [0.256327195431048, 264.862147484364]])

polygons3 = geometry.polygon.Polygon([[0.1881868620686138, 495.787213546976],
                                      [0.2489663569543978, 1068.759507715801],
                                      [278.5237826513765, 1068.771705676713],
                                      [253.1747585693791, 1022.707432912829],
                                      [238.4673614961476, 988.0127920572999],
                                      [223.2773388108611, 950.5412472888196],
                                      [181.6130863080805, 848.1383536684762],
                                      [140.8016980977082, 751.5687640683381],
                                      [116.3598841035946, 686.6205264296495],
                                      [87.60971085508348, 626.5358946968061],
                                      [71.72834900044251, 604.1435562935771],
                                      [41.09257661030145, 546.423052386348],
                                      [8.840672956663411, 504.1418419288791],
                                      [0.1881868620686138, 495.787213546976]])

polygons4 = geometry.polygon.Polygon([[511.6747662706227, 1068.852464590333],
                                      [971.6979382849257, 1068.79988717261],
                                      [971.7381074320145, 832.9758676364308],
                                      [959.3724321757514, 800.0232029873156],
                                      [941.2919969974961, 754.801239252839],
                                      [925.3512303815758, 711.3618892206762],
                                      [908.0193245862978, 675.9805395602916],
                                      [898.5466471498011, 648.474213176897],
                                      [882.8225435370409, 621.7105496080459],
                                      [869.6801729001976, 594.8650467133475],
                                      [852.4571920161416, 577.4546649446074],
                                      [831.5351890703644, 569.0914089904873],
                                      [810.139521279569, 572.9321657819203],
                                      [790.2973687137878, 583.524646115289],
                                      [772.6673853601612, 603.6205529572973],
                                      [757.3071124285267, 626.8243574896944],
                                      [740.0407842579184, 658.5366142184982],
                                      [723.285910634351, 693.8358155691311],
                                      [700.6304344116447, 739.6185634923892],
                                      [677.9243082421169, 782.4971293357705],
                                      [636.0233035161142, 859.787826958076],
                                      [608.8509591448512, 912.3962634589865],
                                      [560.1099247138658, 990.6172670215257],
                                      [526.3753184316984, 1045.388234108946],
                                      [511.6747662706227, 1068.852464590333]])

gdf_geolmap1 = gpd.GeoDataFrame(geometry=[polygons1, polygons2, polygons3, polygons4], crs='EPSG:4326')
gdf_geolmap1['id'] = None
gdf_geolmap1['formation'] = ['Sand1', 'Ton', 'Sand2', 'Sand2']

points1 = geometry.point.Point(695.4667461080886, 3.226225077137428)

points2 = geometry.point.Point(669.2840030245482, 1060.822026058724)

gdf_customsection1 = gpd.GeoDataFrame(geometry=[points1, points2], crs='EPSG:4326')
gdf_customsection1['id'] = None

lines1 = geometry.linestring.LineString([[62.76372633685696, 44.51145167379445],
                                         [641.6436191608124, 1036.876982229147]])

lines2 = geometry.linestring.LineString([[863.8921494414382, 52.26430738125828],
                                         [168.7194210055274, 1021.371270814219]])

gdf_customsection1_lines = gpd.GeoDataFrame(geometry=[lines1, lines2], crs='EPSG:4326')
gdf_customsection1_lines['id'] = None
gdf_customsection1_lines['section'] = ['Section1', 'Section2']

polygon = geometry.polygon.Polygon([[0, 1069],
                                    [972, 1069],
                                    [972, 0],
                                    [0, 0],
                                    [0, 1069]])

gdf_extent = gpd.GeoDataFrame(geometry=[polygon], crs='EPSG:4326')
gdf_extent['id'] = None

points1 = geometry.point.Point(0, 1069)
points2 = geometry.point.Point(972, 1069)
points3 = geometry.point.Point(972, 0)
points4 = geometry.point.Point(0, 0)

gdf_extent_points = gpd.GeoDataFrame(geometry=[points1, points2, points3, points4], crs='EPSG:4326')
gdf_extent_points['id'] = None

lines1 = geometry.linestring.LineString([[1642.839442255191, 2582.579025497104],
                                         [2829.347688067769, 2205.936550313893]])

lines2 = geometry.linestring.LineString([[1705.331691500837, 1759.201417193002],
                                         [2875.794630074668, 1406.767651177023]])

lines3 = geometry.linestring.LineString([[1017.766042705972, 1722.234095004106],
                                         [2979.938220033745, 1137.002625717191]])

lines4 = geometry.linestring.LineString([[99.95647777150653, 1763.423866466358],
                                         [765.8367281794893, 1620.70508102698],
                                         [2963.048422940328, 955.4373069629546]])

lines5 = geometry.linestring.LineString([[1078.14706731494, 1313.501005343405],
                                         [2963.048422940328, 752.759741841945]])

lines6 = geometry.linestring.LineString([[1387.23035412448, 1097.311602547661],
                                         [2968.115362068354, 619.3303448039474]])

lines7 = geometry.linestring.LineString([[1152.462174525977, 806.8070925408806],
                                         [2966.426382359013, 249.4437884581041]])

lines8 = geometry.linestring.LineString([[1103.481762955068, 3493.9738101036],
                                         [2974.871280905723, 2872.429277065837]])

gdf_lines5 = gpd.GeoDataFrame(geometry=[lines1, lines2, lines3, lines4, lines5, lines6, lines7, lines8],
                              crs='EPSG:4326')
gdf_lines5['id'] = None
gdf_lines5['Z'] = [0, 150, 200, 250, 200, 150, 0, 250]
gdf_lines5['formation'] = 'Coal1'

points1 = geometry.point.Point(96.47104121438838, 451.5636209742439)

points2 = geometry.point.Point(172.7610088740548, 661.8765047927839)

points3 = geometry.point.Point(383.0738926925949, 957.7578658512201)

points4 = geometry.point.Point(592.3558310022205, 722.7022898187342)

points5 = geometry.point.Point(766.5856220087561, 348.4690700828027)

points6 = geometry.point.Point(843.906535177337, 167.0226605138662)

points7 = geometry.point.Point(941.8463585242062, 428.8828197781268)

points8 = geometry.point.Point(22.14220797837953, 299.5527565162123)

gdf_orientations1 = gpd.GeoDataFrame(
    geometry=[points1, points2, points3, points4, points5, points6, points7, points8, ], crs='EPSG:4326')
gdf_orientations1['id'] = None
gdf_orientations1['formation'] = 'Ton'
gdf_orientations1['dip'] = 30.5
gdf_orientations1['azimuth'] = 180


# Testing the GemPyData Class
###########################################################
def test_gem_py_data_empty():
    from gemgis import GemPyData
    data = GemPyData()
    assert data.model_name is None
    assert data.crs is None
    assert data.interfaces is None
    assert data.orientations is None
    assert data.extent is None
    assert data.section_dict is None
    assert data.resolution is None
    assert data.dem is None
    assert data.stack is None
    assert data.surface_colors is None
    assert data.is_fault is None
    assert data.geolmap is None
    assert data.faults is None
    assert data.tectonics is None
    assert data.raw_i is None
    assert data.raw_o is None
    assert data.raw_dem is None
    assert data.wms is None
    assert data.slope is None
    assert data.hillshades is None
    assert data.aspect is None
    assert data.basemap is None
    assert data.customsections is None
    assert data.contours is None


@pytest.mark.parametrize("interface_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1']]),
                                          columns=['X', 'Y', 'Z', 'formation'])
                         ])
@pytest.mark.parametrize("orientation_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1', 45, 90, 1]]),
                                          columns=['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'])
                         ])
@pytest.mark.parametrize("gdf_geolmap1", [gdf_geolmap1])
@pytest.mark.parametrize("gdf_interfaces1_lines", [gdf_interfaces1_lines])
def test_gem_py_data(interface_df, orientation_df, gdf_geolmap1, gdf_interfaces1_lines):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1',
                     crs='EPSG:4326',
                     interfaces=interface_df,
                     orientations=orientation_df,
                     extent=[0, 100, 0, 100, 0, 100],
                     resolution=[50, 50, 50],
                     section_dict={'SectionA': ([0, 10], [0, 0], [100, 80])},
                     stack={'Layer1': 'Layer1',
                            'Layer2': ('Layer2', 'Layer3')},
                     dem='path/to/dem.tif',
                     surface_colors={'Layer1': '#FFFFFF',
                                     'Layer2': '#000000',
                                     'Layer3': '#111111'},
                     geolmap=gdf_geolmap1,
                     faults=gdf_interfaces1_lines,
                     is_fault=['Fault1', 'Fault2']
                     )
    assert isinstance(data.model_name, str)
    assert data.model_name == 'Model1'
    assert isinstance(data.crs, str)
    assert data.crs == 'EPSG:4326'
    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interface_df.columns)
    assert isinstance(data.orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(orientation_df.columns)
    assert isinstance(data.extent, list)
    assert all(isinstance(n, (int, float)) for n in data.extent)
    assert data.extent == [0, 100, 0, 100, 0, 100]
    assert all(isinstance(n, (int, float)) for n in data.extent)
    assert isinstance(data.resolution, list)
    assert all(isinstance(n, int) for n in data.resolution)
    assert data.resolution == [50, 50, 50]
    assert all(isinstance(n, (int, float)) for n in data.resolution)
    assert isinstance(data.section_dict, dict)
    assert all(isinstance(n, tuple) for n in [data.section_dict[key] for key in data.section_dict])
    assert data.section_dict == {'SectionA': ([0, 10], [0, 0], [100, 80])}
    assert isinstance(data.stack, dict)
    assert all(isinstance(n, (str, tuple)) for n in [data.stack[key] for key in data.stack])
    assert data.stack == {'Layer1': 'Layer1', 'Layer2': ('Layer2', 'Layer3')}
    assert isinstance(data.dem, str)
    assert data.dem == 'path/to/dem.tif'
    assert isinstance(data.surface_colors, dict)
    assert all(isinstance(n, str) for n in [data.surface_colors[key] for key in data.surface_colors])
    assert data.surface_colors == {'Layer1': '#FFFFFF', 'Layer2': '#000000', 'Layer3': '#111111'}
    assert isinstance(data.geolmap, gpd.geodataframe.GeoDataFrame)
    assert 'geometry' in data.geolmap
    assert isinstance(data.faults, gpd.geodataframe.GeoDataFrame)
    assert 'geometry' in data.faults
    assert isinstance(data.is_fault, list)
    assert all(isinstance(n, str) for n in data.is_fault)


@pytest.mark.parametrize("interface_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1']]),
                                          columns=['X', 'Y', 'Z', 'formation'])
                         ])
@pytest.mark.parametrize("orientation_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1', 45, 90, 1]]),
                                          columns=['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'])
                         ])
@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
def test_gem_py_data_errors(interface_df, orientation_df, gdf_interfaces1_points):
    from gemgis import GemPyData
    with pytest.raises(TypeError):
        GemPyData(model_name=['Model1'])
    with pytest.raises(TypeError):
        GemPyData(crs=['EPSG:4326'])
    with pytest.raises(TypeError):
        GemPyData(interfaces=[interface_df])
    with pytest.raises(TypeError):
        GemPyData(interfaces=[orientation_df])
    with pytest.raises(TypeError):
        GemPyData(extent=(0, 100, 0, 100, 0))
    with pytest.raises(TypeError):
        GemPyData(extent=(0, 100, 0, 100, 0, 100))
    with pytest.raises(ValueError):
        GemPyData(extent=[0, 100, 0, 100])
    with pytest.raises(TypeError):
        GemPyData(extent=[0, 100, 0, 100, 0, '100'])
    with pytest.raises(TypeError):
        GemPyData(resolution=(50, 50, 50))
    with pytest.raises(ValueError):
        GemPyData(resolution=[50, 50, 50, 50])
    with pytest.raises(TypeError):
        GemPyData(resolution=[0, 100, 100.0])
    with pytest.raises(TypeError):
        GemPyData(section_dict=[[0, 100], [0, 100], [0, 100]])
    with pytest.raises(TypeError):
        GemPyData(stack=[[0, 100], [0, 100], [0, 100]])
    with pytest.raises(TypeError):
        GemPyData(dem=['path/to/dem.tif'])
    with pytest.raises(TypeError):
        GemPyData(surface_colors=['#FFFFFF', '#000000', '#111111'])
    with pytest.raises(TypeError):
        GemPyData(geolmap=['#FFFFFF', '#000000', '#111111'])
    with pytest.raises(TypeError):
        GemPyData(geolmap=gdf_interfaces1_points)
    with pytest.raises(TypeError):
        GemPyData(faults=gdf_interfaces1_points)
    with pytest.raises(TypeError):
        GemPyData(is_fault=np.array[['Fault1', 'Fault2']])
    with pytest.raises(TypeError):
        GemPyData(basemap=['Test'])
    with pytest.raises(TypeError):
        GemPyData(faults=['Test'])
    with pytest.raises(TypeError):
        GemPyData(faults=gdf_interfaces1_points)
    with pytest.raises(TypeError):
        GemPyData(is_fault=['Fault1', ['Fault2']])
    with pytest.raises(TypeError):
        GemPyData(customsections=['Test'])
    with pytest.raises(TypeError):
        GemPyData(contours=['Test'])


# Testing data.to_section_dict
###########################################################

@pytest.mark.parametrize("gdf_customsection1", [gdf_customsection1])
def test_to_section_dict_points_data(gdf_customsection1):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_customsection1['section_name'] = 'SectionA'
    data.to_section_dict(gdf_customsection1, 'section_name', [100, 80])

    assert isinstance(gdf_customsection1, gpd.geodataframe.GeoDataFrame)
    assert isinstance('section', str)
    assert isinstance([100, 80], list)
    assert isinstance(data.section_dict, dict)
    assert data.section_dict['SectionA'] == (
        [695.4667461080886, 3.226225077137428], [669.2840030245482, 1060.822026058724],
        [100, 80])
    assert len(data.section_dict) == 1


@pytest.mark.parametrize("gdf_customsection1_lines", [gdf_customsection1_lines])
def test_to_section_dict_lines_data(gdf_customsection1_lines):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.to_section_dict(gdf_customsection1_lines, 'section', [100, 80])

    assert isinstance(gdf_customsection1_lines, gpd.geodataframe.GeoDataFrame)
    assert isinstance('section', str)
    assert isinstance([100, 80], list)
    assert isinstance(data.section_dict, dict)
    assert data.section_dict['Section1'] == (
        [62.76372633685696, 44.51145167379445], [641.6436191608124, 1036.876982229147],
        [100, 80])
    assert data.section_dict['Section2'] == (
        [863.8921494414382, 52.26430738125828], [168.7194210055274, 1021.371270814219],
        [100, 80])
    assert len(data.section_dict) == 2


@pytest.mark.parametrize("gdf_customsection1_lines", [gdf_customsection1_lines])
def test_to_section_dict_error_data(gdf_customsection1_lines):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    with pytest.raises(TypeError):
        data.to_section_dict([gdf_customsection1_lines], 'section', [100, 80])
    with pytest.raises(TypeError):
        data.to_section_dict(gdf_customsection1_lines, ['section'], [100, 80])
    with pytest.raises(TypeError):
        data.to_section_dict(gdf_customsection1_lines, 'section', (100, 80))
    with pytest.raises(ValueError):
        data.to_section_dict(gdf_customsection1_lines, 'section', [100, 80, 50])
    with pytest.raises(ValueError):
        data.to_section_dict(gdf_customsection1_lines, 'section_1', [100, 80, 50])


# Testing data.to_gempy_df
###########################################################

@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_gemgis/raster1.tif')
                         ])
def test_to_gempy_df_points_data(gdf_interfaces1_points, dem):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.to_gempy_df(gdf_interfaces1_points, cat='interfaces', dem=dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf_interfaces1_points
    assert all(gdf_interfaces1_points.geom_type == 'Point')
    assert not {'X', 'Y', 'Z'}.issubset(gdf_interfaces1_points.columns)

    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.3578600758187,
                                                     157.812298994796,
                                                     191.3180280345144]
    assert data.interfaces['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                                     615.9994296460927,
                                                     719.0939805375339]
    assert data.interfaces['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625,
                                                     525.6910400390625,
                                                     597.6325073242188]


    with pytest.raises(TypeError):
        data.to_gempy_df([gdf_interfaces1_points], cat='interfaces', dem=dem)

    with pytest.raises(TypeError):
        data.to_gempy_df(gdf_interfaces1_points, cat=['interfaces'], dem=dem)

    with pytest.raises(TypeError):
        data.to_gempy_df(gdf_interfaces1_points, cat='interfaces', dem=[dem])

    with pytest.raises(ValueError):
        data.to_gempy_df(gdf_interfaces1_points.drop('formation', axis=1), cat='interfaces', dem=dem)


@pytest.mark.parametrize("gdf_interfaces1_lines", [gdf_interfaces1_lines])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_gemgis/raster1.tif')
                         ])
def test_to_gempy_df_lines_data(gdf_interfaces1_lines, dem):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.to_gempy_df(gdf_interfaces1_lines, cat='interfaces', dem=dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf_interfaces1_points
    assert all(gdf_interfaces1_lines.geom_type == 'LineString')
    assert not {'X', 'Y', 'Z'}.issubset(gdf_interfaces1_points.columns)
    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                                     19.150128045807676,
                                                     27.79511673965105]
    assert data.interfaces['Y'].head().to_list() == [264.86214748436396, 276.73370778641777, 289.089821570188,
                                                     293.313485355882,
                                                     310.571692592952]
    assert data.interfaces['Z'].head().to_list() == [353.9727783203125, 359.03631591796875, 364.28497314453125,
                                                     364.994873046875,
                                                     372.81036376953125]


@pytest.mark.parametrize("gdf_interfaces1_lines", [gdf_interfaces1_lines])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_gemgis/raster1.tif')
                         ])
def test_to_gempy_df_lines_xyz_data(gdf_interfaces1_lines, dem):
    from gemgis.vector import extract_xyz
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_xyz = extract_xyz(gdf_interfaces1_lines, dem)
    data.to_gempy_df(gdf_xyz, cat='interfaces')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf_interfaces1_lines
    assert all(gdf_interfaces1_lines.geom_type == 'LineString')
    assert {'X', 'Y', 'Z', 'formation'}.issubset(gdf_xyz.columns)

    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                                     19.150128045807676,
                                                     27.79511673965105]
    assert data.interfaces['Y'].head().to_list() == [264.86214748436396, 276.73370778641777, 289.089821570188,
                                                     293.313485355882,
                                                     310.571692592952]
    assert data.interfaces['Z'].head().to_list() == [353.9727783203125, 359.03631591796875, 364.28497314453125,
                                                     364.994873046875,
                                                     372.81036376953125]


@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_gemgis/raster1.tif')
                         ])
def test_to_gempy_df_points_xyz_data(gdf_interfaces1_points, dem):
    from gemgis.vector import extract_xyz
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_xyz = extract_xyz(gdf_interfaces1_points, dem)
    data.to_gempy_df(gdf_xyz, cat='interfaces')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf_interfaces1_points
    assert all(gdf_interfaces1_points.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'formation'}.issubset(gdf_xyz.columns)

    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.3578600758187,
                                                     157.812298994796,
                                                     191.3180280345144]
    assert data.interfaces['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                                     615.9994296460927,
                                                     719.0939805375339]
    assert data.interfaces['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625,
                                                     525.6910400390625,
                                                     597.6325073242188]


# Testing data.set_extent
###########################################################

def test_set_extent_data():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(0, 100, 0, 100)

    assert isinstance(data.extent, list)
    assert len(data.extent) == 4
    assert data.extent == [0, 100, 0, 100]


def test_set_extent_z_data():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(0, 100, 0, 100, 0, 100)

    assert isinstance(data.extent, list)
    assert len(data.extent) == 6
    assert data.extent == [0, 100, 0, 100, 0, 100]


@pytest.mark.parametrize("gdf_extent", [gdf_extent])
def test_set_extent_z_data(gdf_extent):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(gdf=gdf_extent)

    assert isinstance(gdf_extent, gpd.geodataframe.GeoDataFrame)
    assert all(gdf_extent.geom_type == 'Polygon')
    assert 'geometry' in gdf_extent

    assert isinstance(data.extent, list)
    assert len(data.extent) == 4
    assert data.extent == [-0.0, 972.0, -0.0, 1069.0]


@pytest.mark.parametrize("gdf_extent_points", [gdf_extent_points])
def test_set_extent_z_data(gdf_extent_points):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(gdf=gdf_extent_points)

    assert isinstance(gdf_extent_points, gpd.geodataframe.GeoDataFrame)
    assert all(gdf_extent_points.geom_type == 'Point')
    assert 'geometry' in gdf_extent_points

    assert isinstance(data.extent, list)
    assert len(data.extent) == 6
    assert data.extent == [-0.0, 972.0, -0.0, 1069.0, 0, 0]


@pytest.mark.parametrize("gdf_extent_points", [gdf_extent_points])
def test_set_extent_error_data(gdf_extent_points):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.set_extent(gdf=[gdf_extent_points])
    with pytest.raises(TypeError):
        data.set_extent(0, 1.1, 2, 3, 4, [5])


# Testing set_resolution
###########################################################

def test_set_resolution_go():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_resolution(50, 50, 50)

    assert isinstance(data.resolution, list)
    assert all(isinstance(n, int) for n in data.resolution)
    assert len(data.resolution) == 3
    assert data.resolution == [50, 50, 50]


def test_set_resolution_error():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.set_resolution(50.0, 50, 50)

    with pytest.raises(TypeError):
        data.set_resolution(50, 50.0, 50)

    with pytest.raises(TypeError):
        data.set_resolution(50, 50, 50.0)

    with pytest.raises(TypeError):
        data.set_resolution(50, 50, 50, 50)


# Testing data.to_surface_color_dict
###########################################################

def test_create_surface_color_dict():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    data.to_surface_color_dict('../docs/getting_started/tutorial/data/test_gemgis/style1.qml')

    assert isinstance(data.surface_colors, dict)
    assert data.surface_colors == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.to_surface_color_dict(['../docs/getting_started/tutorial/data/test_gemgis/style1.qml'])

    with pytest.raises(TypeError):
        data.to_surface_color_dict(['../docs/getting_started/tutorial/data/test_gemgis/style1.qml'],
                                   basement=['Basement'])



# Testing calculate_number_of_isopoints
###########################################################
@pytest.mark.parametrize("gdf_lines5", [gdf_lines5])
def test_calculate_number_of_isopoints(gdf_lines5):
    from gemgis.utils import calculate_number_of_isopoints

    number = calculate_number_of_isopoints(gdf_lines5, 50, zcol='Z')
    assert number == 2

# Testing calculate_number_of_isopoints
###########################################################
# @pytest.mark.parametrize("gdf_lines5", [gdf_lines5])
# def test_calculate_lines(gdf_lines5):
#     from gemgis.utils import calculate_lines
#
#     gdf_lines5['X'] = 500
#     gdf_lines5['Y'] = 100
#
#     gdf = gdf_lines5[gdf_lines5.is_valid]
#
#     lines = calculate_lines(gdf, 50, xcol='X', zcol='Z')
#
#     assert isinstance(lines, gpd.geodataframe.GeoDataFrame)
#     assert len(lines) == 4
#     assert lines.crs == 'EPSG:4326'


# Testing interpolate_strike_lines
###########################################################
# @pytest.mark.parametrize("gdf_lines5", [gdf_lines5])
# def test_interpolate_strike_lines(gdf_lines5):
#     from gemgis.utils import interpolate_strike_lines
#
#     lines = interpolate_strike_lines(gdf_lines5, 50)
#
#     assert isinstance(lines, gpd.geodataframe.GeoDataFrame)
#     assert lines.crs == 'EPSG:4326'
#     assert len(lines) == 33
#     assert {'X', 'Y', 'Z'}.issubset(lines.columns)


# Testing extract_boreholes
###########################################################
# @pytest.mark.parametrize("gdf_interfaces1_lines", [gdf_interfaces1_lines])
# @pytest.mark.parametrize("gdf_orientations1", [gdf_orientations1])
# @pytest.mark.parametrize("dem",
#                          [
#                              rasterio.open('../docs/getting_started/tutorial/data/test_gemgis/raster1.tif')
#                          ])
# def test_extract_borehole(gdf_interfaces1_lines, gdf_orientations1, dem):
#     from gemgis.postprocessing import extract_borehole
#
#     geo_data = gg.GemPyData(model_name='Model1',
#                             crs='EPSG:4326')
#
#     geo_data.set_extent(-0.0, 972.0, -0.0, 1069.0, 300, 800)
#     geo_data.set_resolution(50, 50, 50)
#
#     interfaces_coords = gg.vector.extract_xyz(gdf_interfaces1_lines, dem, extent=geo_data.extent)
#     geo_data.to_gempy_df(interfaces_coords, 'interfaces')
#
#     orientations_coords = gg.vector.extract_xyz(gdf_orientations1, dem, extent=geo_data.extent)
#     geo_data.to_gempy_df(orientations_coords, 'orientations')
#
#     geo_data.stack = {"Strat_Series": ('Sand1', 'Ton')}
#
#     geo_model = gp.create_model(geo_data.model_name)
#
#     gp.init_data(geo_model, geo_data.extent, geo_data.resolution,
#                  surface_points_df=geo_data.interfaces,
#                  orientations_df=geo_data.orientations,
#                  default_values=True)
#
#     gp.map_stack_to_surfaces(geo_model,
#                              geo_data.stack,
#                              remove_unused_series=True)
#     geo_model.add_surfaces('basement')
#
#     geo_model.set_topography(
#         source='gdal', filepath='../docs/getting_started/tutorial/data/test_gemgis/raster1.tif')
#
#     gp.set_interpolator(geo_model,
#                         compile_theano=True,
#                         theano_optimizer='fast_compile',
#                         verbose=[],
#                         update_kriging=False
#                         )
#
#     gp.compute_model(geo_model, compute_mesh=True)
#
#     sol, well_model, depth_dict = extract_borehole(geo_model, geo_data, [500, 500])
#
#     assert depth_dict == {1: 460.0, 2: 400.0, 3: 300.0}
