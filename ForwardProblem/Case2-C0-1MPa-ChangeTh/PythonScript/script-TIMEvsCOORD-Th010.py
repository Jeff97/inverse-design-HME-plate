from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os

# Function to process a single ODB file
def process_odb_file(odb_file_name, output_file_name):
    # Open the ODB file
    o1 = session.openOdb(name=odb_path + odb_file_name)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].makeCurrent()

    odb = session.odbs[odb_path + odb_file_name]
    
    # Get the data for the selected node

    pth = session.Path(name='Path-Bottom', type=NODE_LIST, expression=(('PART-1-1', (4260,14595,4230,14522,4200,14449,4170,14376,4140,14303,4110,14230,4080,14157,4050,14084,4020,14011,3990,13938,3960,13865,3930,13792,3900,13719,3870,13646,3840,13573,3810,13500,3780,13427,3750,13354,3720,13281,3690,13208,3660,13135,3630,13062,3600,12989,3570,12916,3540,12843,3510,12770,3480,12697,3450,12624,3420,12551,3390,12478,3360,12405,3330,12332,3300,12259,3270,12186,3240,12113,3210,12040,3180,11967,3150,11894,3120,11821,3090,11748,3060,11675,3030,11602,3000,11529,2970,11456,2940,11383,2910,11310,2880,11237,2850,11164,2820,11091,2790,11018,2760,10945,2730,10872,2700,10799,2670,10726,2640,10653,2610,10580,2580,10507,2550,10434,2520,10361,2490,10288,2460,10215,2430,10142,2400,10069,2370,9996,2340,9923,2310,9850,2280,9777,2250,9704,2220,9631,2190,9558,2160,9485,2130,9412,2100,9339,2070,9266,2040,9193,2010,9120,1980,9047,1950,8974,1920,8901,1890,8828,1860,8755,1830,8682,1800,8609,1770,8536,1740,8463,1710,8390,1680,8317,1650,8244,1620,8171,1590,8098,1560,8025,1530,7952,1500,7879,1470,7806,1440,7733,1410,7660,1380,7587,1350,7514,1320,7441,1290,7368,1260,7295,1230,7222,1200,7149,1170,7076,1140,7003,1110,6930,1080,6857,1050,6784,1020,6711,990,6638,960,6565,930,6492,900,6419,870,6346,840,6273,810,6200,780,6127,750,6054,720,5981,690,5908,660,5835,630,5762,600,5689,570,5616,540,5543,510,5470,480,5397,450,5324,420,5251,390,5178,360,5105,330,5032,300,4959,270,4886,240,4813,210,4740,180,4667,150,4594,120,4521,90,4448,60,4375,30, )) ))


    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 
    'COOR1'), )

    COORDXAtBot = session.XYDataFromPath(path=pth, name='X-Coords',
    includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS,
    numIntervals=10, projectionTolerance=0, shape=DEFORMED, 
    labelType=X_COORDINATE, 
    removeDuplicateXYPairs=True, includeAllElements=False)

    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 
    'COOR3'), )

    COORDZAtBot = session.XYDataFromPath(path=pth, name='Z-Coords',
    includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS,
    numIntervals=10, projectionTolerance=0, shape=DEFORMED, 
    labelType=X_COORDINATE, 
    removeDuplicateXYPairs=True, includeAllElements=False)
    
    # Write the TIMEvsCOORD to a local file
    output_file_path = os.path.join(folder_name, output_file_name)
    exp_data = (COORDXAtBot, COORDZAtBot)
    
    # Use the session's writeXYReport method to save the data to a CSV file
    session.writeXYReport(fileName=output_file_path, xyData=exp_data)
    
    # Close the ODB file to release resources
    odb.close()

# Process each ODB file ===============================================

# Path to your ODB files folder
odb_path = 'T:/Abaqus-Temp/20260110MagnetoBeam/Case2-C0-1MPa-ChangeTh/'

# Define the file path for the output
folder_name = os.path.join(odb_path, 'Output/')
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
process_odb_file('Beam-C0-1-Th010.odb', 'Beam-C0-1-Th010.csv')