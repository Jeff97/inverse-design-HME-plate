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

    pth = session.Path(name='Path-Bottom', type=NODE_LIST, expression=(('PART-1-1', (4498,15370,4472,15307,4446,15244,4420,15181,4394,15118,4368,15055,4342,14992,4316,14929,4290,14866,4264,14803,4238,14740,4212,14677,4186,14614,4160,14551,4134,14488,4108,14425,4082,14362,4056,14299,4030,14236,4004,14173,3978,14110,3952,14047,3926,13984,3900,13921,3874,13858,3848,13795,3822,13732,3796,13669,3770,13606,3744,13543,3718,13480,3692,13417,3666,13354,3640,13291,3614,13228,3588,13165,3562,13102,3536,13039,3510,12976,3484,12913,3458,12850,3432,12787,3406,12724,3380,12661,3354,12598,3328,12535,3302,12472,3276,12409,3250,12346,3224,12283,3198,12220,3172,12157,3146,12094,3120,12031,3094,11968,3068,11905,3042,11842,3016,11779,2990,11716,2964,11653,2938,11590,2912,11527,2886,11464,2860,11401,2834,11338,2808,11275,2782,11212,2756,11149,2730,11086,2704,11023,2678,10960,2652,10897,2626,10834,2600,10771,2574,10708,2548,10645,2522,10582,2496,10519,2470,10456,2444,10393,2418,10330,2392,10267,2366,10204,2340,10141,2314,10078,2288,10015,2262,9952,2236,9889,2210,9826,2184,9763,2158,9700,2132,9637,2106,9574,2080,9511,2054,9448,2028,9385,2002,9322,1976,9259,1950,9196,1924,9133,1898,9070,1872,9007,1846,8944,1820,8881,1794,8818,1768,8755,1742,8692,1716,8629,1690,8566,1664,8503,1638,8440,1612,8377,1586,8314,1560,8251,1534,8188,1508,8125,1482,8062,1456,7999,1430,7936,1404,7873,1378,7810,1352,7747,1326,7684,1300,7621,1274,7558,1248,7495,1222,7432,1196,7369,1170,7306,1144,7243,1118,7180,1092,7117,1066,7054,1040,6991,1014,6928,988,6865,962,6802,936,6739,910,6676,884,6613,858,6550,832,6487,806,6424,780,6361,754,6298,728,6235,702,6172,676,6109,650,6046,624,5983,598,5920,572,5857,546,5794,520,5731,494,5668,468,5605,442,5542,416,5479,390,5416,364,5353,338,5290,312,5227,286,5164,260,5101,234,5038,208,4975,182,4912,156,4849,130,4786,104,4723,78,4660,52,4597,26, )) ))


    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 
    'COOR1'), )

    COORDXAtBot = session.XYDataFromPath(path=pth, name='X-Coords',
    includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS,
    numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, 
    labelType=X_COORDINATE, 
    removeDuplicateXYPairs=True, includeAllElements=False)

    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 
    'COOR3'), )

    COORDZAtBot = session.XYDataFromPath(path=pth, name='Z-Coords',
    includeIntersections=False, projectOntoMesh=False, pathStyle=PATH_POINTS,
    numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, 
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
odb_path = 'T:/Abaqus-Temp/20260110MagnetoBeam/Benchmark-Th-Epsilon003/'

# Define the file path for the output
folder_name = os.path.join(odb_path, 'Output/')
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
process_odb_file('Beam-Th007.odb', 'Beam-Th007.csv')