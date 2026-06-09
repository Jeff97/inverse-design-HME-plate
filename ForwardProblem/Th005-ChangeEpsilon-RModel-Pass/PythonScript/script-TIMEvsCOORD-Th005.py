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

    pth = session.Path(name='Path-Bottom', type=NODE_LIST, expression=(('PART-1-1', (4421, 
    15053, 4399, 15000, 4377, 14947, 4355, 14894, 4333, 14841, 4311, 14788, 
    4289, 14735, 4267, 14682, 4245, 14629, 4223, 14576, 4201, 14523, 4179, 
    14470, 4157, 14417, 4135, 14364, 4113, 14311, 4091, 14258, 4069, 14205, 
    4047, 14152, 4025, 14099, 4003, 14046, 3981, 13993, 3959, 13940, 3937, 
    13887, 3915, 13834, 3893, 13781, 3871, 13728, 3849, 13675, 3827, 13622, 
    3805, 13569, 3783, 13516, 3761, 13463, 3739, 13410, 3717, 13357, 3695, 
    13304, 3673, 13251, 3651, 13198, 3629, 13145, 3607, 13092, 3585, 13039, 
    3563, 12986, 3541, 12933, 3519, 12880, 3497, 12827, 3475, 12774, 3453, 
    12721, 3431, 12668, 3409, 12615, 3387, 12562, 3365, 12509, 3343, 12456, 
    3321, 12403, 3299, 12350, 3277, 12297, 3255, 12244, 3233, 12191, 3211, 
    12138, 3189, 12085, 3167, 12032, 3145, 11979, 3123, 11926, 3101, 11873, 
    3079, 11820, 3057, 11767, 3035, 11714, 3013, 11661, 2991, 11608, 2969, 
    11555, 2947, 11502, 2925, 11449, 2903, 11396, 2881, 11343, 2859, 11290, 
    2837, 11237, 2815, 11184, 2793, 11131, 2771, 11078, 2749, 11025, 2727, 
    10972, 2705, 10919, 2683, 10866, 2661, 10813, 2639, 10760, 2617, 10707, 
    2595, 10654, 2573, 10601, 2551, 10548, 2529, 10495, 2507, 10442, 2485, 
    10389, 2463, 10336, 2441, 10283, 2419, 10230, 2397, 10177, 2375, 10124, 
    2353, 10071, 2331, 10018, 2309, 9965, 2287, 9912, 2265, 9859, 2243, 9806, 
    2221, 9753, 2199, 9700, 2177, 9647, 2155, 9594, 2133, 9541, 2111, 9488, 
    2089, 9435, 2067, 9382, 2045, 9329, 2023, 9276, 2001, 9223, 1979, 9170, 
    1957, 9117, 1935, 9064, 1913, 9011, 1891, 8958, 1869, 8905, 1847, 8852, 
    1825, 8799, 1803, 8746, 1781, 8693, 1759, 8640, 1737, 8587, 1715, 8534, 
    1693, 8481, 1671, 8428, 1649, 8375, 1627, 8322, 1605, 8269, 1583, 8216, 
    1561, 8163, 1539, 8110, 1517, 8057, 1495, 8004, 1473, 7951, 1451, 7898, 
    1429, 7845, 1407, 7792, 1385, 7739, 1363, 7686, 1341, 7633, 1319, 7580, 
    1297, 7527, 1275, 7474, 1253, 7421, 1231, 7368, 1209, 7315, 1187, 7262, 
    1165, 7209, 1143, 7156, 1121, 7103, 1099, 7050, 1077, 6997, 1055, 6944, 
    1033, 6891, 1011, 6838, 989, 6785, 967, 6732, 945, 6679, 923, 6626, 901, 
    6573, 879, 6520, 857, 6467, 835, 6414, 813, 6361, 791, 6308, 769, 6255, 
    747, 6202, 725, 6149, 703, 6096, 681, 6043, 659, 5990, 637, 5937, 615, 
    5884, 593, 5831, 571, 5778, 549, 5725, 527, 5672, 505, 5619, 483, 5566, 
    461, 5513, 439, 5460, 417, 5407, 395, 5354, 373, 5301, 351, 5248, 329, 
    5195, 307, 5142, 285, 5089, 263, 5036, 241, 4983, 219, 4930, 197, 4877, 
    175, 4824, 153, 4771, 131, 4718, 109, 4665, 87, 4612, 65, 4559, 43, 4506, 
    21, )) ))


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
odb_path = 'T:/Abaqus-Temp/20260110MagnetoBeam/ForwardProblem/Th005-ChangeEpsilon-RModel-Pass/'

# Define the file path for the output
folder_name = os.path.join(odb_path, 'Output/')
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
process_odb_file('Beam-Th005-C0-1.odb', 'Beam-Th005-C0-1.csv')
process_odb_file('Beam-Th005-C0-2.odb', 'Beam-Th005-C0-2.csv')
process_odb_file('Beam-Th005-C0-3.odb', 'Beam-Th005-C0-3.csv')
process_odb_file('Beam-Th005-C0-4.odb', 'Beam-Th005-C0-4.csv')
process_odb_file('Beam-Th005-C0-5.odb', 'Beam-Th005-C0-5.csv')