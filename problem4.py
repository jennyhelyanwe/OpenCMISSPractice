# # Problem 4. Pressure loading, tricubic Hermite.
# Author: ZJW
# Date: 28th Oct. 2014

import sys
import os

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))
### Initialise OpenCMISS
from opencmiss import CMISS
from lib import *

### Set problem parameters
height = 1.0
width = 1.0
length = 1.0
dimensions = [height, width, length]

numberGlobalXElements = 1
numberGlobalYElements = 1
numberGlobalZElements = 1

elements = [numberGlobalXElements, numberGlobalYElements, numberGlobalZElements]

numOfXi = 3
numNodes = 8

option = [2, 1]  # Tricubic Hermite, with unit scaling
microstructure = 1  # Homogeneous fibre angles

### Set arbitrary user numbers which are unique to each object.
(coordinateSystemUserNumber,
 regionUserNumber,
 linearBasisUserNumber,
 cubicBasisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 fibreFieldUserNumber,
 materialFieldUserNumber,
 dependentFieldUserNumber,
 deformedFieldUserNumber,
 equationsSetFieldUserNumber,
 equationsSetUserNumber,
 equationsUserNumber,
 problemUserNumber) = range(1, 17)

# Set up region and CS
[numOfCompNodes, compNodeNum, CS, region] = BasicSetUp(regionUserNumber, coordinateSystemUserNumber)

# Set up tricubic Hermite basis functions
basisUserNumber = [linearBasisUserNumber, cubicBasisUserNumber]
[linearBasis, cubicBasis] = BasisFunction(basisUserNumber, numOfXi, option)
#bases = [linearBasis, cubicBasis]
bases = [cubicBasis, linearBasis]
# Set up generated mesh
[generatedMesh, mesh] = GeneratedMesh(generatedMeshUserNumber, meshUserNumber, region, bases, dimensions, elements)

# Set up decomposition for the mesh. 
decomposition = DecompositionSetUp(decompositionUserNumber, mesh, numOfCompNodes)

# Set up geometric field. 
geometricField = GeometricFieldSetUp(geometricFieldUserNumber, region, decomposition, option)

# Initialise geometric field.
generatedMesh.GeometricParametersCalculate(geometricField)

# Export undeformed geometry. 
GeometricFieldExport(region, "problem4_undeformed")

# Set up fibre field.
fibreField = FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure,
                             numNodes, [0, 0, 0])

# Set up material field.
params = [2.0, 6.0]
numParams = 2
materialField = MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, numParams, option,
                                   params)

# Set up dependent field.
dependentField = DependentFieldSetUp(dependentFieldUserNumber, region, decomposition, geometricField, option)

# Initialise dependent field. 
DependentFieldInitialise(dependentField, geometricField, -8.0)

# Set up equations set. 
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.MOONEY_RIVLIN,
                         equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet = EquationsSetUp(equationsSet, materialFieldUserNumber, materialField, dependentFieldUserNumber,
                              dependentField)

# Set up problem
[problem, solverEquations] = ProblemAndSolverSetUp(problemUserNumber, equationsSet, 1e-10)

# Initialise BC
appliedFace = [2, 4, 6, 8]
faceNormal = 1
appliedDirection = 1  # X direction.
increm = 8.0
optionBC = 3  # Pressure.

fixXFace = [1, 3, 5, 7]
fixYFace = [1, 3, 5, 7]
fixZFace = [1, 3, 5, 7]

solverEquations = BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm,
                                   optionBC, fixXFace, fixYFace, fixZFace, numNodes, option)

# Solve problem
problem.Solve()

# Export results
filename = "problem4"
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)

