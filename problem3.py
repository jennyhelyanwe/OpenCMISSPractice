# Problem 3. Force loading, tricubic hermite.
# Author: ZJW
# Date: 22nd Oct. 2014

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

numberGlobalXElements = 1
numberGlobalYElements = 1
numberGlobalZElements = 1
numOfXi = 3

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

# Set up mesh
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber, region, numOfXi)
mesh.NumberOfComponentsSet(2)
mesh.NumberOfElementsSet(1)

nodes = CMISS.Nodes()
nodes.CreateStart(region, 8)
nodes.CreateFinish()

# Cubic Hermite component
cubicMeshComponentNumber = 1
cubicElem = CMISS.MeshElements()
cubicElem.CreateStart(mesh, cubicMeshComponentNumber, cubicBasis)
cubicElem.NodesSet(1, [1, 2, 3, 4, 5, 6, 7, 8])
cubicElem.CreateFinish()

# Linear lagrange component
linearMeshComponentNumber = 2
linearElem = CMISS.MeshElements()
linearElem.CreateStart(mesh, linearMeshComponentNumber, linearBasis)
linearElem.NodesSet(1, [1, 2, 3, 4, 5, 6, 7, 8])
linearElem.CreateFinish()

mesh.CreateFinish()

# Set up decomposition for the mesh. 
decomposition = DecompositionSetUp(decompositionUserNumber, mesh, numOfCompNodes)

# Set up geometric field. 
geometricField = GeometricFieldSetUp(geometricFieldUserNumber, region, decomposition, option)

# Initialise geometric field.
xNodes = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
yNodes = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
zNodes = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
geometricField = GeometricFieldInitialise(xNodes, yNodes, zNodes, geometricField, 8, option)

# Export undeformed geometry. 
GeometricFieldExport(region, "problem3_undeformed")

# Set up fibre field.
fibreField = FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure, 8,
                             [0, 0, 0])

# Set up material field.
params = [2.0, 6.0]
materialField = MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, 2, option, params)

# Set up dependent field.
dependentField = DependentFieldSetUp(dependentFieldUserNumber, region, decomposition, geometricField, option)

# Initialise dependent field. 
DependentFieldInitialise(dependentField, geometricField, -8.0)

# Set up equations set. 
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField, CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.MOONEY_RIVLIN,
                         equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet = EquationsSetUp(equationsSet, materialFieldUserNumber, materialField, dependentFieldUserNumber,
                              dependentField)

# Set up problem
[problem, solverEquations] = ProblemAndSolverSetUp(problemUserNumber, equationsSet, 1e-5)

# Initialise BC
appliedFace = [2, 4, 6, 8]
faceNormal = 1
appliedDirection = 1  # X direction.
increm = 2.0
optionBC = 2  # Force.

fixXFace = [1, 3, 5, 7]
fixYFace = [1, 2, 5, 6]
fixZFace = [1, 2, 3, 4]

solverEquations = BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm,
                                   optionBC, fixXFace, fixYFace, fixZFace, 8, option)

# Solve problem
problem.Solve()

# Export results
filename = "problem3"
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)

