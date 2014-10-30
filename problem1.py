# Problem 1. Uniaxial extension, trilinear.
# Author: ZJW
# Date: 20 Oct. 2014

import sys
import os

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))
# Initialise OpenCMISS
from lib import *
from numpy import array

### Set problem parameters
height = 1.0
width = 1.0
length = 1.0

numberGlobalXElements = 1
numberGlobalYElements = 1
numberGlobalZElements = 1
numOfXi = 3
numNodes = 8

option = [1] # Trilinear
microstructure = 1 # Homogeneous fibre angles. 

### Set arbitrary user numbers which are unique to each object.
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    problemUserNumber,
    fibreFieldUserNumber,
    materialFieldUserNumber,
    deformedFieldUserNumber) = range(1, 15)

# Set up region and CS
[numOfCompNodes, compNodeNum, CS, region] = BasicSetUp(regionUserNumber, coordinateSystemUserNumber)

# Set up trilinear lagrange basis function
basis = BasisFunction(basisUserNumber, numOfXi, option)

# Create a mesh
mesh = CMISS.Mesh()  # Initialise an instance of mesh
mesh.CreateStart(meshUserNumber, region, numOfXi)  # Create mesh in the specified region with 3 elemental directions.
mesh.NumberOfComponentsSet(1)  # Set the number of components of the mesh
mesh.NumberOfElementsSet(1)  # Set the total number of elements in the mesh

nodes = CMISS.Nodes()  # Initialise an instance of nodes object
nodes.CreateStart(region, 8)  # Create 8 nodes
nodes.CreateFinish()  # Finish creation

elem = CMISS.MeshElements()  # Initialise an instance of mesh elements.
elem.CreateStart(mesh, 1, basis)  # Create an element in mesh component 1.
elem.NodesSet(1, [1, 2, 3, 4, 5, 6, 7, 8])  # The nodes in global element 1 are given user node numbers 1 to 8.
elem.CreateFinish()
mesh.CreateFinish()

# Set up decomposition for the mesh. 
decomposition = DecompositionSetUp(decompositionUserNumber, mesh, numOfCompNodes)

# Set up geometric field. 
geometricField = GeometricFieldSetUp(geometricFieldUserNumber, region, decomposition, option)


### Update the geometric field parameters manually
xNodes = [0.0, length, 0.0, length, 0.0, length, 0.0, length]
yNodes = [0.0, 0.0, width, width, 0.0, 0.0, width, width]
zNodes = [0.0, 0.0, 0.0, 0.0, height, height, height, height]
GeometricFieldInitialise(xNodes, yNodes, zNodes, geometricField, 8, option)

# Export undeformed geometry. 
GeometricFieldExport(region, "problem1_undeformed")

# Set up fibre field.
fibreField = FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure,
                             8, [0, 0, 0])

# Set up material field.
params = [2.0, 6.0]  # kPa
materialField = MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, 2, option, params)

# Set up dependent field.
dependentField = DependentFieldSetUp(dependentFieldUserNumber, region, decomposition, geometricField, option)

# Initialise dependent field. 
DependentFieldInitialise(dependentField, geometricField, -8.0)

# Set up equations set
equationsSetField = CMISS.Field()  # Equations are also in a field
equationsSet = CMISS.EquationsSet()  # Initialise an equation set.
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.MOONEY_RIVLIN,
                         equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet = EquationsSetUp(equationsSet, materialFieldUserNumber, materialField, dependentFieldUserNumber,
                              dependentField)

# Set up problem
[problem, solverEquations] = ProblemAndSolverSetUp(problemUserNumber, equationsSet, 1e-5)

# Initialise BC
appliedFace = [2, 4, 6, 8]
faceNormal = 1  # Face normal vector aligns with x direction.
appliedDirection = 1  # X direction.
increm = 1.1
optionBC = 1  # Compression/extension.

fixXFace = [1, 3, 5, 7]
fixYFace = [1, 2, 5, 6]
fixZFace = [1, 2, 3, 4]

BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm,
                                   optionBC, fixXFace, fixYFace, fixZFace, numNodes, option)

# Solve Problem
problem.Solve()

# Export results
filename = "problem1"
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)


