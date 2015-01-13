# Problem 1. Uniaxial extension, trilinear.
# Author: ZJW
# Date: 20 Oct. 2014

import sys
import os
import numpy as np

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
cellMLOption = [False]

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
    deformedFieldUserNumber,
    strainFieldUserNumber) = range(1, 16)

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
params = [2.0, 6.0]
materialField = MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, params, option,
                                   cellMLOption)

# Set up equations set
equationsSetField = CMISS.Field()  # Equations are also in a field
equationsSet = CMISS.EquationsSet()  # Initialise an equation set.
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY, CMISS.EquationsSetSubtypes.MOONEY_RIVLIN,
                         equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()
print "----> Set up equations set <---\n"

# Set up material field in equations set.
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()

# Set up dependent field.
[dependentField, equationsSet] = DependentFieldSetUp(dependentFieldUserNumber, equationsSet, option, cellMLOption)

# Initialise dependent field.
DependentFieldInitialise(dependentField, geometricField, -8.0)


# Set up strain field
strainField = StrainFieldSetUp(strainFieldUserNumber, region, decomposition, geometricField, equationsSet)

# Set up problem
[problem, solverEquations] = EquationsProblemSolverSetUp(problemUserNumber, equationsSet, cellMLOption)

# Initialise BC
appliedFace = [2, 4, 6, 8]
faceNormal = 1  # Face normal vector aligns with x direction.
appliedDirection = 1  # X direction.
increm = 0.1
optionBC = 1  # Compression/extension.

fixXFace = [1, 3, 5, 7]
fixYFace = [1, 2, 5, 6]
fixZFace = [1, 2, 3, 4]

BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm,
                                   optionBC, fixXFace, fixYFace, fixZFace, numNodes, option)

# Solve Problem
problem.Solve()

elementNumbers = [1]
xiPositions = [[0.25, 0.25, 0.25],
    [0.25, 0.75, 0.25],
    [0.25, 0.25, 0.75],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.25],
    [0.75, 0.75, 0.25],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.75]]
xiPositions = np.array(xiPositions)
ExportStressStrain(elementNumbers,xiPositions, equationsSet, 'problem1_strain.exdata', 'problem1_stress.exdata',
                   'Strain', 'Stress', cellMLOption)


xiPosition = [0.5, 0.5, 0.5]
elementNumber = 1
[calculatedStrain, PK2_tensor] = equationsSet.StrainInterpolateXi(elementNumber, xiPosition)

calculatedStrainTensor = matrixFromSymmetricComponents(calculatedStrain)
calculated2PKTensor = matrixFromSymmetricComponents(PK2_tensor)
print("Calculated Cauchy strain at gauss point 0.5, 0.5, 0.5:")
print(calculatedStrainTensor)
print("Calculated 2nd Piola-Kirchhoff stress gauss point 0.5, 0.5, 0.5:")
print(calculated2PKTensor)

# Export results
filename = "problem1"
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)


