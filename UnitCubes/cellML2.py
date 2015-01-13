# # Problem 6. Pressure loading, transversely isotropic, heterogeneous fibre.
# Author: ZJW
# Date: 29th Oct. 2014

import sys
import os

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))
### Initialise OpenCMISS
from opencmiss import CMISS
from math import pi
from numpy import array
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

option = [1]  # Tricubic Hermite, with unit scaling
microstructure = 1  # Homogeneous fibre angles
cellMLOption = [True]
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
 problemUserNumber,
 cellMLUserNumber,
 cellMLModelsFieldUserNumber,
 cellMLParametersFieldUserNumber,
 cellMLIntermediateFieldUserNumber,
 strainFieldUserNumber) = range(1, 22)

# Set up region and CS
[numOfCompNodes, compNodeNum, CS, region] = BasicSetUp(regionUserNumber, coordinateSystemUserNumber)

# Set up tricubic Hermite basis functions
basisUserNumber = linearBasisUserNumber
basis = BasisFunction(basisUserNumber, numOfXi, option)
#bases = [linearBasis, cubicBasis]

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

# Initialise geometric field.
xNodes = [0.0, length, 0.0, length, 0.0, length, 0.0, length]
yNodes = [0.0, 0.0, width, width, 0.0, 0.0, width, width]
zNodes = [0.0, 0.0, 0.0, 0.0, height, height, height, height]
GeometricFieldInitialise(xNodes, yNodes, zNodes, geometricField, 8, option)
# Export undeformed geometry. 
GeometricFieldExport(region, "cellML2_undeformed")

# Set up fibre field.
fibreAngle =[0*pi/180, 0, 0]
fibreField = FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure,
                             numNodes, fibreAngle)

# Set up material field
params = [1.0, 5.0, 10.0, 0.0, 5.0]
materialField = MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, params, option,
                                   cellMLOption)

# Set up equations set.
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField, CMISS.EquationsSetClasses.ELASTICITY,
                         CMISS.EquationsSetTypes.FINITE_ELASTICITY,
                         CMISS.EquationsSetSubtypes.CONSTITUTIVE_LAW_IN_CELLML_EVALUATE, equationsSetFieldUserNumber,
                         equationsSetField)
equationsSet.CreateFinish()
print "----> Set up equations set <---\n"

# Set up material field in equations set.
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()

# Set up dependent field.
[dependentField, equationsSet] = DependentFieldSetUp(dependentFieldUserNumber, equationsSet, option, cellMLOption)

# Initialise dependent field.
DependentFieldInitialise(dependentField, geometricField, -2.5)

# Set up CellML constitutive model linkage
cellMLfilename = "transversely_isotropic_exponential.cellml"
parameter_names = ["C1", "C2", "C3", "C4", "C5"]
[cellML, CellMLModelsField, CellMLParametersField, CellMLIntermediateField] = CellMLSetUp(cellMLUserNumber,
                                                                                          cellMLModelsFieldUserNumber,
                                                                                          cellMLParametersFieldUserNumber,
                                                                                          cellMLIntermediateFieldUserNumber,
                                                                                          region, materialField,
                                                                                          dependentField, parameter_names,
                                                                                          cellMLfilename, option)


# Set up problem
cellMLOptionFull = [cellMLOption, cellML]
[problem, solverEquations] = EquationsProblemSolverSetUp(problemUserNumber, equationsSet, cellMLOptionFull)

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

# Solve problem
problem.Solve()

cellML.Destroy()
CellMLIntermediateField.Destroy()
CellMLParametersField.Destroy()
CellMLModelsField.Destroy()

# Extract strain field at gauss points from dependent field variable type U1.
print "----> Setting up strain field <----\n"
strainField = StrainFieldSetUp(strainFieldUserNumber, region, decomposition, geometricField, equationsSet, option)

elementNumbers = [1]
gaussPoints = [1,8,64]
ExportStressStrain(elementNumbers,gaussPoints, equationsSet, 'cellML2_strain.exdata', 'cellML2_stress.exdata',
                   'Strain', 'Stress', cellMLOption)

[strain, stress] = equationsSet.StressStrainEvaluateGaussCellML(elementNumbers[0], 1)
calculatedStrainTensor = matrixFromSymmetricComponents(strain)
calculated2PKTensor = matrixFromSymmetricComponents(stress)
print("Calculated Cauchy strain at gauss point 1:")
print(calculatedStrainTensor)
print("Calculated 2nd Piola-Kirchhoff stress gauss point 1:")
print(calculated2PKTensor)

# Export results
filename = "cellML2"
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)

