# Author: ZJW
# Date: 25th Nov. 2014

import sys
import os
import numpy as np

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

option = [2, 1]  # Tricubic Hermite, with unit scaling
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
 strainFieldUserNumber,
 stressFieldUserNumber) = range(1, 23)

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
GeometricFieldExport(region, "cellML1_undeformed")

# Set up fibre field.
fibreField = FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure,
                             8, [0, 0, 0])

# Set up material field
params = [2.0, 6.0]
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
DependentFieldInitialise(dependentField, geometricField, -14.0)

# Set up CellML constitutive model linkage
cellMLfilename = "mooney_rivlin.xml"
parameter_names = ["c1", "c2"]
[cellML, CellMLModelsField, CellMLParametersField, CellMLIntermediateField] = CellMLSetUp(cellMLUserNumber,
                                                                                          cellMLModelsFieldUserNumber,
                                                                                          cellMLParametersFieldUserNumber,
                                                                                          cellMLIntermediateFieldUserNumber,
                                                                                          region, materialField,
                                                                                          dependentField, parameter_names,
                                                                                          cellMLfilename, option)

# Set up problem
cellMLOptionFull = [cellMLOption, cellML]
print "----> Setting up equations and problem <----\n"
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

solverEquations = BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm,
                                   optionBC, fixXFace, fixYFace, fixZFace, numNodes, option)

# Solve problem
problem.Solve()

# Extract strain field at gauss points from dependent field variable type U1.
print "----> Setting up strain field <----\n"
strainField = StrainFieldSetUp(strainFieldUserNumber, region, decomposition, geometricField, equationsSet, option)

# Export stresses and strains
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

#ExportStressStrain(elementNumbers, xiPositions,cellML, equationsSet, 'problem1_tricubic_CellML_strain.exdata',
                   #'problem1_tricubic_CellML_stress2PK.exdata','problem1_tricubic_CellML_stressCauchy.exdata',
                   #'Strain', 'Stress', [True])

[strain, stress2PK, stressCauchy] = equationsSet.StrainInterpolateXi(1, [0.5, 0.5,0.5], cellML)

calculated2PKTensor = matrixFromSymmetricComponents(stress2PK)
print(calculated2PKTensor)
"""
[strain, stress] = equationsSet.StressStrainEvaluateGaussCellML(elementNumbers[0], 1)
calculatedStrainTensor = matrixFromSymmetricComponents(strain)
calculated2PKTensor = matrixFromSymmetricComponents(stress)
print("Calculated Cauchy strain at gauss point 0.5, 0.5, 0.5:")
print(calculatedStrainTensor)
print("Calculated 2nd Piola-Kirchhoff stress gauss point 0.5, 0.5, 0.5:")
print(calculated2PKTensor)
"""
print "----> Copied over strains from dependent field to strain field.  <----\n"
cellML.Destroy()
CellMLIntermediateField.Destroy()
CellMLParametersField.Destroy()
CellMLModelsField.Destroy()
# Export results
filename = "cellML1"
ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option)

