# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# This python file is a personal library of snippets of python code which allows
# re-use of common set-up commands for solving a problem in OpenCMISS. 

# Each snippet is written with a range of input options and calls the appropriate
# OpenCMISS linked commands to set up the problem. It is essentially a personal high 
# level library that will allow cleaner scripting for solving cardiac mechanics
# simulations, making it easier to debug.  

# Author: Zhinuo Jenny Wang
# Date: 20th October 2014

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

from opencmiss import CMISS
#=================================================================================#
def BasicSetUp(regionUserNumber, coordinateSystemUserNumber):
    # This snippet sets up the world region, 3D CS, parallel computing nodes, and
    # diagnostics. These are oft-repeated commands.

    # Set up diagnostics/debug
    CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.IN, [1, 2, 3, 4, 5], "Diagnostics",
                           ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

    # Get computational node information for parallel computing
    numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
    computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

    #  Set up 3D RC coordinate system
    coordinateSystem = CMISS.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.dimension = 3
    coordinateSystem.CreateFinish()

    # Create world region
    region = CMISS.Region()
    region.CreateStart(regionUserNumber, CMISS.WorldRegion)
    region.label = "Region"
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()

    # Output for diagnostics
    print "----> Set up coordinate system and world region <----\n"

    return numberOfComputationalNodes, computationalNodeNumber, coordinateSystem, region


#=================================================================================#

#=================================================================================#
def BasisFunction(basisUserNumber, numOfXi, option):
    # This snippet sets up the basis function depending on the option given.
    if option[0] == 1:
        # Trilinear Lagrange basis
        basis = CMISS.Basis()
        basis.CreateStart(basisUserNumber)
        basis.numberOfXi = numOfXi
        basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
        basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numOfXi
        basis.quadratureNumberOfGaussXi = [2] * numOfXi
        basis.CreateFinish()
        # Output for diagnostics
        print "----> Set up trilinear basis functions <----\n"
        return basis
    elif option[0] == 2:
        # Tricubic Hermite basis
        linearBasis = CMISS.Basis()
        linearBasis.CreateStart(basisUserNumber[0])
        linearBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numOfXi)
        linearBasis.QuadratureNumberOfGaussXiSet([4] * numOfXi)
        linearBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        linearBasis.CreateFinish()

        cubicBasis = CMISS.Basis()  # For geometry.
        cubicBasis.CreateStart(basisUserNumber[1])
        cubicBasis.InterpolationXiSet([CMISS.BasisInterpolationSpecifications.CUBIC_HERMITE] * numOfXi)
        cubicBasis.QuadratureNumberOfGaussXiSet([4] * numOfXi)
        cubicBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        cubicBasis.CreateFinish()
        # Output for diagnostics
        print "----> Set up trilinear and tricubic hermite basis functions <----\n"
        return linearBasis, cubicBasis


#=================================================================================#

#=================================================================================#
def GeneratedMesh(generatedMeshUserNumber, meshUserNumber, region, bases, dimensions, elements):
    generatedMesh = CMISS.GeneratedMesh()
    generatedMesh.CreateStart(generatedMeshUserNumber, region)
    generatedMesh.TypeSet(CMISS.GeneratedMeshTypes.REGULAR)
    generatedMesh.BasisSet(bases)
    generatedMesh.ExtentSet(dimensions)
    generatedMesh.NumberOfElementsSet(elements)
    mesh = CMISS.Mesh()
    generatedMesh.CreateFinish(meshUserNumber, mesh)

    return generatedMesh, mesh


#=================================================================================#

#=================================================================================#
def DecompositionSetUp(decompositionUserNumber, mesh, numberOfComputationalNodes):
    # Set up decomposition
    decomposition = CMISS.Decomposition()
    decomposition.CreateStart(decompositionUserNumber, mesh)
    decomposition.type = CMISS.DecompositionTypes.CALCULATED
    decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    decomposition.CalculateFacesSet(True)
    decomposition.CreateFinish()

    # Output for diagnostics
    print "----> Set up decomposition <----\n"
    return decomposition


#=================================================================================#

#=================================================================================#
def GeometricFieldSetUp(geometricFieldUserNumber, region, decomposition, option):
    # Sets up the geometry field - linear, same interpolation in all 3 directions.
    geometricField = CMISS.Field()  # Initialise
    geometricField.CreateStart(geometricFieldUserNumber, region)
    geometricField.MeshDecompositionSet(decomposition)
    geometricField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
    geometricField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Geometry")
    geometricField.NumberOfVariablesSet(1)
    geometricField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)

    # Always use the first mesh component for geometry.
    for i in [1, 2, 3]:
        geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)

    if option[0] == 2:
        # Tricubic Hermite
        if option[1] == 1:
            geometricField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
            # Output for diagnostics
            print "----> Set up tricubic Hermite geometric field with unit scaling <----\n"
        elif option[1] == 2:
            geometricField.ScalingTypeSet(CMISS.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output for diagnostics
            print "----> Set up tricubic Hermite geometric field with arithmetic mean scaling <----\n"

    geometricField.CreateFinish()

    return geometricField


#=================================================================================#

#=================================================================================#
def GeometricFieldInitialise(xNodes, yNodes, zNodes, geometricField, numNodes, option):
    # Initialise nodal values.
    for node, value in enumerate(xNodes, 1):
        geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,
                                                CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1, value)
    for node, value in enumerate(yNodes, 1):
        geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,
                                                CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2, value)
    for node, value in enumerate(zNodes, 1):
        geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 1,
                                                CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3, value)

    # Output
    print "----> Initialised geometric nodal values <----\n"

    return geometricField


#=================================================================================#

#=================================================================================#
def GeometricFieldExport(region, filename):
    exportField = CMISS.Fields()
    exportField.CreateRegion(region)
    exportField.NodesExport(filename, "FORTRAN")
    exportField.ElementsExport(filename, "FORTRAN")
    exportField.Finalise()

    # Output
    print "----> Export undeformed geometry <----\n"


#=================================================================================#

#=================================================================================#
def FibreFieldSetUp(fibreFieldUserNumber, region, decomposition, geometricField, option, microstructure, numNodes,
                    angles):
    # Sets up the fibre field.
    fibreField = CMISS.Field()
    fibreField.CreateStart(fibreFieldUserNumber, region)
    fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Fibre")

    if option[0] == 2:
        # Tricubic Hermite interpolation
        fibreField.NumberOfVariablesSet(1)
        fibreField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 3)

        if option[1] == 1:
            fibreField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
            # Output
            print "----> Set up tricubic hermite fibre field with unit scaling <----\n"
        elif option[1] == 2:
            fibreField.ScalingTypeSet(CMISS.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output
            print "----> Set up tricubic hermite fibre field with arithmetic mean scaling <----\n"

        if microstructure == 1:
            # Homogeneous fibre field.
            for component in [1, 2, 3]:
                fibreField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, component,
                                                     CMISS.FieldInterpolationTypes.CONSTANT)
        elif microstructure == 2:
            # Heterogeneous fibre field using linear interpolation.
            for component in [1, 2, 3]:
                fibreField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, component, 1)

    fibreField.CreateFinish()

    #####################################################
    if option[0] == 2:
        # Initialise
        if microstructure == 1:
            # Homogeneous fibre field.
            for component, fibreAngle in enumerate(angles, 1):
                fibreField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,
                                                       component, fibreAngle)
            # Output
            print "----> Initialised homogeneous fibre angles <----\n"
        elif microstructure == 2:
            # Inhomogeneous fibre field using linear interpolation.
            for node in range(1, numNodes + 1):
                for component in [1, 2, 3]:
                    fibreField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES,
                                                        1, CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node,
                                                        component, angles[node - 1, component - 1])
            # Output
            print "----> Initialised heterogeneous fibre angles <----\n"

    return fibreField


#=================================================================================#

#=================================================================================#
def MaterialFieldSetUp(materialFieldUserNumber, region, decomposition, geometricField, numParams, option, params):
    # Sets up material field, and apply field to mesh component.
    materialField = CMISS.Field()
    materialField.CreateStart(materialFieldUserNumber, region)
    materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
    materialField.MeshDecompositionSet(decomposition)
    materialField.GeometricFieldSet(geometricField)
    materialField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Material")

    materialField.NumberOfVariablesSet(1)
    materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, numParams)

    if option[0] == 1:
        # Trilinear mesh.
        for i in range(numParams):
            materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i + 1, 1)
            # Output
            print "----> Set up trilinear material field <----\n"
    elif option[0] == 2:
        # Tricubic Hermite
        for i in range(numParams):
            materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i + 1, 1)

            print "Set up material field for component " + str(i + 1) + "\n"
        # Output
        print "----> Set up tricubic hermite material field <----\n"

        # Tricubic Hermite mesh.
        if option[1] == 1:
            materialField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
            # Output
            print "----> Set up tricubic hermite material field with unit scaling <----\n"
        elif option[1] == 2:
            materialField.ScalingTypeSet(CMISS.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output
            print "----> Set up tricubic hermite material field with arithmetic mean scaling <----\n"
    materialField.CreateFinish()

    #########################################################################
    # Initialise parameter values.
    for i in range(numParams):
        materialField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, i + 1,
                                                params[i])


    # Output
    print "----> Initialised " + str(numParams) + " material parameters <----\n"
    return materialField


#=================================================================================#

#=================================================================================#
def DependentFieldSetUp(dependentFieldUserNumber, region, decomposition, geometricField, option):
    # Standard set up dependent field with first derivatives.
    dependentField = CMISS.Field()
    dependentField.CreateStart(dependentFieldUserNumber, region)
    dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)
    dependentField.MeshDecompositionSet(decomposition)
    dependentField.GeometricFieldSet(geometricField)
    dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U, "Dependent")
    dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT)

    dependentField.NumberOfVariablesSet(2)
    dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U, 4)
    dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN, 4)

    if option[0] == 1:
        # Trilinear
        for i in [1, 2, 3]:
            dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
            dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i, 1)

        dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U, 4,
                                                 CMISS.FieldInterpolationTypes.ELEMENT_BASED)
        dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.DELUDELN, 4,
                                                 CMISS.FieldInterpolationTypes.ELEMENT_BASED)
        # Output
        print "----> Interpolate element based hydrostatic pressure <----\n"

    elif option[0] == 2:
        # Tricubic Hermite
        for i in [1, 2, 3]:
            dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, i, 1)
            dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, i, 1)

        dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, 4, 2)
        dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN, 4, 2)
        # Output
        print "----> Interpolate hydrostatic pressure linearly <----\n"
        if option[1] == 1:
            dependentField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
            # Output
            print "----> Set up dependent field with unit scaling <----\n"
        elif option[1] == 2:
            dependentField.ScalingTypeSet(CMISS.FieldScalingTypes.ARITHMETIC_MEAN)
            # Output
            print "----> Set up dependent field with arithmetic mean scaling <----\n"

    dependentField.CreateFinish()

    return dependentField


#=================================================================================#

#=================================================================================#
def DependentFieldInitialise(dependentField, geometricField, hydroInit):
    # Copy over undeformed geometry to initialise dependent field.
    for i in [1, 2, 3]:
        CMISS.Field.ParametersToFieldParametersComponentCopy(geometricField, CMISS.FieldVariableTypes.U,
                                                             CMISS.FieldParameterSetTypes.VALUES, i, dependentField,
                                                             CMISS.FieldVariableTypes.U,
                                                             CMISS.FieldParameterSetTypes.VALUES, i)

    # Output
    print "----> Initialised dependent field with undeformed geometry <----\n"

    # Set hydrostatic pressure initial guess.
    dependentField.ComponentValuesInitialise(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, 4,
                                             hydroInit)

    # Output
    print "----> Initialised hydrostatic pressure guess of " + str(hydroInit) + " <----\n"


#=================================================================================#

#=================================================================================#
def EquationsSetUp(equationsSet, materialFieldUserNumber, materialField, dependentFieldUserNumber, dependentField):
    # Add material and fibre fields to equations set.
    # Add the material field to the equation
    equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
    equationsSet.MaterialsCreateFinish()

    # Add the dependent field to the equation
    equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
    equationsSet.DependentCreateFinish()

    # Create equations from the equation set
    equations = CMISS.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
    equations.outputType = CMISS.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    # Output
    print "----> Set up equations set <----\n"
    return equationsSet


#=================================================================================#

#=================================================================================#
def ProblemAndSolverSetUp(problemUserNumber, equationsSet, TOL):
    # Set up standard options for problem and solvers.
    # Define problem
    problem = CMISS.Problem()
    problem.CreateStart(problemUserNumber)
    problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY, CMISS.ProblemTypes.FINITE_ELASTICITY,
                             CMISS.ProblemSubTypes.NONE)
    problem.CreateFinish()
    # Output
    print "----> Set up problem <----\n"
    # Create control loops
    problem.ControlLoopCreateStart()
    controlLoop = CMISS.ControlLoop()
    problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE], controlLoop)
    #controlLoop.MaximumIterationsSet(2)
    problem.ControlLoopCreateFinish()
    # Output
    print "----> Set up control loop <----\n"
    # Create nonlinear numerical solver
    linearSolver = CMISS.Solver()
    nonLinearSolver = CMISS.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
    nonLinearSolver.OutputTypeSet(CMISS.SolverOutputTypes.PROGRESS)
    nonLinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.EQUATIONS)
    #nonLinearSolver.NewtonRelativeToleranceSet(TOL)
    nonLinearSolver.NewtonLinearSolverGet(linearSolver)
    linearSolver.LinearTypeSet(CMISS.LinearSolverTypes.DIRECT)
    problem.SolversCreateFinish()
    # Output
    print "----> Set up linear and nonlinear solvers <----\n"
    # Add solver equations sets which encompass the physics
    solverEquations = CMISS.SolverEquations()
    solver = CMISS.Solver()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE], 1, solver)
    solver.SolverEquationsGet(solverEquations)
    solverEquations.SparsityTypeSet(CMISS.SolverEquationsSparsityTypes.SPARSE)
    solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()
    # Output
    print "----> Set up solver with equations <----\n"
    return problem, solverEquations


#=================================================================================#

#=================================================================================#
def BCCubeSingleFace(solverEquations, dependentField, appliedFace, faceNormal, appliedDirection, increm, optionBC,
                     fixXFace, fixYFace, fixZFace, numNodes, option):
    # Set up
    boundaryConditions = CMISS.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    # Initialise fixed faces node values.
    for node in fixXFace:
        boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,
                                   CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 1,
                                   CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    for node in fixYFace:
        boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,
                                   CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 2,
                                   CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    for node in fixZFace:
        boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,
                                   CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, 3,
                                   CMISS.BoundaryConditionsTypes.FIXED, 0.0)

    if option[0] == 2:
        # Fix derivatives
        if faceNormal == 1:
            derivFix = [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                        CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
        elif faceNormal == 2:
            derivFix = [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                        CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
        else:
            derivFix = [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                        CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2]

        for i in range(1, numNodes + 1):
            for j in derivFix:
                boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1, j, i, appliedDirection,
                                           CMISS.BoundaryConditionsTypes.FIXED, 0.0)

        # Fix all second and third derivatives.
        for i in range(1, numNodes + 1):
            for j in [CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                      CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                      CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,
                      CMISS.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]:
                for k in [1, 2, 3]:
                    boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1, j, i, k,
                                               CMISS.BoundaryConditionsTypes.FIXED, 0.0)
    # Output
    print "----> Implemented fixed boundary conditions <----\n"

    # Initialise applied faces.
    if optionBC == 1:
        # Option 1: Compression/extension
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.U, 1,
                                       CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, appliedDirection,
                                       CMISS.BoundaryConditionsTypes.FIXED, increm)
        # Output
        print "----> Implemented compression/extension boundary condition of " + str(increm) + " <----\n"

    elif optionBC == 2:
        # Option 2: Force
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,
                                       CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, appliedDirection,
                                       CMISS.BoundaryConditionsTypes.FIXED_INCREMENTED, increm)
        # Output
        print "----> Implemented force boundary condition of " + str(increm) + "N <----\n"
    elif optionBC == 3:
        # Option 3: Pressure
        for node in appliedFace:
            boundaryConditions.AddNode(dependentField, CMISS.FieldVariableTypes.DELUDELN, 1,
                                       CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV, node, faceNormal,
                                       CMISS.BoundaryConditionsTypes.PRESSURE_INCREMENTED, increm)
        # Output
        print "----> Implemented pressure boundary condition of " + str(increm) + " kPa <----\n"

    solverEquations.BoundaryConditionsCreateFinish()


#=================================================================================#

#=================================================================================#
def ExportResults(dependentField, deformedFieldUserNumber, decomposition, region, filename, option):
    deformedField = CMISS.Field()
    deformedField.CreateStart(deformedFieldUserNumber, region)
    deformedField.MeshDecompositionSet(decomposition)
    deformedField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
    deformedField.VariableLabelSet(CMISS.FieldVariableTypes.U, "DeformedGeometry")

    if option[0] == 1:
        # Trilinear.
        for component in [1, 2, 3, ]:
            deformedField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, component, 1)

    elif option[0] == 2:
        # Tricubic hermite. Geometry interpolated using cubic hermite basis (2nd mesh component).
        for component in [1, 2, 3, ]:
            deformedField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U, component, 1)
        if option[1] == 1:
            deformedField.ScalingTypeSet(CMISS.FieldScalingTypes.UNIT)
        elif option[1] == 2:
            deformedField.ScalingTypeSet(CMISS.FieldScalingTypes.ARITHMETIC_MEAN)

    deformedField.CreateFinish()
    for component in [1, 2, 3]:
        dependentField.ParametersToFieldParametersComponentCopy(CMISS.FieldVariableTypes.U,
                                                                CMISS.FieldParameterSetTypes.VALUES, component,
                                                                deformedField, CMISS.FieldVariableTypes.U,
                                                                CMISS.FieldParameterSetTypes.VALUES, component)


    # Export deformation.
    fields = CMISS.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(filename, "FORTRAN")
    fields.ElementsExport(filename, "FORTRAN")
    fields.Finalise()

    # Output
    print "----> Export solutions <----\n"

#=================================================================================#
