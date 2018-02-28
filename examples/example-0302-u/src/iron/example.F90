!> \file
!> \author Andreas Hessenthaler
!> \brief This is an example program to solve a lid-driven cavity problem using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> Main program
PROGRAM LidDrivenCavity

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: VelocityBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: VelocityMeshComponent=1
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponent=2
  REAL(CMISSRP),      PARAMETER :: PI=4.0_CMISSRP*DATAN(1.0_CMISSRP)

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg)    :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=255)    :: COMMAND_ARGUMENT,Filename
  INTEGER(CMISSIntg)    :: NUMBER_OF_GAUSS_XI
  INTEGER(CMISSIntg)    :: MeshRefinementLevel,NumberOfTimeSteps,TimeStep
  REAL(CMISSRP)         :: StartTime,StopTime,TimeStepSize,LoopStartTime,LoopStopTime
  REAL(CMISSRP)         :: Density,Viscosity,VelocityBC
  INTEGER(CMISSIntg)    :: SolverType,OutputFrequency

  !CMISS variables

  TYPE(cmfe_BasisType) :: VelocityBasis,PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: VelocityElements
  TYPE(cmfe_MeshElementsType) :: PressureElements
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,NonlinearSolver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg)    :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg)    :: EquationsSetIndex
  INTEGER(CMISSIntg)    :: NumberOfDimensions,NumberOfNodes,NumberOfElements
  INTEGER(CMISSIntg)    :: NumberOfNodesPerElement, NumberOfBoundaryPatches
  INTEGER(CMISSIntg)    :: NodeNumber,NodeDomain,NumberOfPatchIDs
  INTEGER(CMISSIntg)    :: NodeIdx,ComponentIdx,ElementIdx,CurrentPatchID,PatchIdx,StartIdx,StopIdx
  INTEGER(CMISSIntg)    :: Err
  LOGICAL               :: FileExists,BoundaryFound
  REAL(CMISSRP)         :: vx,vy,vz,x,y,z

  REAL(CMISSRP),        ALLOCATABLE :: NodesImport(:,:)           !< The coordinates of the mesh nodes
  INTEGER(CMISSIntg),   ALLOCATABLE :: ElementsImport(:,:)        !< The node IDs for each element
  INTEGER(CMISSIntg),   ALLOCATABLE :: BoundaryPatchesImport(:) !< The boundary patch labels for all boundary nodes

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS == 8) THEN
    ! get NumberOfDimensions
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberOfDimensions
    IF(NumberOfDimensions<2) CALL HANDLE_ERROR("Invalid number of dimensions.")
    ! get MeshRefinementLevel
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) MeshRefinementLevel
    IF(MeshRefinementLevel<1) CALL HANDLE_ERROR("Invalid mesh refinement level.")
    ! get StartTime, StopTime, TimeStepSize
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) StartTime
    IF(StartTime<0.0_CMISSRP) CALL HANDLE_ERROR("Invalid starting time.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) StopTime
    IF(StopTime<0.0_CMISSRP) CALL HANDLE_ERROR("Invalid stopping time.")
    CALL GET_COMMAND_ARGUMENT(5,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 5.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) TimeStepSize
    IF(TimeStepSize<=0.0_CMISSRP) CALL HANDLE_ERROR("Invalid time step size.")
    ! Material parameters for incompressible Newtonian fluid model
    ! get Viscosity, Density
    CALL GET_COMMAND_ARGUMENT(6,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 6.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) Viscosity
    IF(Viscosity<0.0_CMISSRP) CALL HANDLE_ERROR("Invalid viscosity.")
    CALL GET_COMMAND_ARGUMENT(7,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 7.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) Density
    IF(Density<0.0_CMISSRP) CALL HANDLE_ERROR("Invalid density.")
    ! get SolverType (direct, iterative)
    CALL GET_COMMAND_ARGUMENT(8,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 8.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SolverType
    IF((SolverType<0).OR.(SolverType>1)) CALL HANDLE_ERROR("Invalid solver type specification.")
  ELSE
    !If there are not enough arguments default the problem specification 
    NumberOfDimensions          = 2
    MeshRefinementLevel         = 1
    StartTime                   = 0.0_CMISSRP
    StopTime                    = 1.0_CMISSRP
    TimeStepSize                = 0.001_CMISSRP
    Density                     = 1.0_CMISSRP
    Viscosity                   = 0.025_CMISSRP
    SolverType                  = 0
  ENDIF
  NumberOfTimeSteps=(StopTime-StartTime)/TimeStepSize
  OutputFrequency=INT(NumberOfTimeSteps/10)

  ! Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  CALL cmfe_RandomSeedsSet(9999,Err)

  ! Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  ! Create a rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfDimensions,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  ! Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  ! Define basis functions for velocity
  CALL cmfe_Basis_Initialise(VelocityBasis,Err)
  CALL cmfe_Basis_CreateStart(VelocityBasisUserNumber,VelocityBasis,Err)
  CALL cmfe_Basis_TypeSet(VelocityBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(VelocityBasis,NumberOfDimensions,Err)
  SELECT CASE(NumberOfDimensions)
  CASE(2)
    CALL cmfe_Basis_InterpolationXiSet(VelocityBasis, &
      & [CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      &  CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION],Err)
  CASE(3)
    CALL cmfe_Basis_InterpolationXiSet(VelocityBasis, &
      & [CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      &  CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
      &  CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION],Err)
  CASE DEFAULT
    ! Do nothing
  END SELECT
  CALL cmfe_Basis_CreateFinish(VelocityBasis,Err)

  ! Define basis for pressure
  CALL cmfe_Basis_Initialise(PressureBasis,Err)
  CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(PressureBasis,NumberOfDimensions,Err)
  SELECT CASE(NumberOfDimensions)
  CASE(2)
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis, &
      & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
      &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
  CASE(3)
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis, &
      & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
      &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
      &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
  CASE DEFAULT
    ! Do nothing
  END SELECT
  CALL cmfe_Basis_CreateFinish(PressureBasis,Err)

  ! get user-defined mesh file name
  SELECT CASE(MeshRefinementLevel)
  CASE(1)
    WRITE(Filename, "(A18,I1,A40)") "src/cheart/meshes/",NumberOfDimensions,"D_MeshRefinementLevel_001/domain_quad_FE"
  CASE(2)
    WRITE(Filename, "(A18,I1,A40)") "src/cheart/meshes/",NumberOfDimensions,"D_MeshRefinementLevel_002/domain_quad_FE"
  CASE(3)
    WRITE(Filename, "(A18,I1,A40)") "src/cheart/meshes/",NumberOfDimensions,"D_MeshRefinementLevel_003/domain_quad_FE"
  CASE DEFAULT
    WRITE(Filename, "(A)") "src/cheart/meshes/MeshRefinementLevel_001/domain_quad_FE"
  END SELECT
  ! Check whether file exists
  INQUIRE(FILE=trim(Filename)//".X",EXIST=FileExists)
  IF(.NOT.FileExists) THEN
    CALL HANDLE_ERROR("File does not exist: "//trim(Filename)//".X")
  ENDIF
  INQUIRE(FILE=trim(Filename)//".T",EXIST=FileExists)
  IF(.NOT.FileExists) THEN
    CALL HANDLE_ERROR("File does not exist: "//trim(Filename)//".T")
  ENDIF
  INQUIRE(FILE=trim(Filename)//".B",EXIST=FileExists)
  IF(.NOT.FileExists) THEN
    CALL HANDLE_ERROR("File does not exist: "//trim(Filename)//".B")
  ENDIF
  ! Read CHeart mesh based on the given command line arguments
  WRITE(*,*) "Reading CHeart mesh data file "//TRIM(Filename)
  CALL cmfe_ReadMesh(trim(Filename), NodesImport, ElementsImport, BoundaryPatchesImport, "CHeart", Err)
  NumberOfNodes             = SIZE(NodesImport,1)
  NumberOfDimensions        = SIZE(NodesImport,2)
  NumberOfElements          = SIZE(ElementsImport,1)
  NumberOfNodesPerElement   = SIZE(ElementsImport,2)
  NumberOfBoundaryPatches   = BoundaryPatchesImport(1)
  WRITE(*,*) "...done"

  !== Set up mesh from imported mesh data
  ! Set up nodes
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  ! Set up mesh with velocity and pressure components
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfDimensions,Mesh,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)
  ! Set up mesh elements for velocity and pressure
  CALL cmfe_MeshElements_Initialise(VelocityElements,Err)
  CALL cmfe_MeshElements_Initialise(PressureElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,VelocityMeshComponent,VelocityBasis,VelocityElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,PressureMeshComponent,PressureBasis,PressureElements,Err)
  ! Set element connectivity
  DO ElementIdx=1,NumberOfElements
    CALL cmfe_MeshElements_NodesSet(VelocityElements,ElementIdx,ElementsImport(ElementIdx,:),Err)
    CALL cmfe_MeshElements_NodesSet(PressureElements,ElementIdx,ElementsImport(ElementIdx,1:NumberOfDimensions+1),Err)
  END DO
  ! Finish mesh elements
  CALL cmfe_MeshElements_CreateFinish(VelocityElements,Err)
  CALL cmfe_MeshElements_CreateFinish(PressureElements,Err)
  ! Finish mesh
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  ! Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  ! Create a geometric field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions,Err)  
  DO NodeIdx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,NodeIdx,VelocityMeshComponent,Err)
  END DO
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !== Update the geometric field parameters from the imported node coordinates
  ! For all node IDs
  DO NodeIdx=1,NumberOfNodes
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeIdx,1,NodeDomain,Err)
    ! Check if node is present in this processors' computational domain
    IF(NodeDomain==ComputationalNodeNumber) THEN
      ! For all components
      DO ComponentIdx=1,NumberOfDimensions
        ! Update coordinates
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,CMFE_NO_GLOBAL_DERIV,NodeIdx,ComponentIdx, &
          & NodesImport(NodeIdx,ComponentIdx),Err)
      END DO
    END IF
  END DO

  ! Create the equations set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    &  CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    &  CMFE_EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  ! Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"VectorField",Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"VectorField (derivative)",Err)
  DO ComponentIdx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,ComponentIdx,VelocityMeshComponent,Err)
    CALL cmfe_Field_ComponentMeshComponentSet( &
      & DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,ComponentIdx,VelocityMeshComponent,Err)
  END DO
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)
  DO ComponentIdx=1,NumberOfDimensions+1
    CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx,0.0_CMISSRP,Err)
  END DO
  CALL cmfe_Field_ParameterSetUpdateStart(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  ! Auto-create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialFieldUserNumber,MaterialField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  ! Set viscosity and density
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,Viscosity,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,Density,Err)

  ! Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  ! Start the creation of a problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber, &
    & [CMFE_PROBLEM_FLUID_MECHANICS_CLASS, &
    &  CMFE_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    &  CMFE_PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE], &
    & Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  ! Start the creation of the problem control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,StartTime,StopTime,TimeStepSize,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
 
  ! Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_DynamicThetaSet(Solver,0.5_CMISSRP,Err)
  CALL cmfe_Solver_DynamicNonlinearSolverGet(Solver,NonlinearSolver,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolver, &
    & CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonLineSearchTypeSet(NonlinearSolver,CMFE_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  CALL cmfe_Solver_NewtonLineSearchAlphaSet(NonlinearSolver,0.5_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLineSearchMonitorOutputSet(NonlinearSolver,.TRUE.,Err)  
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(NonlinearSolver,1.0E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(NonlinearSolver,1.0E-6_CMISSRP,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(NonlinearSolver,50,Err)
  CALL cmfe_Solver_NewtonMaximumFunctionEvaluationsSet(NonlinearSolver,50,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  CALL cmfe_Solver_OutputTypeSet(LinearSolver,CMFE_SOLVER_NO_OUTPUT,Err)
  IF(SolverType==0) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolver,1.0E-12_CMISSRP,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolver,1.0E-12_CMISSRP,Err)
  ENDIF
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  ! Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  ! Prescribe boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  ! BC correspond to lid-driven cavity in 2D or 3D
  ! in both cases, the y=1 boundary is the lid with shear flow, no-slip otherwise
  IF (NumberOfDimensions==2) THEN
      !=== Inflow boundary
      ! Get index in boundary file
      CurrentPatchID=4
      CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
      ! Now, set boundary condition
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_Field_ParameterSetGetNode( &
            & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
            & 1,1,NodeNumber,1,x,Err)
          vy = 0.0_CMISSRP
          vx = 1.0_CMISSRP*DSIN(PI*x)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_FIXED,vx,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 2,CMFE_BOUNDARY_CONDITION_FIXED,vy,Err)
        END IF
      END DO
      !=== No-slip boundary
      DO PatchIdx=1,4
        IF (PatchIdx==4) CYCLE
        ! Get index in boundary file
        CurrentPatchID=PatchIdx
        CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
        ! Now, set boundary condition
        DO NodeIdx=StartIdx,StopIdx
          NodeNumber=BoundaryPatchesImport(NodeIdx)
          CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
          IF(NodeDomain==ComputationalNodeNumber) THEN
            vy = 0.0_CMISSRP
            vx = 0.0_CMISSRP
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 1,CMFE_BOUNDARY_CONDITION_FIXED,vx,Err)
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 2,CMFE_BOUNDARY_CONDITION_FIXED,vy,Err)
          END IF
        END DO
      END DO
  ELSE
    !=== Inflow boundary
    ! Get index in boundary file
    CurrentPatchID=4
    CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
    ! Now, set boundary condition
    DO NodeIdx=StartIdx,StopIdx
      NodeNumber=BoundaryPatchesImport(NodeIdx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetGetNode( &
          & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
          & 1,1,NodeNumber,1,x,Err)
        vx = 1.0_CMISSRP*DSIN(PI*x)
        vy = 0.0_CMISSRP
        vz = 0.0_CMISSRP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & 1,CMFE_BOUNDARY_CONDITION_FIXED,vx,Err)
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & 2,CMFE_BOUNDARY_CONDITION_FIXED,vy,Err)
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & 3,CMFE_BOUNDARY_CONDITION_FIXED,vz,Err)
      END IF
    END DO
    !=== No-slip boundary
    DO PatchIdx=1,6
      IF (PatchIdx==4) CYCLE
      ! Get index in boundary file
      CurrentPatchID=PatchIdx
      CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
      ! Now, set boundary condition
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          vx = 0.0_CMISSRP
          vy = 0.0_CMISSRP
          vz = 0.0_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_FIXED,vx,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 2,CMFE_BOUNDARY_CONDITION_FIXED,vy,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 3,CMFE_BOUNDARY_CONDITION_FIXED,vz,Err)
        END IF
      END DO
    END DO
  END IF
  ! Finish BC solver equations
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  ! Solve multiple timesteps
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  WRITE(Filename, "(A20,I1,A24,I1,A8)") "results/current_run/",NumberOfDimensions, &
    & "D_MeshRefinementLevel_00",MeshRefinementLevel,"/Example"
  CALL cmfe_Fields_ElementsExport(Fields,Filename,"FORTRAN",Err)
  CALL cmfe_Fields_NodesExport(Fields,Filename,"FORTRAN",Err)
  DO TimeStep=1,NumberOfTimeSteps
    LoopStartTime=StartTime+FLOAT(TimeStep-1_CMISSIntg)*TimeStepSize
    LoopStopTime=StartTime+FLOAT(TimeStep)*TimeStepSize
    WRITE(*,*) "TIME STEP ",TimeStep,"/",NumberOfTimeSteps
    WRITE(*,*) " from ",LoopStartTime," to ",LoopStopTime," with dt ",TimeStepSize
    CALL cmfe_ControlLoop_TimesSet(ControlLoop,LoopStartTime,LoopStopTime,TimeStepSize,Err)
    ! Update BC - same as above, we are just doing this for demonstrating
    !             how to modify BC values in time
    IF (NumberOfDimensions==2) THEN
      !=== Inflow boundary
      ! Get index in boundary file
      CurrentPatchID=4
      CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
      ! Now, set boundary condition
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_Field_ParameterSetGetNode( &
            & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
            & 1,1,NodeNumber,1,x,Err)
          vy = 0.0_CMISSRP
          vx = 1.0_CMISSRP*DSIN(PI*x)
          CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & 1,1,NodeNumber,1,vx,Err)
          CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & 1,1,NodeNumber,2,vy,Err)
        END IF
      END DO
      !=== No-slip boundary
      DO PatchIdx=1,4
        IF (PatchIdx==4) CYCLE
        ! Get index in boundary file
        CurrentPatchID=PatchIdx
        CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
        ! Now, set boundary condition
        DO NodeIdx=StartIdx,StopIdx
          NodeNumber=BoundaryPatchesImport(NodeIdx)
          CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
          IF(NodeDomain==ComputationalNodeNumber) THEN
            vy = 0.0_CMISSRP
            vx = 0.0_CMISSRP
            CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
              & 1,1,NodeNumber,1,vx,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
              & 1,1,NodeNumber,2,vy,Err)
          END IF
        END DO
      END DO
    ELSE
      !=== Inflow boundary
      ! Get index in boundary file
      CurrentPatchID=4
      CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
      ! Now, set boundary condition
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_Field_ParameterSetGetNode( &
            & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
            & 1,1,NodeNumber,1,x,Err)
          vx = 1.0_CMISSRP*DSIN(PI*x)
          vy = 0.0_CMISSRP
          vz = 0.0_CMISSRP
          CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & 1,1,NodeNumber,1,vx,Err)
          CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & 1,1,NodeNumber,2,vy,Err)
          CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
            & 1,1,NodeNumber,3,vz,Err)
        END IF
      END DO
      !=== No-slip boundary
      DO PatchIdx=1,6
        IF (PatchIdx==4) CYCLE
        ! Get index in boundary file
        CurrentPatchID=PatchIdx
        CALL cmfe_ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
        ! Now, set boundary condition
        DO NodeIdx=StartIdx,StopIdx
          NodeNumber=BoundaryPatchesImport(NodeIdx)
          CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
          IF(NodeDomain==ComputationalNodeNumber) THEN
            vx = 0.0_CMISSRP
            vy = 0.0_CMISSRP
            vz = 0.0_CMISSRP
            CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
              & 1,1,NodeNumber,1,vx,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
              & 1,1,NodeNumber,2,vy,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
              & 1,1,NodeNumber,3,vz,Err)
          END IF
        END DO
      END DO
    END IF
    CALL cmfe_Field_ParameterSetUpdateStart(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    ! The actual solve
    CALL cmfe_Problem_Solve(Problem,Err)
    ! Export solution
    IF (MOD(TimeStep,OutputFrequency)/=0) CYCLE
    IF(TIMESTEP.LE.9) THEN
      WRITE(filename, "(A20,I1,A24,I1,A9,I1)") "results/current_run/",NumberOfDimensions, &
        & "D_MeshRefinementLevel_00",MeshRefinementLevel,"/Example_",TIMESTEP
      filename=trim(filename)
    ELSEIF(TIMESTEP.LE.99) THEN
      WRITE(filename, "(A20,I1,A24,I1,A9,I2)") "results/current_run/",NumberOfDimensions, &
        & "D_MeshRefinementLevel_00",MeshRefinementLevel,"/Example_",TIMESTEP
      filename=trim(filename)
    ELSEIF(TIMESTEP.LE.999) THEN
      WRITE(filename, "(A20,I1,A24,I1,A9,I3)") "results/current_run/",NumberOfDimensions, &
        & "D_MeshRefinementLevel_00",MeshRefinementLevel,"/Example_",TIMESTEP
      filename=trim(filename)
    ELSE
      WRITE(filename, "(A20,I1,A24,I1,A9,I4)") "results/current_run/",NumberOfDimensions, &
        & "D_MeshRefinementLevel_00",MeshRefinementLevel,"/Example_",TIMESTEP
      filename=trim(filename)
    END IF
    CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  END DO

  !== Clean house
  CALL cmfe_Fields_Finalise(Fields,Err)
  ! Deallocate variables that store mesh data as soon as we don't need them anymore
  IF(ALLOCATED(NodesImport))            DEALLOCATE(NodesImport)
  IF(ALLOCATED(ElementsImport))         DEALLOCATE(ElementsImport)
  IF(ALLOCATED(BoundaryPatchesImport))  DEALLOCATE(BoundaryPatchesImport)
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM LidDrivenCavity
