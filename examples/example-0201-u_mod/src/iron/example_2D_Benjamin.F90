!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a uniaxial extension equation using OpenCMISS calls.
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
!> Contributor(s): Andreas Hessenthaler, Thomas Klotz
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
PROGRAM FORTRANEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE meshReader
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
  REAL(CMISSRP) :: Width
  REAL(CMISSRP) :: Length
  REAL(CMISSRP) :: Height
  REAL(CMISSRP) :: MooneyRivlin1,MooneyRivlin2
  ! for Neo-Hookean solid, set MooneyRivlin2==0 and MooneyRivlin1 consistent
  ! with the 1st Piola-Kirchhoff stress tensor:
  ! P = 2.0 * MooneyRivlin1 / (det F)^(2/NumberOfDimensions) * (F - (F : F) / NumberOfDimensions * F^{-T}) - pressure * F^{-T}
  REAL(CMISSRP) :: BCDISP_MAX,NeumannBCvalue,lambda,areaChange,GaussXi,NodalForce
  INTEGER(CMISSIntg) :: DisplacementInterpolationType
  INTEGER(CMISSIntg) :: PressureInterpolationType
  INTEGER(CMISSIntg) :: PressureMeshComponent
  INTEGER(CMISSIntg) :: NumberOfGaussXi
  INTEGER(CMISSIntg) :: NumberOfDimensions
  INTEGER(CMISSIntg) :: ScalingType,GaussPoint
  INTEGER(CMISSIntg) :: SolverIsDirect,JACOBIAN_FD
  REAL(CMISSRP), PARAMETER :: Density=9.0E-4_CMISSRP !in g mm^-3
  REAL(CMISSRP), PARAMETER :: Gravity(3)=[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !in m s^-2
  INTEGER(CMISSIntg) :: NumberOfLoadIncrements=2
  INTEGER(CMISSIntg) :: useGeneratedMesh,bcDirichlet
  INTEGER(CMISSIntg), ALLOCATABLE :: nodeNumbers(:)
  REAL(CMISSRP),      ALLOCATABLE :: nodalWeights(:)

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DisplacementBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldSourceUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: deformedFieldUserNumber = 7
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program variables
  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,node_idx_2,component_idx,deriv_idx,NumberOfNodes
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:),RightSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi,RightNormalXi
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument,filename


  ! consistent nodal forces stuff
  INTEGER(CMISSIntg)  :: numNodesX,numNodesY,numNodesZ,numNodesXY,numNodesXZ,numNodesYZ
  REAL(CMISSRP)       :: weight

  ! CMISS variables
  TYPE(cmfe_BasisType)                  :: DisplacementBasis
  TYPE(cmfe_BasisType)                  :: PressureBasis
  TYPE(cmfe_BoundaryConditionsType)     :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType)       :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_DecompositionType)          :: Decomposition
  TYPE(cmfe_EquationsType)              :: Equations
  TYPE(cmfe_EquationsSetType)           :: EquationsSet
  TYPE(cmfe_FieldType)                  :: GeometricField
  TYPE(cmfe_FieldType)                  :: FibreField
  TYPE(cmfe_FieldType)                  :: MaterialField
  TYPE(cmfe_FieldType)                  :: DependentField
  TYPE(cmfe_FieldType)                  :: SourceField
  TYPE(cmfe_FieldType)                  :: deformedField
  TYPE(cmfe_FieldType)                  :: EquationsSetField
  TYPE(cmfe_FieldsType)                 :: Fields
  TYPE(cmfe_MeshElementsType)           :: DisplacementElements,PressureElements
  TYPE(cmfe_GeneratedMeshType)          :: GeneratedMesh
  TYPE(cmfe_MeshType)                   :: Mesh
  TYPE(cmfe_NodesType)                  :: Nodes
  TYPE(cmfe_ProblemType)                :: Problem
  TYPE(cmfe_RegionType)                 :: Region,WorldRegion
  TYPE(cmfe_SolverType)                 :: Solver
  TYPE(cmfe_SolverType)                 :: LinearSolver
  TYPE(cmfe_SolverEquationsType)        :: SolverEquations
  TYPE(cmfe_ControlLoopType)            :: ControlLoop

  ! generic CMISS variables
  INTEGER(CMISSIntg)                    :: Err

  ! local variables
  LOGICAL                               :: directory_exists = .FALSE.
  LOGICAL                               :: FileExists

  ! local variables w.r.t. mesh import
  REAL(CMISSRP),        ALLOCATABLE     :: NodesImport(:,:)           !< The coordinates of the mesh nodes
  INTEGER(CMISSIntg),   ALLOCATABLE     :: ElementsImport(:,:)        !< The node IDs for each element
  INTEGER(CMISSIntg),   ALLOCATABLE     :: BoundaryPatchesImport(:) !< The boundary patch labels for all boundary nodes
  INTEGER(CMISSIntg)                    :: NumberOfElements
  INTEGER(CMISSIntg)                    :: NumberOfNodesPerElement
  INTEGER(CMISSIntg)                    :: NumberOfBoundaryPatches
  INTEGER(CMISSIntg)                    :: ElementIdx,NodeIdx,ComponentIdx,CurrentPatchID,StartIdx,StopIdx
  REAL(CMISSRP)                         :: x,y,z

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  CALL cmfe_OutputSetOn("Example",Err)

  DisplacementInterpolationType = CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  ScalingType                   = CMFE_FIELD_ARITHMETIC_MEAN_SCALING

!!!!!!!! Command Line Interface !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 13) THEN
    ! get extents of spatial domain
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) WIDTH
    IF(WIDTH<=0) CALL HANDLE_ERROR("Invalid width.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) HEIGHT
    IF(HEIGHT<=0) CALL HANDLE_ERROR("Invalid height.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) LENGTH
    IF(LENGTH<0) CALL HANDLE_ERROR("Invalid length.")
    ! number of elements in each coordinate direction
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberGlobalXElements
    IF(NumberGlobalXElements<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(5,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 5.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberGlobalYElements
    IF(NumberGlobalYElements<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(6,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 6.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NumberGlobalZElements
    IF(NumberGlobalZElements<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    ! solver type (0: direct linear solve,1: linear iterative solve)
    CALL GET_COMMAND_ARGUMENT(7,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 7.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SolverIsDirect
    IF((SolverIsDirect<0).OR.(SolverIsDirect>1)) CALL HANDLE_ERROR("Invalid solver type.")
    ! Jacobian type (0: analytical Jacobian, 1: finite difference Jacobian)
    CALL GET_COMMAND_ARGUMENT(8,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 8.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) JACOBIAN_FD
    IF((JACOBIAN_FD<0).OR.(JACOBIAN_FD>1)) CALL HANDLE_ERROR("Invalid Jacobian type.")
    ! Material Parameter -> Mooney-Rivlin parameter 1
    CALL GET_COMMAND_ARGUMENT(9,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 9.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) MooneyRivlin1
    IF(MooneyRivlin1<=0) CALL HANDLE_ERROR("Invalid Mooney-Rivlin specification.")
    ! Material Parameter -> Mooney-Rivlin parameter 2
    CALL GET_COMMAND_ARGUMENT(10,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 10.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) MooneyRivlin2
    IF(MooneyRivlin2<0) CALL HANDLE_ERROR("Invalid MooneyRivlin specification.")
    ! Is mesh user-defined?
    CALL GET_COMMAND_ARGUMENT(11,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 11.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) useGeneratedMesh
    IF((useGeneratedMesh<0).OR.(useGeneratedMesh>1)) CALL HANDLE_ERROR("Invalid mesh type.")
    ! BC -> Maximum Displacement in percent
    CALL GET_COMMAND_ARGUMENT(12,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 12.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) BCDISP_MAX
    IF(BCDISP_MAX<=-1.0) CALL HANDLE_ERROR("Invalid BC specification.")
    ! Do boundary condition as Dirichlet?
    CALL GET_COMMAND_ARGUMENT(13,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 13.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) bcDirichlet
    IF((bcDirichlet<0).OR.(bcDirichlet>1)) CALL HANDLE_ERROR("Invalid BC type.")
  ELSE
    ! defaults for input arguments
    WIDTH                         = 1.0_CMISSRP ! x
    HEIGHT                        = 1.0_CMISSRP ! y
    LENGTH                        = 0.0_CMISSRP ! z
    NumberGlobalXElements         = 1
    NumberGlobalYElements         = 1
    NumberGlobalZElements         = 0
    SolverIsDirect                = 1 ! direct solver by default
    JACOBIAN_FD                   = 1 ! finite-difference Jacobian by default
    MooneyRivlin1                 = 1.0_CMISSRP !0.5_CMISSRP*35.7_CMISSRP ! see note above for Neo-Hookean solid
    MooneyRivlin2                 = 0.0_CMISSRP  ! If MooneyRivlin2 == 0 --> Neo-Hookean solid
    useGeneratedMesh              = 1 ! generated mesh by default
    BCDISP_MAX                    = 0.2_CMISSRP
    bcDirichlet                   = 0 ! Neumann/Dirichlet BC by default
  ENDIF
  IF(NumberGlobalZElements >= 1) THEN
    NumberOfDimensions  = 3
  ELSE
    NumberOfDimensions  = 2
  END IF
  NumberOfNodes       = (2*NumberGlobalXElements+1)*(2*NumberGlobalYElements+1)*(2*NumberGlobalZElements+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NumberOfGaussXi=3
  PressureMeshComponent=2
  PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  WRITE(*,'("Interpolation: ", i3)') DisplacementInterpolationType
  WRITE(*,'("Elements: ", 3 i3)') NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  WRITE(*,'("Scaling type: ", i3)') ScalingType

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfDimensions,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis function for displacement
  CALL cmfe_Basis_Initialise(DisplacementBasis,Err)
  CALL cmfe_Basis_CreateStart(DisplacementBasisUserNumber,DisplacementBasis,Err)
  SELECT CASE(DisplacementInterpolationType)
  CASE(2)
    WRITE(*,*) "Lagrange-Hermite TP type"
    CALL cmfe_Basis_TypeSet(DisplacementBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(8)
    WRITE(*,*) "Simplex type"
    CALL cmfe_Basis_TypeSet(DisplacementBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL cmfe_Basis_NumberOfXiSet(DisplacementBasis,NumberOfDimensions,Err)
  IF(NumberOfDimensions == 3) THEN
    CALL cmfe_Basis_InterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType,DisplacementInterpolationType, &
        & DisplacementInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType, &
        & DisplacementInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  END IF
  CALL cmfe_Basis_CreateFinish(DisplacementBasis,Err)

  !Basis for pressure
  CALL cmfe_Basis_Initialise(PressureBasis,Err)
  CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(PressureInterpolationType)
  CASE(1)
    WRITE(*,*) "Lagrange-Hermite TP type"
    CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7)
    WRITE(*,*) "Simplex type"
    CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL cmfe_Basis_NumberOfXiSet(PressureBasis,NumberOfDimensions,Err)
  IF(NumberOfDimensions == 3) THEN
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  END IF
  CALL cmfe_Basis_CreateFinish(PressureBasis,Err)

  CALL cmfe_Mesh_Initialise(Mesh,Err)
  IF (useGeneratedMesh==1) THEN
    ! Start the creation of a generated mesh in the region
    CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis,PressureBasis],Err)
    IF(NumberOfDimensions==2) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height,Length],Err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
        & NumberGlobalZElements],Err)
    ENDIF
    CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  ELSE
    ! get user-defined mesh file name
    WRITE(Filename, "(A26,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A3)") &
      & "src/cheart/meshes/domain_l", &
      & INT(WIDTH),"x",INT(HEIGHT),"x",INT(LENGTH), &
      & "_n", &
      & NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements, &
      & "_i",DisplacementInterpolationType,"_FE"
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
    CALL ReadMesh(trim(Filename), NodesImport, ElementsImport, BoundaryPatchesImport, "CHeart", Err)
    WRITE(*,*) NumberOfNodes
    NumberOfNodes             = SIZE(NodesImport,1)
    NumberOfDimensions        = SIZE(NodesImport,2)
    NumberOfElements          = SIZE(ElementsImport,1)
    NumberOfNodesPerElement   = SIZE(ElementsImport,2)
    NumberOfBoundaryPatches   = BoundaryPatchesImport(1)
    WRITE(*,*) NumberOfNodes,NumberOfDimensions,NumberOfElements,NumberOfNodesPerElement,NumberOfBoundaryPatches
    WRITE(*,*) "...done"

    ! create nodes
    CALL cmfe_Nodes_Initialise(Nodes,Err)
    CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
    CALL cmfe_Nodes_CreateFinish(Nodes,Err)
    !create mesh
    CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NumberOfDimensions,Mesh,Err)
    CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)
    CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,2,Err)
    ! Set up mesh elements for displacement and pressure
    CALL cmfe_MeshElements_Initialise(DisplacementElements,Err)
    CALL cmfe_MeshElements_Initialise(PressureElements,Err)
    CALL cmfe_MeshElements_CreateStart(Mesh,1,DisplacementBasis,DisplacementElements,Err)
    CALL cmfe_MeshElements_CreateStart(Mesh,PressureMeshComponent,PressureBasis,PressureElements,Err)
    ! Set element connectivity
    DO ElementIdx=1,NumberOfElements
      CALL cmfe_MeshElements_NodesSet(DisplacementElements,ElementIdx, &
        & ElementsImport(ElementIdx,:),Err)
      CALL cmfe_MeshElements_NodesSet(PressureElements,ElementIdx, &
        & ElementsImport(ElementIdx,(/1,3,7,9,19,21,25,27/)),Err)
    END DO
    ! Finish mesh elements
    CALL cmfe_MeshElements_CreateFinish(DisplacementElements,Err)
    CALL cmfe_MeshElements_CreateFinish(PressureElements,Err)
    ! Finish mesh
    CALL cmfe_Mesh_CreateFinish(Mesh,Err)
  END IF

  ! create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  ! define geometry
  IF (useGeneratedMesh==1) THEN
    ! create a field to put the geometry (default is geometry)
    CALL cmfe_Field_Initialise(GeometricField,Err)
    CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
    CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
    CALL cmfe_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
    CALL cmfe_Field_CreateFinish(GeometricField,Err)
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  ELSE
    ! create a field to put the geometry (default is geometry)
    CALL cmfe_Field_Initialise(GeometricField,Err)
    CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
    CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
    CALL cmfe_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
    CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions,Err)
    DO NodeIdx=1,NumberOfDimensions
      CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,NodeIdx,1,Err)
    END DO
    CALL cmfe_Field_CreateFinish(GeometricField,Err)
    !== update the geometric field parameters from the imported node coordinates
    ! for all node IDs
    DO NodeIdx=1,NumberOfNodes
      WRITE(*,*) NodeIdx
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeIdx,1,NodeDomain,Err)
      ! check if node is present in this processors' computational domain
      IF(NodeDomain==ComputationalNodeNumber) THEN
        ! for all components
        DO ComponentIdx=1,NumberOfDimensions
          ! update coordinates
          CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
            & CMFE_FIELD_VALUES_SET_TYPE,1,CMFE_NO_GLOBAL_DERIV,NodeIdx,ComponentIdx, &
            & NodesImport(NodeIdx,ComponentIdx),Err)
        END DO
      END IF
    END DO
    WRITE(*,*) "User-defined geometric field done."
  END IF

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_ScalingTypeSet(FibreField,ScalingType,Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"Derivative",Err)
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ScalingTypeSet(DependentField,ScalingType,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,MooneyRivlin1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,MooneyRivlin2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,Density,Err)

  ! NO source 
  !Create the source field with the gravity vector
  !CALL cmfe_Field_Initialise(SourceField,Err)
  !CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,FieldSourceUserNumber,SourceField,Err)
  !CALL cmfe_Field_ScalingTypeSet(SourceField,ScalingType,Err)
  !CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)
  !DO component_idx=1,NumberOfDimensions
  !  CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
  !    & component_idx,Gravity(component_idx),Err)
  !ENDDO

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
!  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err) 
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE, &
      & component_idx,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,component_idx,Err)
  END DO
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,-8.0_CMISSRP,Err) !0.0_CMISSRP
  CALL cmfe_Field_ParameterSetUpdateStart(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  ! Create a deformed geometry field, as cmgui doesn't like displaying
  ! deformed fibres from the dependent field because it isn't a geometric field.
  CALL cmfe_Field_Initialise(deformedField,Err)
  CALL cmfe_Field_CreateStart(deformedFieldUserNumber,region,deformedField,Err)
  CALL cmfe_Field_MeshDecompositionSet(deformedField,decomposition,Err)
  CALL cmfe_Field_TypeSet(deformedField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_VariableLabelSet(deformedField,CMFE_FIELD_U_VARIABLE_TYPE,"DeformedGeometry",Err)
  DO componentIdx=1,numberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(deformedField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,Err)
  END DO
  !IF (interpolationType==4) CALL cmfe_Field_ScalingTypeSet(deformedField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  CALL cmfe_Field_CreateFinish(deformedField,Err)


  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber, &
    & [CMFE_PROBLEM_ELASTICITY_CLASS, &
    &  CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    &  CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,NumberOfLoadIncrements,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  IF(JACOBIAN_FD==1) THEN
    CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  ELSE
    CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  END IF
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(Solver,1.0E-12_CMISSRP,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  ! Chose Solver Type
  IF(SolverIsDirect==1) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err) !<Direct linear solver type.
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err) !<Iterative linear solver type.
  ENDIF
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  CALL cmfe_BoundaryConditions_NeumannSparsityTypeSet(BoundaryConditions,CMFE_BOUNDARY_CONDITION_SPARSE_MATRICES,Err)

  IF(useGeneratedMesh==1) THEN
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)

    ! Dirichlet BC at x=0
    DO node_idx=1,SIZE(LeftSurfaceNodes,1)
      NodeNumber=LeftSurfaceNodes(node_idx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        ! constrain x-direction
        CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & 1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      ENDIF
    ENDDO
    ! 2D: BC at x=0, y=0
    CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1, &
      & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)

    ! Dirichlet BC at x=0, y=ly/2, z=z/2
    !NodeNumber=(NumberOfNodes/2)+1-((2*NumberGlobalXElements+1)/2) ! implicit FLOOR() !!!
    !CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    !IF(NodeDomain==ComputationalNodeNumber) THEN
    !  DO component_idx=2,NumberOfDimensions
    !    CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
    !      & component_idx,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    !  END DO
    !END IF

    ! stretch
    lambda  = (1.0_CMISSRP+BCDISP_MAX)
    IF(bcDirichlet==1) THEN
      WRITE(*,*) "Applying displacement BC.."
      ! Dirichlet BC at x=lx
      DO node_idx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          ! constrain x-direction
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_FIXED,WIDTH*lambda,Err)
        ENDIF
      ENDDO
    ELSE
      WRITE(*,*) "Applying traction BC.."
      ! corresponding traction value
      ! For Neumann-integrated (get weights, only works for 3D!!!):
      !NeumannBCvalue = 2.0_CMISSRP * MooneyRivlin1 * (lambda - 1.0_CMISSRP / lambda / lambda) * HEIGHT * LENGTH
      ! For Neumann-point:
      NeumannBCvalue =0.5/ HEIGHT
      !WRITE(*,*) "  getting nodal weights"
      ! compute consistent nodal weights (mesh_reader function!)
      !CALL GeneratedMesh_SurfaceWeightsGet(nodalWeights,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE, &
      !  & numberGlobalXelements,numberGlobalYelements,numberGlobalZelements,Err)
      WRITE(*,*) "  applying consistent nodal forces"
      ! apply consistent nodal forces based on nodal weights
      ! Neumann BC at x=lx
      DO node_idx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(node_idx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          !NodalForce=NeumannBCvalue*nodalWeights(node_idx)
          !CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField, &
          !  & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
          !  & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_INTEGRATED,NodalForce,Err)
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,nodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_POINT,NeumannBCvalue,Err)          
        ENDIF
      ENDDO
    END IF
  ELSE
    !=== Wall boundary
    ! Get index in boundary file
    CurrentPatchID=1
    CALL ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
    ! Now, set boundary condition
    DO NodeIdx=StartIdx,StopIdx
      NodeNumber=BoundaryPatchesImport(NodeIdx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetGetNode( &
          & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
          & 1,1,NodeNumber,1,x,Err)
        CALL cmfe_Field_ParameterSetGetNode( &
          & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
          & 1,1,NodeNumber,2,y,Err)
        CALL cmfe_Field_ParameterSetGetNode( &
          & GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
          & 1,1,NodeNumber,3,z,Err)
        CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & 1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        ! constrain center-node also in y- and z-direction
        IF((ABS(y-0.5*HEIGHT)<1.0e-12).AND.(ABS(z-0.5*LENGTH)<1.0e-12)) THEN
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        END IF
      END IF
    END DO
    WRITE(*,*) "Wall BC done."
    !=== displacement/traction boundary
    ! Get index in boundary file
    CurrentPatchID=2
    CALL ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
    
    ! Now, set boundary condition
    node_idx  = 1
    DO NodeIdx=StartIdx,StopIdx
      NodeNumber=BoundaryPatchesImport(NodeIdx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        ! stretch
        lambda  = (1.0_CMISSRP+BCDISP_MAX)
        ! applyt Dirichlet or Neumann BC
        IF(bcDirichlet==1) THEN
          WRITE(*,*) "Applying displacement BC.."
          ! Dirichlet BC at x=lx
          CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
          IF(NodeDomain==ComputationalNodeNumber) THEN
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
              & 1,CMFE_BOUNDARY_CONDITION_FIXED,WIDTH*lambda,Err)
          END IF
        ELSE
          WRITE(*,*) "Applying traction BC.."
          ! corresponding traction value
          NeumannBCvalue = 2.0_CMISSRP * MooneyRivlin1 * (lambda - 1.0_CMISSRP / lambda / lambda) * HEIGHT * LENGTH
          WRITE(*,*) "  getting nodal weights"
          ! compute consistent nodal weights
          CALL ImportedMesh_SurfaceWeightsGet(nodalWeights,NumberOfNodesPerElement,Err)
          WRITE(*,*) "  applying consistent nodal forces"
          ! apply consistent nodal forces based on nodal weights
          ! Neumann BC at x=lx
          NodalForce=NeumannBCvalue*nodalWeights(node_idx)/(NumberGlobalYElements*NumberGlobalZElements)
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_INTEGRATED,NodalForce,Err)
        END IF
      END IF
      node_idx  = node_idx + 1
      IF(node_idx>9) node_idx = 1
    END DO
  END IF

  ! finish BC
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  ! Copy deformed geometry into deformed field
  DO componentIdx=1,numberOfDimensions
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,componentIdx, &
      & deformedField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,componentIdx,Err)
  END DO

  ! convenience switch for creating reference results
  ! TODO remove before release
  IF(.FALSE.) THEN
    ! Set export file name
    WRITE(filename, "(A24,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A2,I1,A3,I1,A3,I1,A3,I1)") &
      & "results/reference/iron/l", &
      & INT(WIDTH),"x",INT(HEIGHT),"x",INT(LENGTH), &
      & "_n", &
      & NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements, &
      & "_i",displacementInterpolationType,"_s",SolverIsDirect, &
      & "_fd",JACOBIAN_FD,"_gm",useGeneratedMesh,"_bc",bcDirichlet
    ! make sure directories exist
    INQUIRE(file="./results/", exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir ./results/")
    END IF
    INQUIRE(file="./results/reference/", exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir ./results/reference/")
    END IF
    INQUIRE(file="./results/reference/iron/", exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir ./results/reference/iron/")
    END IF
    INQUIRE(file=trim(filename), exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir "//trim(filename))
    END IF
  ELSE
    ! Set export file name
    WRITE(filename, "(A21,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A2,I1,A3,I1,A3,I1,A3,I1)") &
      & "results/current_run/l", &
      & INT(WIDTH),"x",INT(HEIGHT),"x",INT(LENGTH), &
      & "_n", &
      & NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements, &
      & "_i",displacementInterpolationType,"_s",SolverIsDirect, &
      & "_fd",JACOBIAN_FD,"_gm",useGeneratedMesh,"_bc",bcDirichlet
    ! make sure directories exist
    INQUIRE(file="./results/", exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir ./results/")
    END IF
    INQUIRE(file="./results/current_run/", exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir ./results/current_run/")
    END IF
    INQUIRE(file=trim(filename), exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir "//trim(filename))
    END IF
  END IF

  ! Export solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  filename=trim(filename)//trim("/Example")
  CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  IF(ALLOCATED(nodalWeights))           DEALLOCATE(nodalWeights)
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

END PROGRAM FORTRANEXAMPLE

