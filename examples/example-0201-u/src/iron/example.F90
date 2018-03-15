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


  ! test program parameters
  REAL(CMISSRP)                         :: Width
  REAL(CMISSRP)                         :: Length
  REAL(CMISSRP)                         :: Height
  REAL(CMISSRP)                         :: MooneyRivlin1
  REAL(CMISSRP)                         :: MooneyRivlin2
  ! for Neo-Hookean solid, set MooneyRivlin2==0 and MooneyRivlin1 consistent
  ! with the 1st Piola-Kirchhoff stress tensor:
  ! P = 2.0 * MooneyRivlin1 / (det F)^(2/NumberOfDimensions) * (F - (F : F) / NumberOfDimensions * F^{-T}) + pressure * F^{-T}
  REAL(CMISSRP)                         :: BCDISP_MAX
  REAL(CMISSRP)                         :: NeumannBCvalue
  REAL(CMISSRP)                         :: lambda
  REAL(CMISSRP)                         :: NodalForce
  INTEGER(CMISSIntg)                    :: DisplacementInterpolationType
  INTEGER(CMISSIntg)                    :: PressureInterpolationType
  INTEGER(CMISSIntg)                    :: PressureMeshComponent
  INTEGER(CMISSIntg)                    :: NumberOfGaussXi
  INTEGER(CMISSIntg)                    :: NumberOfDimensions
  INTEGER(CMISSIntg)                    :: ScalingType
  INTEGER(CMISSIntg)                    :: SolverIsDirect
  INTEGER(CMISSIntg)                    :: JACOBIAN_FD
  REAL(CMISSRP),      PARAMETER         :: Density                        = 0.0E-4_CMISSRP !in g mm^-3
  REAL(CMISSRP),      PARAMETER         :: Gravity(3)                     = [0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP] !in m s^-2
  INTEGER(CMISSIntg)                    :: NumberOfLoadIncrements
  INTEGER(CMISSIntg)                    :: useGeneratedMesh
  INTEGER(CMISSIntg)                    :: bcType
  INTEGER(CMISSIntg), ALLOCATABLE       :: nodeNumbers(:)
  REAL(CMISSRP),      ALLOCATABLE       :: nodalWeights(:)

  ! user numbers
  INTEGER(CMISSIntg), PARAMETER         :: CoordinateSystemUserNumber     = 1
  INTEGER(CMISSIntg), PARAMETER         :: RegionUserNumber               = 1
  INTEGER(CMISSIntg), PARAMETER         :: DisplacementBasisUserNumber    = 1
  INTEGER(CMISSIntg), PARAMETER         :: PressureBasisUserNumber        = 2
  INTEGER(CMISSIntg), PARAMETER         :: GeneratedMeshUserNumber        = 1
  INTEGER(CMISSIntg), PARAMETER         :: MeshUserNumber                 = 1
  INTEGER(CMISSIntg), PARAMETER         :: DecompositionUserNumber        = 1
  INTEGER(CMISSIntg), PARAMETER         :: FieldGeometryUserNumber        = 1
  INTEGER(CMISSIntg), PARAMETER         :: FieldFibreUserNumber           = 2
  INTEGER(CMISSIntg), PARAMETER         :: FieldMaterialUserNumber        = 3
  INTEGER(CMISSIntg), PARAMETER         :: FieldDependentUserNumber       = 4
  INTEGER(CMISSIntg), PARAMETER         :: FieldSourceUserNumber          = 5
  INTEGER(CMISSIntg), PARAMETER         :: EquationsSetFieldUserNumber    = 6
  INTEGER(CMISSIntg), PARAMETER         :: EquationSetUserNumber          = 1
  INTEGER(CMISSIntg), PARAMETER         :: ProblemUserNumber              = 1

  ! program variables
  INTEGER(CMISSIntg)                    :: NumberGlobalXElements
  INTEGER(CMISSIntg)                    :: NumberGlobalYElements
  INTEGER(CMISSIntg)                    :: NumberGlobalZElements
  INTEGER(CMISSIntg)                    :: EquationsSetIndex
  INTEGER(CMISSIntg)                    :: NumberOfComputationalNodes
  INTEGER(CMISSIntg)                    :: NumberOfDomains
  INTEGER(CMISSIntg)                    :: ComputationalNodeNumber
  INTEGER(CMISSIntg)                    :: NodeNumber
  INTEGER(CMISSIntg)                    :: NodeDomain
  INTEGER(CMISSIntg)                    :: NodeIdx
  INTEGER(CMISSIntg)                    :: ComponentIdx
  INTEGER(CMISSIntg)                    :: NumberOfNodes
  INTEGER(CMISSIntg), ALLOCATABLE       :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE       :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg)                    :: LeftNormalXi
  INTEGER(CMISSIntg)                    :: RightNormalXi
  INTEGER(CMISSIntg)                    :: NumberOfArguments
  INTEGER(CMISSIntg)                    :: ArgumentLength
  INTEGER(CMISSIntg)                    :: ArgumentStatus
  CHARACTER(LEN=255)                    :: CommandArgument
  CHARACTER(LEN=255)                    :: Filename

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
  INTEGER(CMISSIntg)                    :: ElementIdx
  INTEGER(CMISSIntg)                    :: CurrentPatchID
  INTEGER(CMISSIntg)                    :: StartIdx
  INTEGER(CMISSIntg)                    :: StopIdx
  REAL(CMISSRP)                         :: x
  REAL(CMISSRP)                         :: y
  REAL(CMISSRP)                         :: z

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
  PressureInterpolationType     = CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  PressureMeshComponent         = 2
  ScalingType                   = CMFE_FIELD_ARITHMETIC_MEAN_SCALING
  NumberOfGaussXi               = 3

  ! CLI
  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 14) THEN
    ! get extents of spatial domain
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) Width
    IF(Width<=0) CALL HandleError("Invalid width.")
    CALL GET_COMMAND_ARGUMENT(2,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 2.")
    READ(CommandArgument(1:ArgumentLength),*) Height
    IF(Height<=0) CALL HandleError("Invalid height.")
    CALL GET_COMMAND_ARGUMENT(3,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 3.")
    READ(CommandArgument(1:ArgumentLength),*) Length
    IF(Length<0) CALL HandleError("Invalid length.")
    ! number of elements in each coordinate direction
    CALL GET_COMMAND_ARGUMENT(4,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 4.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalXElements
    IF(NumberGlobalXElements<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(5,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 5.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalYElements
    IF(NumberGlobalYElements<=0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(6,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 6.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalZElements
    IF(NumberGlobalZElements<0) CALL HandleError("Invalid number of Z elements.")
    ! solver type (0: direct linear solve,1: linear iterative solve)
    CALL GET_COMMAND_ARGUMENT(7,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 7.")
    READ(CommandArgument(1:ArgumentLength),*) SolverIsDirect
    IF((SolverIsDirect<0).OR.(SolverIsDirect>1)) CALL HandleError("Invalid solver type.")
    ! Jacobian type (0: analytical Jacobian, 1: finite difference Jacobian)
    CALL GET_COMMAND_ARGUMENT(8,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 8.")
    READ(CommandArgument(1:ArgumentLength),*) JACOBIAN_FD
    IF((JACOBIAN_FD<0).OR.(JACOBIAN_FD>1)) CALL HandleError("Invalid Jacobian type.")
    ! Material Parameter -> Mooney-Rivlin parameter 1
    CALL GET_COMMAND_ARGUMENT(9,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 9.")
    READ(CommandArgument(1:ArgumentLength),*) MooneyRivlin1
    IF(MooneyRivlin1<=0) CALL HandleError("Invalid Mooney-Rivlin specification.")
    ! Material Parameter -> Mooney-Rivlin parameter 2
    CALL GET_COMMAND_ARGUMENT(10,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 10.")
    READ(CommandArgument(1:ArgumentLength),*) MooneyRivlin2
    IF(MooneyRivlin2<0) CALL HandleError("Invalid MooneyRivlin specification.")
    ! Is mesh user-defined?
    CALL GET_COMMAND_ARGUMENT(11,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 11.")
    READ(CommandArgument(1:ArgumentLength),*) useGeneratedMesh
    IF((useGeneratedMesh<0).OR.(useGeneratedMesh>1)) CALL HandleError("Invalid mesh type.")
    ! BC -> Maximum Displacement in percent
    CALL GET_COMMAND_ARGUMENT(12,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 12.")
    READ(CommandArgument(1:ArgumentLength),*) BCDISP_MAX
    IF(BCDISP_MAX<=-1.0) CALL HandleError("Invalid BC specification.")
    ! BC -> Maximum Displacement in percent
    CALL GET_COMMAND_ARGUMENT(13,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 13.")
    READ(CommandArgument(1:ArgumentLength),*) NumberOfLoadIncrements
    IF(NumberOfLoadIncrements<1) CALL HandleError("Invalid number of load increments.")
    ! Do boundary condition as Dirichlet?
    CALL GET_COMMAND_ARGUMENT(14,CommandArgument,ArgumentLength,ArgumentStatus)
    IF(ArgumentStatus>0) CALL HandleError("Error for command argument 14.")
    READ(CommandArgument(1:ArgumentLength),*) bcType
    IF((bcType<0).OR.(bcType>2)) CALL HandleError("Invalid BC type.")
    WRITE(*,*) "!=============================================================="
    WRITE(*,*) "User arguments:"
    WRITE(*,*) "  Width                 = ",Width
    WRITE(*,*) "  Height                = ",Height
    WRITE(*,*) "  Length                = ",Length
    WRITE(*,*) "  NumberGlobalXElements = ",NumberGlobalXElements
    WRITE(*,*) "  NumberGlobalYElements = ",NumberGlobalYElements
    WRITE(*,*) "  NumberGlobalZElements = ",NumberGlobalZElements
    WRITE(*,*) "  SolverIsDirect        = ",SolverIsDirect
    WRITE(*,*) "  JACOBIAN_FD           = ",JACOBIAN_FD
    WRITE(*,*) "  MooneyRivlin1         = ",MooneyRivlin1
    WRITE(*,*) "  MooneyRivlin2         = ",MooneyRivlin2
    WRITE(*,*) "  useGeneratedMesh      = ",useGeneratedMesh
    WRITE(*,*) "  BCDISP_MAX            = ",BCDISP_MAX
    WRITE(*,*) "  bcType                = ",bcType
    WRITE(*,*) "!=============================================================="
  ELSE
    ! defaults for input arguments
    Width                         = 1.0_CMISSRP
    Height                        = 0.2_CMISSRP
    Length                        = 0.2_CMISSRP
    NumberGlobalXElements         = 2
    NumberGlobalYElements         = 2
    NumberGlobalZElements         = 2
    SolverIsDirect                = 1 ! direct solver by default
    JACOBIAN_FD                   = 1 ! finite-difference Jacobian by default
    MooneyRivlin1                 = 0.5_CMISSRP*35.7_CMISSRP ! see note above for Neo-Hookean solid
    MooneyRivlin2                 = 0.0_CMISSRP  ! If MooneyRivlin2 == 0 --> Neo-Hookean solid
    useGeneratedMesh              = 1 ! generated mesh by default
    BCDISP_MAX                    = 0.2_CMISSRP
    NumberOfLoadIncrements        = 1 ! one load increment by default
    bcType                        = 0 ! 0 - Dirichlet BC by default; else: 1 - Neumann_integrated, 2 - Neumann_point
    WRITE(*,*) "!=============================================================="
    WRITE(*,*) "Default arguments:"
    WRITE(*,*) "  Width                 = ",Width
    WRITE(*,*) "  Height                = ",Height
    WRITE(*,*) "  Length                = ",Length
    WRITE(*,*) "  NumberGlobalXElements = ",NumberGlobalXElements
    WRITE(*,*) "  NumberGlobalYElements = ",NumberGlobalYElements
    WRITE(*,*) "  NumberGlobalZElements = ",NumberGlobalZElements
    WRITE(*,*) "  SolverIsDirect        = ",SolverIsDirect
    WRITE(*,*) "  JACOBIAN_FD           = ",JACOBIAN_FD
    WRITE(*,*) "  MooneyRivlin1         = ",MooneyRivlin1
    WRITE(*,*) "  MooneyRivlin2         = ",MooneyRivlin2
    WRITE(*,*) "  useGeneratedMesh      = ",useGeneratedMesh
    WRITE(*,*) "  BCDISP_MAX            = ",BCDISP_MAX
    WRITE(*,*) "  bcType                = ",bcType
    WRITE(*,*) "!=============================================================="
  END IF
  IF(NumberGlobalZElements >= 1) THEN
    NumberOfDimensions  = 3
  ELSE
    NumberOfDimensions  = 2
  END IF
  NumberOfNodes       = (2*NumberGlobalXElements+1)*(2*NumberGlobalYElements+1)*(2*NumberGlobalZElements+1)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains = NumberOfComputationalNodes

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
    END IF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType, &
        & DisplacementInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    END IF
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
    END IF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    END IF
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
    END IF
    CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  ELSE
    ! get user-defined mesh file name
    WRITE(Filename, "(A26,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A3)") &
      & "src/cheart/meshes/domain_l", &
      & INT(Width),"x",INT(Height),"x",INT(Length), &
      & "_n", &
      & NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements, &
      & "_i",DisplacementInterpolationType,"_FE"
    ! Check whether file exists
    INQUIRE(FILE=trim(Filename)//".X",EXIST=FileExists)
    IF(.NOT.FileExists) THEN
      CALL HandleError("File does not exist: "//trim(Filename)//".X")
    END IF
    INQUIRE(FILE=trim(Filename)//".T",EXIST=FileExists)
    IF(.NOT.FileExists) THEN
      CALL HandleError("File does not exist: "//trim(Filename)//".T")
    END IF
    INQUIRE(FILE=trim(Filename)//".B",EXIST=FileExists)
    IF(.NOT.FileExists) THEN
      CALL HandleError("File does not exist: "//trim(Filename)//".B")
    END IF
    ! Read CHeart mesh based on the given command line arguments
    WRITE(*,*) "Reading CHeart mesh data file "//TRIM(Filename)
    CALL ReadMesh(trim(Filename), NodesImport, ElementsImport, BoundaryPatchesImport, "CHeart", Err)
    NumberOfNodes             = SIZE(NodesImport,1)
    NumberOfDimensions        = SIZE(NodesImport,2)
    NumberOfElements          = SIZE(ElementsImport,1)
    NumberOfNodesPerElement   = SIZE(ElementsImport,2)
    NumberOfBoundaryPatches   = BoundaryPatchesImport(1)
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
    CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Undeformed",Err)
    CALL cmfe_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
    CALL cmfe_Field_CreateFinish(GeometricField,Err)
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  ELSE
    ! create a field to put the geometry (default is geometry)
    CALL cmfe_Field_Initialise(GeometricField,Err)
    CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
    CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Undeformed",Err)
    CALL cmfe_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
    CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions,Err)
    DO NodeIdx=1,NumberOfDimensions
      CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,NodeIdx,1,Err)
    END DO
    CALL cmfe_Field_CreateFinish(GeometricField,Err)
    !== update the geometric field parameters from the imported node coordinates
    ! for all node IDs
    DO NodeIdx=1,NumberOfNodes
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
  DO ComponentIdx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,ComponentIdx,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,ComponentIdx,1,Err)
  END DO
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

  !Create the source field with the gravity vector
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,FieldSourceUserNumber,SourceField,Err)
  CALL cmfe_Field_ScalingTypeSet(SourceField,ScalingType,Err)
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)
  DO ComponentIdx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & ComponentIdx,Gravity(ComponentIdx),Err)
  END DO

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  DO ComponentIdx=1,NumberOfDimensions
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE, &
      & ComponentIdx,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,ComponentIdx,Err)
  END DO
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

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
  CALL cmfe_ControlLoop_OutputTypeSet(CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
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
  CALL cmfe_Solver_NewtonLineSearchTypeSet(Solver,CMFE_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  CALL cmfe_Solver_NewtonLineSearchMonitorOutputSet(Solver,.TRUE.,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  ! Chose Solver Type
  IF(SolverIsDirect==1) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err) !<Direct linear solver type.
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err) !<Iterative linear solver type.
    !NO,JACOBI,BLOCK_JACOBI,SOR,INCOMPLETE_CHOLESKY,INCOMPLETE_LU,ADDITIVE_SCHWARZ
    CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(LinearSolver,CMFE_SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolver,1.0E-12_CMISSRP,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolver,1.0E-12_CMISSRP,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolver,1.0E+8_CMISSRP,Err)
  END IF
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

  ! generated mesh
  IF(useGeneratedMesh==1) THEN
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)

    ! Dirichlet BC at x=0
    DO NodeIdx=1,SIZE(LeftSurfaceNodes,1)
      NodeNumber=LeftSurfaceNodes(NodeIdx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        ! constrain x-direction
        CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & 1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      END IF
    END DO
    ! Dirichlet BC at x=0, y=ly/2, z=z/2
    NodeNumber=(NumberOfNodes/2)+1-((2*NumberGlobalXElements+1)/2) ! implicit FLOOR() !!!
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO ComponentIdx=2,NumberOfDimensions
        CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & ComponentIdx,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      END DO
    END IF
    ! Dirichlet BC at x=0, y=0, z=z/2 to avoid rotations
    NodeNumber=(2*NumberGlobalXElements+1)*NumberGlobalYElements+1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
        & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    END IF

    ! stretch
    lambda  = (1.0_CMISSRP+BCDISP_MAX)
    IF(bcType==0) THEN
      WRITE(*,*) "Applying displacement BC.."
      ! Dirichlet BC at x=lx
      DO NodeIdx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          ! constrain x-direction
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,Width*BCDISP_MAX,Err)
        END IF
      END DO
    ELSE IF(bcType==1) THEN
      WRITE(*,*) "Applying traction BC.."
      ! corresponding traction value
      NeumannBCvalue = 2.0_CMISSRP * MooneyRivlin1 * (lambda - 1.0_CMISSRP / lambda / lambda) * Height * Length
      WRITE(*,*) "  getting nodal weights"
      ! compute consistent nodal weights
      CALL GeneratedMesh_SurfaceWeightsGet(nodalWeights,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE, &
        & numberGlobalXelements,numberGlobalYelements,numberGlobalZelements,Err)
      WRITE(*,*) "  applying consistent nodal forces"
      ! apply consistent nodal forces based on nodal weights
      ! Neumann BC at x=lx
      DO NodeIdx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          NodalForce=NeumannBCvalue*nodalWeights(NodeIdx)
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_INTEGRATED,NodalForce,Err)
        END IF
      END DO
    ELSE IF(bcType==2) THEN
      WRITE(*,*) "Applying traction BC.."
      ! corresponding traction value
      NeumannBCvalue = 2.0_CMISSRP * MooneyRivlin1 * (lambda - 1.0_CMISSRP / lambda / lambda)
      WRITE(*,*) "  applying point forces"
      ! apply point forces, Neumann BC at x=lx
      DO NodeIdx=1,SIZE(RightSurfaceNodes,1)
        NodeNumber=RightSurfaceNodes(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_POINT,NeumannBCvalue,Err)
        END IF
      END DO
      CALL HandleError("Warning: BC type 2 not fully implemented for generated mesh.")
    END IF
  ! user-defined mesh
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
        IF((ABS(y-0.5*Height)<1.0e-12).AND.(ABS(z-0.5*Length)<1.0e-12)) THEN
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        END IF
        ! constrain bottom/center-node also in y-direction
        IF((ABS(y-0.5*Height)<1.0e-12).AND.(ABS(z-0.0*Length)<1.0e-12)) THEN
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        END IF
      END IF
    END DO
    WRITE(*,*) "Wall BC done."
    !=== displacement/traction boundary
    ! stretch
    lambda  = (1.0_CMISSRP+BCDISP_MAX)
    ! Get index in boundary file
    CurrentPatchID=2
    CALL ImportedMesh_SurfaceGet(BoundaryPatchesImport,CurrentPatchID,StartIdx,StopIdx,Err)
    ! apply Dirichlet BC
    IF(bcType==0) THEN
      WRITE(*,*) "Applying displacement BC.."
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          ! Dirichlet BC at x=lx
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,Width*lambda,Err)
        END IF
      END DO
    ! apply Neumann integrated BC
    ELSE IF(bcType==1) THEN
      WRITE(*,*) "Applying traction BC.."
      NodeIdx  = 1
      ! corresponding traction value
      NeumannBCvalue = 2.0_CMISSRP * MooneyRivlin1 * (lambda - 1.0_CMISSRP / lambda / lambda) * Height * Length
      WRITE(*,*) "  getting nodal weights"
      ! compute consistent nodal weights
      CALL ImportedMesh_SurfaceWeightsGet(nodalWeights,NumberOfNodesPerElement,Err)
      WRITE(*,*) "  applying consistent nodal forces"
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          ! apply consistent nodal forces based on nodal weights
          ! Neumann BC at x=lx
          NodalForce=NeumannBCvalue*nodalWeights(NodeIdx)/(NumberGlobalYElements*NumberGlobalZElements)
          CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_INTEGRATED,NodalForce,Err)
        END IF
        NodeIdx  = NodeIdx + 1
        IF(NodeIdx>9) NodeIdx = 1
      END DO
    ! apply Neumann point BC
    ELSE IF(bcType==2) THEN
      WRITE(*,*) "Applying traction BC.."
      ! corresponding traction value
      NeumannBCvalue = 2.0_CMISSRP * MooneyRivlin1 * (lambda - 1.0_CMISSRP / lambda / lambda)
      DO NodeIdx=StartIdx,StopIdx
        NodeNumber=BoundaryPatchesImport(NodeIdx)
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          ! Neumann BC at x=lx
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NodeNumber, &
            & 1,CMFE_BOUNDARY_CONDITION_NEUMANN_POINT,NeumannBCvalue,Err)
        END IF
      END DO
      CALL HandleError("Warning: BC type 2 not fully implemented for user-defined mesh.")
    END IF
  END IF

  ! finish BC
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  ! Set export file name
  IF(NumberGlobalXElements<10) THEN
    WRITE(Filename, "(A21,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A2,I1,A3,I1,A3,I1,A3,I1)") &
      & "results/current_run/l", &
      & INT(Width),"x",INT(Height),"x",INT(Length), &
      & "_n", &
      & NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements, &
      & "_i",displacementInterpolationType,"_s",SolverIsDirect, &
      & "_fd",JACOBIAN_FD,"_gm",useGeneratedMesh,"_bc",bcType
  ELSE IF((NumberGlobalXElements>9).AND.(NumberGlobalXElements<100)) THEN
    WRITE(Filename, "(A21,I1,A1,I1,A1,I1,A2,I2,A1,I1,A1,I1,A2,I1,A2,I1,A3,I1,A3,I1,A3,I1)") &
      & "results/current_run/l", &
      & INT(Width),"x",INT(Height),"x",INT(Length), &
      & "_n", &
      & NumberGlobalXElements,"x",NumberGlobalYElements,"x",NumberGlobalZElements, &
      & "_i",displacementInterpolationType,"_s",SolverIsDirect, &
      & "_fd",JACOBIAN_FD,"_gm",useGeneratedMesh,"_bc",bcType
  END IF
  ! make sure directories exist
  INQUIRE(file="./results/", exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir -p ./results/")
  END IF
  INQUIRE(file="./results/current_run/", exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir -p ./results/current_run/")
  END IF
  INQUIRE(file=trim(Filename), exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir -p "//trim(Filename))
  END IF

  ! Export solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  Filename=trim(Filename)//trim("/Example")
  CALL cmfe_Fields_NodesExport(Fields,Filename,"FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,Filename,"FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  IF(ALLOCATED(nodalWeights))           DEALLOCATE(nodalWeights)
  IF(ALLOCATED(NodesImport))            DEALLOCATE(NodesImport)
  IF(ALLOCATED(ElementsImport))         DEALLOCATE(ElementsImport)
  IF(ALLOCATED(BoundaryPatchesImport))  DEALLOCATE(BoundaryPatchesImport)
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(ErrorString)
    CHARACTER(LEN=*), INTENT(IN) :: ErrorString

    WRITE(*,'(">>ERROR: ",A)') ErrorString(1:LEN_TRIM(ErrorString))
    STOP

  END SUBROUTINE HandleError

END PROGRAM FORTRANEXAMPLE

