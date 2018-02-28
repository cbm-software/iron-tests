!> \file
!> \author MM
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls.
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
!> Contributor(s): Andreas Hessenthaler
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

!> \example ClassicalField/Laplace/GeneralisedLaplace/Fortran/src/FortranExample.f90
!! Example program to solve a Generalised Laplace equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/Laplace/GeneralisedLaplace/history.html
!!
!<

!> Main program
PROGRAM GENERALISEDLAPLACEEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7

  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=13

  INTEGER(CMISSIntg), PARAMETER :: FibreFieldNumberOfVariables=1

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11

  INTEGER(CMISSIntg), PARAMETER :: DerivativeUserNumber=1

  REAL(CMISSRP), PARAMETER :: PI=4.0_CMISSRP*DATAN(1.0_CMISSRP)

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg)    :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg)    :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg)    :: INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI,SOLVER_TYPE

  INTEGER(CMISSIntg)    :: node_idx,component_idx,TotalNumberOfNodes

  INTEGER(CMISSIntg)    :: FibreFieldNumberOfComponents

  CHARACTER(LEN=255)    :: COMMAND_ARGUMENT,Filename
  REAL(CMISSRP)         :: WIDTH,HEIGHT,LENGTH,SIGMA11,SIGMA22,SIGMA33
  REAL(CMISSRP)         :: FibreFieldAngles(3),ANGLE1,ANGLE2,ANGLE3

  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,FibreField,MaterialsField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain,NodeDomain
  INTEGER(CMISSIntg) :: Err
  
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
  IF(NUMBER_OF_ARGUMENTS >= 14) THEN
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
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(5,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 5.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Y_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(6,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 6.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Z_ELEMENTS<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    ! interpolation type (linear, quadratic
    CALL GET_COMMAND_ARGUMENT(7,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 7.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid interpolation specification.")
    ! solver type (direct, iterative)
    CALL GET_COMMAND_ARGUMENT(8,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 8.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SOLVER_TYPE
    IF((SOLVER_TYPE<0).OR.(SOLVER_TYPE>1)) CALL HANDLE_ERROR("Invalid solver type specification.")
    ! sigma_11
    CALL GET_COMMAND_ARGUMENT(9,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 9.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SIGMA11
    ! sigma_22
    CALL GET_COMMAND_ARGUMENT(10,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 10.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SIGMA22
    ! sigma_33
    CALL GET_COMMAND_ARGUMENT(11,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 11.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) SIGMA33
    ! angle 1
    CALL GET_COMMAND_ARGUMENT(12,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 12.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) ANGLE1
    ! angle 2
    CALL GET_COMMAND_ARGUMENT(13,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 13.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) ANGLE2
    ! angle 3
    CALL GET_COMMAND_ARGUMENT(14,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 14.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) ANGLE3
  ELSE
    !If there are not enough arguments default the problem specification
    WIDTH                       = 2.0_CMISSRP
    HEIGHT                      = 1.0_CMISSRP
    LENGTH                      = 1.0_CMISSRP
    NUMBER_GLOBAL_X_ELEMENTS    = 2
    NUMBER_GLOBAL_Y_ELEMENTS    = 1
    NUMBER_GLOBAL_Z_ELEMENTS    = 1
    INTERPOLATION_TYPE          = CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    SOLVER_TYPE                 = 0
    SIGMA11                     = 1.0_CMISSRP
    SIGMA22                     = 1.0_CMISSRP
    SIGMA33                     = 1.0_CMISSRP
    ANGLE1                      = 0.0_CMISSRP
    ANGLE2                      = 0.0_CMISSRP
    ANGLE3                      = 0.0_CMISSRP
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  CALL cmfe_RandomSeedsSet(9999,Err)
  
!  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "GeneralisedLaplace",NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,INTERPOLATION_TYPE
  
  CALL cmfe_OutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
    
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    NUMBER_OF_GAUSS_XI=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)
   
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
 
  !Destory the mesh now that we have decomposed it
  !CALL cmfe_Mesh_Destroy(Mesh,Err)
 
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Start to create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FibreFieldNumberOfVariables,Err)

  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    FibreFieldNumberOfComponents=2
  ELSE
    FibreFieldNumberOfComponents=3
  ENDIF
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FibreFieldNumberOfComponents,Err)

  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  ! Rotation Angles (in radiant!!)
  ! in 2D an entry in Angle(1) means rotated x-axis,
  !          entry in Angle(2) doesn't make sense, as rotates out of surface ...
  ! in 3D an entry in Angle(1) means rotated around z-axis,
  !          entry in Angle(2) means rotated around y-axis
  !          entry in Angle(3) means rotated around x-axis => no change
  ! 45° equivalent to pi/4, 90° equivalent to pi/2

  ! Set the fibre field angle
  FibreFieldAngles=(/ANGLE1,ANGLE2,ANGLE3/)   !! ToDo: 2D/3D ??

  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,TotalNumberOfNodes,Err)
  DO node_idx=1,TotalNumberOfNodes
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,FibreFieldNumberOfComponents
        CALL cmfe_Field_ParameterSetUpdateNode(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
          & DerivativeUserNumber,node_idx,component_idx,FibreFieldAngles(component_idx),Err)
      ENDDO
    ENDIF
  ENDDO

  !Start to create material field on the region
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_Field_CreateStart(MaterialsFieldUserNumber,Region,MaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_TypeSet(MaterialsField,CMFE_FIELD_MATERIAL_TYPE,Err)

  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(MaterialsField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialsField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialsField,1,Err)

  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    ! symmetric 2x2 tensor => 3 different entries
    ! 1 - 11, 2 - 22, 3 - 12=21
    CALL cmfe_Field_NumberOfComponentsSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  ELSE
    ! symmetric 3x3 tensor => 6 different entries
    ! 1 - 11, 2 - 22, 3 - 33, 4 - 12=21, 5 - 23=32, 6 - 13=31
    CALL cmfe_Field_NumberOfComponentsSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,6,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  ENDIF

  !Finish creating the field
  CALL cmfe_Field_CreateFinish(MaterialsField,Err)

  ! Set the material parameters, i.e. the conductivity values
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,SIGMA11,Err)  ! 11
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 2,SIGMA22,Err)  ! 22
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 3,0.0_CMISSRP,Err)  ! 12=21
  ELSE
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,SIGMA11,Err)  ! 11
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 2,SIGMA22,Err)  ! 22
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 3,SIGMA33,Err)  ! 33
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 4,0.0_CMISSRP,Err)  ! 12=21
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 5,0.0_CMISSRP,Err)  ! 23=32
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 6,0.0_CMISSRP,Err)  !13=31
  ENDIF

  !Create the Generalised Laplace Equations set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"ScalarField",Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"ScalarField (derivative)",Err)
  !Set the DOFs to be contiguous across components
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP, &
    & Err)

  !Create the equations set material field variables
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)
  
  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMFE_PROBLEM_GENERALISED_LAPLACE_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  
  IF(SOLVER_TYPE==0) THEN
    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)
  ENDIF
  
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Export results
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  WRITE(filename, "(A21,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A2,I1,A8)") &
    & "results/current_run/l", &
    & INT(WIDTH),"x",INT(HEIGHT),"x",INT(LENGTH), &
    & "_n", &
    & NUMBER_GLOBAL_X_ELEMENTS,"x",NUMBER_GLOBAL_Y_ELEMENTS,"x",NUMBER_GLOBAL_Z_ELEMENTS, &
    & "_i",INTERPOLATION_TYPE,"_s",SOLVER_TYPE,"/Example"
  filename=trim(filename)
  CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM GENERALISEDLAPLACEEXAMPLE
