module pflotran_model_module
  use petscsys
  use Simulation_Base_class
  use Realization_Base_class
  use Timestepper_Base_class
  use Option_module
  !use Input_module
  !use Init_module
  use Logging_module
  !use Stochastic_module
  !use Stochastic_Aux_module
  use Waypoint_module
  use Units_module
  use PFLOTRAN_Constants_module
  use Mapping_module

!#include "definitions.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscviewer.h"
!#include "finclude/petscvec.h90"

  use petscsys
  use petscvec

  implicit none


  ! Note:
  !
  ! ELM has the following:
  !   (i) 3D subsurface grid (ELM_SUB);
  !   (ii) 2D surface grid (ELM_SRF).
  ! ELM decomposes the 3D subsurface grid across processors in a 2D (i.e.
  ! cells in Z are not split across processors). Thus, the surface cells of
  ! 3D subsurface grid are on the same processors as the 2D surface grid.
  !
  ! PFLOTRAN has the following:
  !   (i) 3D subsurface grid (PF_SUB);
  !   (ii) surface control volumes of 3D subsurface grid (PF_2DSUB);
  !   (iii) 2D surface grid (PF_SRF).
  ! In PFLOTRAN, control volumes in PF_2DSUB and PF_SRF may reside on different
  ! processors. PF_SUB and PF_2DSUB are derived from simulation%realization;
  ! while PF_SRF refers to simulation%surf_realization.

  ! map level constants
  PetscInt, parameter, public :: ELM_SUB_TO_PF_SUB           = 1 ! 3D --> 3D
  PetscInt, parameter, public :: ELM_SUB_TO_PF_EXTENDED_SUB  = 2 ! 3D --> extended 3D
  PetscInt, parameter, public :: ELM_SRF_TO_PF_2DSUB         = 3 ! 2D --> SURF of 3D grid
  PetscInt, parameter, public :: PF_SUB_TO_ELM_SUB           = 5 ! 3D --> 3D

  ! mesh ids
  PetscInt, parameter, public :: ELM_SUB_MESH   = 1
  PetscInt, parameter, public :: ELM_SRF_MESH   = 2
  PetscInt, parameter, public :: PF_SUB_MESH    = 3
  PetscInt, parameter, public :: PF_2DSUB_MESH  = 4
  PetscInt, parameter, public :: PF_SRF_MESH    = 5

  type, public :: inside_each_overlapped_cell
     PetscInt           :: id
     PetscInt           :: ocell_count
     PetscInt,  pointer :: ocell_id(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_overlapped_cell

  type, public :: pflotran_model_type
    class(simulation_base_type),  pointer :: simulation
    type(option_type),      pointer :: option
    PetscReal :: pause_time_1
    PetscReal :: pause_time_2
    type(inside_each_overlapped_cell), pointer :: pf_cells(:)
    type(inside_each_overlapped_cell), pointer :: elm_cells(:)
    type(mapping_type),                pointer :: map_elm_sub_to_pf_sub
    type(mapping_type),                pointer :: map_elm_sub_to_pf_extended_sub
    type(mapping_type),                pointer :: map_elm_srf_to_pf_2dsub
    type(mapping_type),                pointer :: map_elm_srf_to_pf_srf
    type(mapping_type),                pointer :: map_pf_sub_to_elm_sub
    type(mapping_type),                pointer :: map_pf_srf_to_elm_srf

    PetscLogDouble :: timex(4), timex_wall(4)

  end type pflotran_model_type

  public::pflotranModelCreate,               &
       pflotranModelInitMapping,             &
       pflotranModelSetSoilProp,             &
       pflotranModelGetSoilProp,             &
       pflotranModelSetICs,                  &
       pflotranModelUpdateFlowConds,         &
       pflotranModelGetUpdatedData,          &
       pflotranModelStepperRunInit,          &
       pflotranModelStepperRunTillPauseTime, &
       pflotranModelInsertWaypoint,          &
       pflotranModelDeleteWaypoint,          &
       pflotranModelSetupRestart,            &
       pflotranModelStepperRunFinalize,      &
       pflotranModelStepperCheckpoint,       &
       pflotranModelNSurfCells3DDomain,      &
       pflotranModelGetTopFaceArea,          &
       pflotranModelDestroy

  private :: &
       pflotranModelSetupMappingFiles


contains

! ************************************************************************** !

  function pflotranModelCreate(mpicomm, pflotran_prefix)
  !
  ! Allocates and initializes the pflotranModel object.
  ! It performs the same sequence of commands as done in pflotran.F90
  ! before model integration is performed by the call to StepperRun()
  ! routine
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    use Communicator_Aux_module
    use Option_module
    use Simulation_Base_class
    use Simulation_Subsurface_class
    use Factory_PFLOTRAN_module
    use Factory_Subsurface_module, only : FactorySubsurfaceInitialize
    use Factory_Geomechanics_module
    use Factory_Forward_module
    use Driver_class

    implicit none

#include "petsc/finclude/petscsys.h"

    PetscInt, intent(in) :: mpicomm
    character(len=256), intent(in) :: pflotran_prefix

    type(pflotran_model_type), pointer :: pflotranModelCreate

    type(pflotran_model_type),      pointer :: model
    class(driver_type), pointer :: driver
    PetscErrorCode :: ierr

    allocate(model)

    nullify(model%simulation)
    nullify(model%option)

    driver => DriverCreate()
    call CommInitPetsc(driver%comm,mpicomm)

    ! NOTE(bja, 2013-07-19) GB's Hack to get communicator correctly
    ! setup on mpich/mac. should be generally ok, but may need an
    ! apple/mpich ifdef if it cause problems elsewhere.
    PETSC_COMM_SELF = MPI_COMM_SELF
    PETSC_COMM_WORLD = MPI_COMM_WORLD

    ! [yx 20240603] FactoryPFLOTRANInitialize will reset driver%input_filename to 'pflotran.in'
    ! a possible way to make pflotran_prefix work is to copy and modify code from FactoryPFLOTRANInitialize() here.
    ! ! NOTE(bja) 2013-06-25 : external driver must provide an input
    ! ! prefix string. If the driver wants to use pflotran.in, then it
    ! ! should explicitly request that with 'pflotran'.
    ! if (len(trim(pflotran_prefix)) > 1) then
    !   driver%input_prefix = trim(pflotran_prefix)
    !   driver%input_filename = trim(driver%input_prefix) // '.in'
    !   driver%global_prefix = driver%input_prefix
    ! else
    !   if (driver%IsIORank()) then
    !     print *, 'The external driver must provide the pflotran input &
    !              &file prefix.'
    !   endif
    !   stop
    ! end if

    call FactoryPFLOTRANInitialize(driver,model%simulation)
    select type(s=>model%simulation)
      class is(simulation_subsurface_type)
        model%option => s%option
    end select

    ! TODO(bja, 2013-07-15) this needs to be left alone for pflotran
    ! to deal with, or we need a valid unit number from the driver as
    ! a function parameter.
!!$    model%option%fid_out = 16

    model%pause_time_1 = -1.0d0
    model%pause_time_2 = -1.0d0

    ! NOTE(bja, 2013-07-15) needs to go before InitializeRun()...?
    call pflotranModelSetupMappingFiles(model)

    call model%simulation%InitializeRun()

    pflotranModelCreate => model

  end function pflotranModelCreate

! ************************************************************************** !

  subroutine pflotranModelSetupMappingFiles(model)
  !
  ! pflotranModelSetupMappingFiles
  ! create the mapping objects, reopen the input file and read the file names
  ! before model integration is performed by the call to StepperRun()
  ! routine
  ! NOTE(bja, 2013-07) this really needs to be moved out of pflotran
  ! ELM should be responsible for passing data in the correct
  ! format. That may require pflotran to provide a call back function
  ! for grid info.
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    use String_module
    use Option_module
    use Input_Aux_module
    use Mapping_module

    implicit none

#include "petsc/finclude/petscsys.h"

    type(pflotran_model_type), pointer, intent(inout) :: model
    type(input_type), pointer :: input

    type(mapping_type), pointer        :: map

    PetscBool :: elm2pf_flux_file
    PetscBool :: elm2pf_soil_file
    PetscBool :: elm2pf_gflux_file
    PetscBool :: elm2pf_rflux_file
    PetscBool :: pf2elm_flux_file
    PetscBool :: pf2elm_surf_file
    character(len=MAXSTRINGLENGTH) :: string
    character(len=MAXWORDLENGTH) :: word

    nullify(model%pf_cells)
    nullify(model%elm_cells)
    nullify(model%map_elm_sub_to_pf_sub)
    nullify(model%map_elm_sub_to_pf_extended_sub)
    nullify(model%map_elm_srf_to_pf_2dsub)
    nullify(model%map_elm_srf_to_pf_srf)
    nullify(model%map_pf_sub_to_elm_sub)
    nullify(model%map_pf_srf_to_elm_srf)

    model%map_elm_sub_to_pf_sub          => MappingCreate()
    model%map_elm_sub_to_pf_extended_sub => MappingCreate()
    model%map_elm_srf_to_pf_2dsub        => MappingCreate()
    model%map_elm_srf_to_pf_srf          => MappingCreate()
    model%map_pf_sub_to_elm_sub          => MappingCreate()
    model%map_pf_srf_to_elm_srf          => MappingCreate()

    input => InputCreate(15, &
                    model%option%input_filename, model%option)

    ! Read names of mapping file
    elm2pf_flux_file   = PETSC_FALSE
    elm2pf_soil_file   = PETSC_FALSE
    elm2pf_gflux_file  = PETSC_FALSE
    elm2pf_rflux_file  = PETSC_FALSE
    pf2elm_flux_file   = PETSC_FALSE
    pf2elm_surf_file   = PETSC_FALSE

    string = "MAPPING_FILES"
    call InputFindStringInFile(input,model%option,string)
    call InputPushBlock(input,'MAPPING_FILES',model%option)
    do
      call InputReadPflotranString(input, model%option)
      if (InputCheckExit(input, model%option)) exit
      if (input%ierr /= 0) exit

      call InputReadCard(input, model%option, word)
      call InputErrorMsg(input, model%option, 'keyword', 'MAPPING_FILES')
      call StringToUpper(word)

      select case(trim(word))
        case('ELM2PF_FLUX_FILE')
          map => model%map_elm_sub_to_pf_sub
          call InputReadNChars(input, model%option, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, model%option, 'type', 'MAPPING_FILES')
          map%filename     = trim(map%filename)//CHAR(0)
          elm2pf_flux_file = PETSC_TRUE

        case('ELM2PF_SOIL_FILE')
          map => model%map_elm_sub_to_pf_extended_sub
          call InputReadNChars(input, model%option, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, model%option, 'type', 'MAPPING_FILES')
          map%filename     = trim(map%filename)//CHAR(0)
          elm2pf_soil_file = PETSC_TRUE

        case('ELM2PF_GFLUX_FILE')
          map => model%map_elm_srf_to_pf_2dsub
          call InputReadNChars(input, model%option, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, model%option, 'type', 'MAPPING_FILES')
          map%filename      = trim(map%filename)//CHAR(0)
          elm2pf_gflux_file = PETSC_TRUE

        case('ELM2PF_RFLUX_FILE')
          map => model%map_elm_srf_to_pf_srf
          call InputReadNChars(input, model%option, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, model%option, 'type', 'MAPPING_FILES')
          map%filename      = trim(map%filename)//CHAR(0)
          elm2pf_rflux_file = PETSC_TRUE

        case('PF2ELM_SURF_FILE')
          map => model%map_pf_srf_to_elm_srf
          call InputReadNChars(input, model%option, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, model%option, 'type', 'MAPPING_FILES')
          map%filename     = trim(map%filename)//CHAR(0)
          pf2elm_surf_file = PETSC_TRUE

        case('PF2ELM_FLUX_FILE')
          map => model%map_pf_sub_to_elm_sub
          call InputReadNChars(input, model%option, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE)
          call InputErrorMsg(input, model%option, 'type', 'MAPPING_FILES')
          map%filename     = trim(map%filename)//CHAR(0)
          pf2elm_flux_file = PETSC_TRUE

        case default
          model%option%io_buffer='Keyword ' // trim(word) // &
            ' in input file not recognized'
          call PrintErrMsg(model%option)
      end select

      ! Read mapping file
      if (index(map%filename, '.h5') > 0) then
        call MappingReadHDF5(map, map%filename, model%option)
      else
        call MappingReadTxtFile(map, map%filename, model%option)
      endif

    enddo
    call InputPopBlock(input,model%option)
    call InputDestroy(input)

    if ((.not. elm2pf_soil_file) .or. (.not. elm2pf_flux_file) .or. &
        (.not. pf2elm_flux_file) ) then
      model%option%io_buffer='One of the mapping files not found'
      call PrintErrMsg(model%option)
    endif

    if(model%option%iflowmode==TH_MODE.and.(.not.elm2pf_gflux_file)) then
      model%option%io_buffer='Running in TH_MODE without a ELM2PF_GFLUX_FILE'
      call PrintErrMsg(model%option)
    endif

  end subroutine pflotranModelSetupMappingFiles

! ************************************************************************** !

  subroutine pflotranModelStepperRunInit(model)
  !
  ! It performs the same execution of commands
  ! that are carried out in StepperRun() before the model integration
  ! begins over the entire simulation time interval
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    type(pflotran_model_type), pointer, intent(inout) :: model

    call model%simulation%InitializeRun()

  end subroutine pflotranModelStepperRunInit

! ************************************************************************** !

  subroutine pflotranModelStepperCheckpoint(model, id_stamp)
  !
  ! wrapper around StepperCheckpoint
  ! NOTE(bja, 2013-06-27) : the date stamp from elm is 32 characters
  !

    use Option_module
    use Simulation_Subsurface_class

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXSTRINGLENGTH), intent(in) :: id_stamp
    PetscViewer :: viewer

    select type(sim => model%simulation)
      class is(simulation_subsurface_type)
        call sim%process_model_coupler_list%CheckpointBinary(viewer,id_stamp)
    end select

  end subroutine pflotranModelStepperCheckpoint

! ************************************************************************** !

subroutine pflotranModelSetICs(pflotran_model)
  !
  ! Set initial conditions
  !
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  !

    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use Field_module
    use elm_pflotran_interface_data
    use Global_Aux_module
    use Discretization_module
    use Richards_module
    use TH_module
    use Option_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer           :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(richards_auxvar_type), pointer       :: aux_var
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(simulation_base_type), pointer :: simulation

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: xx_loc_p(:)

    PetscScalar, pointer :: press_pf_loc(:) ! Pressure [Pa]

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelSetICs only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field
    global_aux_vars  => patch%aux%Global%auxvars

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_sub, &
                                    elm_pf_idata%press_elm, &
                                    elm_pf_idata%press_pf)

    if (pflotran_model%option%iflowmode .ne. RICHARDS_MODE) then
        pflotran_model%option%io_buffer='pflotranModelSetICs ' // &
          'not implmented for this mode.'
        call PrintErrMsg(pflotran_model%option)
    endif

    call VecGetArrayF90(field%flow_xx,xx_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%press_pf,press_pf_loc, &
                        ierr);CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (associated(patch%imat)) then
          if (patch%imat(ghosted_id) <= 0) cycle
       endif
       xx_loc_p(ghosted_id)=press_pf_loc(local_id)
    enddo

    call VecRestoreArrayF90(field%flow_xx,xx_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%press_pf,press_pf_loc, &
                            ierr);CHKERRQ(ierr)

    ! update dependent vectors: Saturation
    call DiscretizationGlobalToLocal(realization%discretization, field%flow_xx, &
         field%flow_xx_loc, NFLOWDOF)
    call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
      case default
        pflotran_model%option%io_buffer='pflotranModelSetICs ' // &
          'not implmented for this mode.'
        call PrintErrMsg(pflotran_model%option)
    end select

end subroutine pflotranModelSetICs

! ************************************************************************** !

  subroutine pflotranModelSetSoilProp(pflotran_model)
  !
  ! Converts hydraulic properties from ELM units
  ! into PFLOTRAN units.
  ! #ifdef ELM_PFLOTRAN
  !
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  !

    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use TH_Aux_module
    use Field_module
    use Option_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use elm_pflotran_interface_data
    use Mapping_module
    use Material_module
    use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer           :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(richards_auxvar_type), pointer       :: rich_aux_var
    type(th_auxvar_type), pointer             :: th_aux_vars(:)
    type(th_auxvar_type), pointer             :: th_aux_var
    type(simulation_base_type), pointer :: simulation

    PetscErrorCode     :: ierr
    PetscInt           :: ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: bc_lambda, bc_alpha
    Vec                :: porosity_loc, perm_xx_loc, perm_yy_loc, perm_zz_loc

    PetscScalar, pointer :: hksat_x_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw_pf_loc(:)     ! Clapp and Hornberger "b"

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/S^2]

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelSetSoilProp only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field

    select case(pflotran_model%option%iflowmode)
      case(RICHARDS_MODE)
        rich_aux_vars   => patch%aux%Richards%auxvars
      case(TH_MODE)
        th_aux_vars   => patch%aux%TH%auxvars
      case default
        call PrintErrMsg(pflotran_model%option, &
          'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp')
    end select

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_extended_sub, &
                                    elm_pf_idata%hksat_x_elm, &
                                    elm_pf_idata%hksat_x_pf)

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_extended_sub, &
                                    elm_pf_idata%hksat_y_elm, &
                                    elm_pf_idata%hksat_y_pf)

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_extended_sub, &
                                    elm_pf_idata%hksat_z_elm, &
                                    elm_pf_idata%hksat_z_pf)

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_extended_sub, &
                                    elm_pf_idata%sucsat_elm, &
                                    elm_pf_idata%sucsat_pf)

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_extended_sub, &
                                    elm_pf_idata%bsw_elm, &
                                    elm_pf_idata%bsw_pf)

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_extended_sub, &
                                    elm_pf_idata%watsat_elm, &
                                    elm_pf_idata%watsat_pf)

    call VecGetArrayF90(elm_pf_idata%hksat_x_pf,hksat_x_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%hksat_y_pf,hksat_y_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%hksat_z_pf,hksat_z_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%sucsat_pf,sucsat_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%watsat_pf,watsat_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%bsw_pf,bsw_pf_loc,ierr);CHKERRQ(ierr)

    call VecDuplicate(field%work_loc,porosity_loc,ierr);CHKERRQ(ierr)
    call VecDuplicate(field%work_loc,perm_xx_loc,ierr);CHKERRQ(ierr)
    call VecDuplicate(field%work_loc,perm_yy_loc,ierr);CHKERRQ(ierr)
    call VecDuplicate(field%work_loc,perm_zz_loc,ierr);CHKERRQ(ierr)

    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,porosity_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,perm_xx_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,perm_yy_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,perm_zz_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)

    call VecGetArrayF90(porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(perm_yy_loc,perm_yy_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(perm_zz_loc,perm_zz_loc_p,ierr);CHKERRQ(ierr)

    do ghosted_id = 1, grid%ngmax

      ! bc_alpha [1/Pa]; while sucsat [mm of H20]
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      bc_alpha = 1.d0/(sucsat_pf_loc(ghosted_id)*grav)

      ! bc_lambda = 1/bsw
      bc_lambda = 1.d0/bsw_pf_loc(ghosted_id)

      select case(pflotran_model%option%iflowmode)
        case(RICHARDS_MODE)
          rich_aux_var => rich_aux_vars(ghosted_id)
          rich_aux_var%bc_alpha = bc_alpha
          rich_aux_var%bc_lambda = bc_lambda
        case(TH_MODE)
          th_aux_var => th_aux_vars(ghosted_id)
          th_aux_var%bc_alpha = min(bc_alpha,10.d-4)
          th_aux_var%bc_lambda = bc_lambda
      end select

      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      perm_xx_loc_p(ghosted_id) = hksat_x_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
      perm_yy_loc_p(ghosted_id) = hksat_y_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0
      perm_zz_loc_p(ghosted_id) = hksat_z_pf_loc(ghosted_id)*vis/(den*grav)/1000.d0

      porosity_loc_p(ghosted_id) = watsat_pf_loc(ghosted_id)

    enddo

    call VecRestoreArrayF90(elm_pf_idata%hksat_x_pf,hksat_x_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%hksat_y_pf,hksat_y_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%hksat_z_pf,hksat_z_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%sucsat_pf,sucsat_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%watsat_pf,watsat_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%bsw_pf,bsw_pf_loc, &
                            ierr);CHKERRQ(ierr)

    call VecRestoreArrayF90(porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(perm_yy_loc,perm_yy_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(perm_zz_loc,perm_zz_loc_p,ierr);CHKERRQ(ierr)

    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,porosity_loc, &
                                 POROSITY,ZERO_INTEGER)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,perm_xx_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,perm_yy_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,perm_zz_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)

    call VecDestroy(porosity_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(perm_xx_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(perm_yy_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(perm_zz_loc,ierr);CHKERRQ(ierr)

  end subroutine pflotranModelSetSoilProp

! ************************************************************************** !

  subroutine pflotranModelGetSoilProp(pflotran_model)
  !
  ! Converts hydraulic properties from PFLOTRAN units for ELM units
  !
  ! Author: Gautam Bisht
  ! Date: 10/27/2014
  !

    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Richards_Aux_module
    use TH_Aux_module
    use Field_module
    use Option_module

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use elm_pflotran_interface_data
    use Mapping_module
    use Material_module
    use Material_Aux_module
    use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z
    use Saturation_Function_module
    use Characteristic_Curves_module
    use Characteristic_Curves_Common_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(field_type), pointer                 :: field
    type(richards_auxvar_type), pointer       :: rich_aux_vars(:)
    type(richards_auxvar_type), pointer       :: rich_aux_var
    type(th_auxvar_type), pointer             :: th_aux_vars(:)
    type(th_auxvar_type), pointer             :: th_aux_var
    type(simulation_base_type), pointer       :: simulation
    type(saturation_function_type)            :: saturation_function
    class(characteristic_curves_type), pointer:: characteristic_curve

    PetscErrorCode     :: ierr
    PetscInt           :: ghosted_id, local_id
    PetscInt , parameter :: liq_iphase = 1
    PetscReal          :: den, vis, grav, Sr
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: bc_lambda, bc_alpha
    Vec                :: porosity_loc, perm_xx_loc, perm_yy_loc, perm_zz_loc

    PetscScalar, pointer :: hksat_x2_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y2_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z2_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat2_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat2_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw2_pf_loc(:)     ! Clapp and Hornberger "b"
    PetscScalar, pointer :: thetares2_pf_loc(:)! residual soil mosture = sat_res * por

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/S^2]

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: pflotranModelSetSoilProp only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select

    patch           => realization%patch
    grid            => patch%grid
    field           => realization%field

    select case(pflotran_model%option%iflowmode)
      case(RICHARDS_MODE)
        rich_aux_vars   => patch%aux%Richards%auxvars
      case(TH_MODE)
        th_aux_vars   => patch%aux%TH%auxvars
      case default
        call PrintErrMsg(pflotran_model%option, &
          'Current PFLOTRAN mode not supported by pflotranModelSetSoilProp')
    end select

    call VecDuplicate(field%work_loc,porosity_loc,ierr);CHKERRQ(ierr)
    call VecDuplicate(field%work_loc,perm_xx_loc,ierr);CHKERRQ(ierr)
    call VecDuplicate(field%work_loc,perm_yy_loc,ierr);CHKERRQ(ierr)
    call VecDuplicate(field%work_loc,perm_zz_loc,ierr);CHKERRQ(ierr)

    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,porosity_loc, &
                                 POROSITY,POROSITY_CURRENT)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,perm_xx_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,perm_yy_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,perm_zz_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)

    call VecGetArrayF90(elm_pf_idata%hksat_x2_pf,hksat_x2_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%hksat_y2_pf,hksat_y2_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%hksat_z2_pf,hksat_z2_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%sucsat2_pf,sucsat2_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%watsat2_pf,watsat2_pf_loc, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%bsw2_pf,bsw2_pf_loc,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%thetares2_pf,thetares2_pf_loc, &
                        ierr);CHKERRQ(ierr)

    call VecGetArrayF90(porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(perm_yy_loc,perm_yy_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(perm_zz_loc,perm_zz_loc_p,ierr);CHKERRQ(ierr)

    hksat_x2_pf_loc(:) = 0.d0
    hksat_y2_pf_loc(:) = 0.d0
    hksat_z2_pf_loc(:) = 0.d0
    sucsat2_pf_loc(:) = 0.d0
    watsat2_pf_loc(:) = 0.d0
    bsw2_pf_loc(:) = 0.d0
    thetares2_pf_loc(:) = 0.d0

    do ghosted_id = 1, grid%ngmax

      if (patch%cc_id(ghosted_id) < 1) cycle

      local_id = grid%nG2L(ghosted_id)

      select case(pflotran_model%option%iflowmode)
      case(RICHARDS_MODE)
         characteristic_curve => patch%characteristic_curves_array(patch%cc_id(ghosted_id))%ptr

         select type(sf => characteristic_curve%saturation_function)
         class is(sat_func_VG_type)
            bc_lambda = sf%m
            bc_alpha  = sf%alpha
            Sr        = sf%Sr
         class is(sat_func_BC_type)
            bc_lambda = sf%lambda
            bc_alpha  = sf%alpha
            Sr        = sf%Sr
         class default
            pflotran_model%option%io_buffer = 'ELM-PFLOTRAN only supports ' // &
                 'sat_func_VG_type and sat_func_BC_type'
            call PrintErrMsg(pflotran_model%option)
         end select
      case(TH_MODE)
         saturation_function = patch%saturation_function_array(patch%cc_id(ghosted_id))%ptr
         th_aux_var => th_aux_vars(ghosted_id)
         bc_alpha  = saturation_function%alpha
         select case(saturation_function%saturation_function_itype)
         case(VAN_GENUCHTEN)
            bc_lambda = saturation_function%m
         case(BROOKS_COREY)
            bc_lambda = saturation_function%lambda
         end select
         Sr = saturation_function%Sr(liq_iphase)
      end select

      ! bc_alpha [1/Pa]; while sucsat [mm of H20]
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      sucsat2_pf_loc(local_id) = 1.d0/(bc_alpha*grav)

      ! bc_lambda = 1/bsw
      bsw2_pf_loc(local_id) = 1.d0/bc_lambda

      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      hksat_x2_pf_loc(local_id) = perm_xx_loc_p(ghosted_id)/vis*(den*grav)*1000.d0
      hksat_y2_pf_loc(local_id) = perm_yy_loc_p(ghosted_id)/vis*(den*grav)*1000.d0
      hksat_z2_pf_loc(local_id) = perm_zz_loc_p(ghosted_id)/vis*(den*grav)*1000.d0

      watsat2_pf_loc(local_id) = porosity_loc_p(ghosted_id)


      thetares2_pf_loc(local_id) = porosity_loc_p(ghosted_id)*Sr

   enddo

    call VecRestoreArrayF90(elm_pf_idata%hksat_x2_pf,hksat_x2_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%hksat_y2_pf,hksat_y2_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%hksat_z2_pf,hksat_z2_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%sucsat2_pf,sucsat2_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%watsat2_pf,watsat2_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%bsw2_pf,bsw2_pf_loc, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%thetares2_pf,thetares2_pf_loc, &
                            ierr);CHKERRQ(ierr)

    call VecRestoreArrayF90(porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(perm_xx_loc,perm_xx_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(perm_yy_loc,perm_yy_loc_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(perm_zz_loc,perm_zz_loc_p,ierr);CHKERRQ(ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%hksat_x2_pf, &
                                    elm_pf_idata%hksat_x2_elm)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%hksat_y2_pf, &
                                    elm_pf_idata%hksat_y2_elm)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%hksat_z2_pf, &
                                    elm_pf_idata%hksat_z2_elm)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%sucsat2_pf, &
                                    elm_pf_idata%sucsat2_elm)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%bsw2_pf, &
                                    elm_pf_idata%bsw2_elm)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%watsat2_pf, &
                                    elm_pf_idata%watsat2_elm)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%thetares2_pf, &
                                    elm_pf_idata%thetares2_elm)

    call VecDestroy(porosity_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(perm_xx_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(perm_yy_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(perm_zz_loc,ierr);CHKERRQ(ierr)

  end subroutine pflotranModelGetSoilProp

! ************************************************************************** !

  subroutine pflotranModelInitMapping(pflotran_model,  &
                                      grid_elm_cell_ids_nindex, &
                                      grid_elm_npts_local, &
                                      map_id)
  !
  ! #endif
  ! Initialize mapping between the two model grid
  ! (ELM and PFLTORAN)
  !
  ! Author: Gautam Bisht
  ! Date: 03/24/2011
  !

    use Input_Aux_module
    use Option_module
    use Realization_Subsurface_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_elm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_elm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename

    select case (map_id)
      case (ELM_SUB_TO_PF_SUB, ELM_SUB_TO_PF_EXTENDED_SUB, PF_SUB_TO_ELM_SUB)
        call pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_elm_cell_ids_nindex, &
                                      grid_elm_npts_local, &
                                      map_id)
      case (ELM_SRF_TO_PF_2DSUB)
        call pflotranModelInitMapSrfTo2DSub(pflotran_model,  &
                                            grid_elm_cell_ids_nindex, &
                                            grid_elm_npts_local, &
                                            map_id)
      case default
        pflotran_model%option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMapping'
        call PrintErrMsg(pflotran_model%option)
    end select

  end subroutine pflotranModelInitMapping

! ************************************************************************** !

  subroutine pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_elm_cell_ids_nindex, &
                                      grid_elm_npts_local, &
                                      map_id)
  !
  ! Initialize mapping between 3D subsurface
  ! ELM grid and 3D subsurface PFLOTRAN grid.
  !
  ! Author: Gautam Bisht
  ! Date: 11/09/2013
  !

    use Input_Aux_module
    use Option_module
    use Realization_Subsurface_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_elm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_elm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename

    ! local
    PetscInt                           :: local_id, ghosted_id
    PetscInt                           :: grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_elm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_elm_local_nindex(:)

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_subsurface_type), pointer    :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    type(connection_set_type), pointer :: cur_connection_set

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    option          => pflotran_model%option
    patch           => realization%patch
    grid            => patch%grid

    !
    ! Mapping to/from entire PFLOTRAN 3D subsurface domain
    !

    ! Choose the appriopriate map
    select case(map_id)
      case(ELM_SUB_TO_PF_SUB)
        map => pflotran_model%map_elm_sub_to_pf_sub
        source_mesh_id = ELM_SUB_MESH
        dest_mesh_id = PF_SUB_MESH
      case(ELM_SUB_TO_PF_EXTENDED_SUB)
        map => pflotran_model%map_elm_sub_to_pf_extended_sub
        source_mesh_id = ELM_SUB_MESH
        dest_mesh_id = PF_SUB_MESH
      case(PF_SUB_TO_ELM_SUB)
        map => pflotran_model%map_pf_sub_to_elm_sub
        source_mesh_id = PF_SUB_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to pflotranModelInitMapping'
        call PrintErrMsg(option)
    end select

    grid_elm_npts_ghost=0

    ! Allocate memory to identify if ELM cells are local or ghosted.
    ! Note: Presently all ELM cells are local
    allocate(grid_elm_local_nindex(grid_elm_npts_local))
    do local_id = 1, grid_elm_npts_local
      grid_elm_local_nindex(local_id) = 1 ! LOCAL
    enddo

    ! Find cell IDs for PFLOTRAN grid
    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = grid%ngmax - grid%nlmax

    allocate(grid_pf_local_nindex(grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) == 0) then
        grid_pf_local_nindex(ghosted_id) = 0 ! GHOST
      else
        grid_pf_local_nindex(ghosted_id) = 1 ! LOCAL
      endif
    enddo

    select case(source_mesh_id)
      case(ELM_SUB_MESH)

        allocate(grid_pf_cell_ids_nindex(grid%ngmax))
        do ghosted_id = 1, grid%ngmax
          grid_pf_cell_ids_nindex(ghosted_id) = grid%nG2A(ghosted_id)-1
        enddo

        call MappingSetSourceMeshCellIds(map, grid_elm_npts_local, &
                                         grid_elm_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SUB_MESH)

        allocate(grid_pf_cell_ids_nindex(grid%nlmax))
        do ghosted_id = 1, grid%ngmax
          local_id = grid%nG2L(ghosted_id)
          if (local_id > 0) then
            grid_pf_cell_ids_nindex(local_id) = grid%nG2A(ghosted_id)-1
          endif
        enddo

        call MappingSetSourceMeshCellIds(map, grid_pf_npts_local, &
                                        grid_pf_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, grid_elm_npts_local, &
                                              grid_elm_npts_ghost, &
                                              grid_elm_cell_ids_nindex, &
                                              grid_elm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to pflotranModelInitMapping'
        call PrintErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_elm_local_nindex)

    call MappingDecompose(map, option%mycomm)
    call MappingFindDistinctSourceMeshCellIds(map)
    call MappingCreateWeightMatrix(map, option%myrank)
    call MappingCreateScatterOfSourceMesh(map, option%mycomm)

  end subroutine pflotranModelInitMappingSub2Sub

! ************************************************************************** !

  subroutine pflotranModelInitMapSrfTo2DSub(pflotran_model,  &
                                            grid_elm_cell_ids_nindex, &
                                            grid_elm_npts_local, &
                                            map_id)
  !
  ! This routine maps ELM surface grid onto surface of PFLOTRAN 3D subsurface
  ! grid.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/09/13
  !

    use Input_Aux_module
    use Option_module
    use Realization_Subsurface_class
    use Grid_module
    use Patch_module
    use Coupler_module
    use Connection_module
    use String_module
    use elm_pflotran_interface_data
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_elm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_elm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename

    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_elm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_elm_local_nindex(:)
    PetscInt, pointer                  :: grid_elm_cell_ids_nindex_copy(:)
    PetscInt                           :: count
    PetscInt                           :: sum_connection
    PetscInt                           :: ghosted_id
    PetscInt                           :: iconn
    PetscInt                           :: istart
    PetscInt, pointer                  :: int_array(:)
    PetscBool                          :: found
    PetscScalar,pointer                :: v_loc(:)
    PetscErrorCode                     :: ierr

    Vec                                :: surf_ids
    Vec                                :: surf_ids_loc
    IS                                 :: is_from
    IS                                 :: is_to
    VecScatter                         :: vec_scat

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_subsurface_type), pointer    :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(coupler_type), pointer        :: boundary_condition
    type(coupler_type), pointer        :: source_sink
    type(connection_set_type), pointer :: cur_connection_set

    option          => pflotran_model%option

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select

    allocate(grid_elm_cell_ids_nindex_copy(grid_elm_npts_local))
    grid_elm_cell_ids_nindex_copy = grid_elm_cell_ids_nindex

    ! Choose the appriopriate map
    select case(map_id)
      case(ELM_SRF_TO_PF_2DSUB)
        map => pflotran_model%map_elm_srf_to_pf_2dsub
        source_mesh_id = ELM_SUB_MESH
        dest_mesh_id = PF_2DSUB_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMappingSurf3D'
        call PrintErrMsg(option)
    end select

    grid_elm_npts_ghost=0

    ! Allocate memory to identify if ELM cells are local or ghosted.
    ! Note: Presently all ELM cells are local
    allocate(grid_elm_local_nindex(grid_elm_npts_local))
    do local_id = 1, grid_elm_npts_local
      grid_elm_local_nindex(local_id) = 1 ! LOCAL
    enddo

    ! Mapping to/from surface of PFLOTRAN domain
    found=PETSC_FALSE
    grid_pf_npts_local = 0
    grid_pf_npts_ghost = 0

    select case (dest_mesh_id)

      case(PF_2DSUB_MESH)

        patch => realization%patch
        grid => patch%grid

        ! Destination mesh is PF_2DSUB_MESH
        boundary_condition => patch%boundary_condition_list%first
        sum_connection = 0
        do
          if (.not.associated(boundary_condition)) exit
          cur_connection_set => boundary_condition%connection_set

          if(StringCompare(boundary_condition%name,'elm_gflux_bc')) then

            found=PETSC_TRUE

            ! Allocate memory to save cell ids and flag for local cells
            allocate(grid_pf_cell_ids_nindex(cur_connection_set%num_connections))
            allocate(grid_pf_local_nindex(cur_connection_set%num_connections))
            grid_pf_npts_local = cur_connection_set%num_connections

            ! Save cell ids in application order 0-based
            do iconn=1,cur_connection_set%num_connections
              sum_connection = sum_connection + 1
              local_id = cur_connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) cycle
              grid_pf_cell_ids_nindex(iconn) = grid%nG2A(ghosted_id) - 1
              grid_pf_local_nindex(iconn) = 1
            enddo
          else
            sum_connection = sum_connection + cur_connection_set%num_connections
          endif
          boundary_condition => boundary_condition%next
        enddo

        ! Setting the number of cells constituting the surface of the 3D
        ! subsurface domain for each model.
        elm_pf_idata%nlelm_2dsub = grid_elm_npts_local
        elm_pf_idata%ngelm_2dsub = grid_elm_npts_local
        elm_pf_idata%nlpf_2dsub  = grid_pf_npts_local
        elm_pf_idata%ngpf_2dsub  = grid_pf_npts_local

      case default
        option%io_buffer='Unknown source mesh'
        call PrintErrMsg(option)

    end select

    if(.not.found) then
      pflotran_model%option%io_buffer = 'elm_gflux_bc not found in boundary conditions'
      call PrintErrMsg(pflotran_model%option)
    endif

    !
    ! Step-1: Find surface cells-ids of PFLOTRAN subsurface domain
    !
    allocate(v_loc(grid_pf_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF,grid_pf_npts_local,surf_ids_loc, &
                      ierr);CHKERRQ(ierr)
    call VecCreateMPI(option%mycomm,grid%nlmax,PETSC_DECIDE,surf_ids, &
                      ierr);CHKERRQ(ierr)
    call VecSet(surf_ids,-1.d0,ierr);CHKERRQ(ierr)

    ! Set 1.0 to all cells that make up surface of PFLOTRAN subsurface domain
    call VecSetValues(surf_ids,grid_pf_npts_local,grid_pf_cell_ids_nindex, &
                      v_loc,INSERT_VALUES,ierr);CHKERRQ(ierr)
    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(surf_ids,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids,v_loc,ierr);CHKERRQ(ierr)
    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                    option%mycomm,ierr);CHKERRQ(ierr)

    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids,v_loc,ierr);CHKERRQ(ierr)

    !
    allocate(int_array(grid_pf_npts_local))
    do iconn = 1, grid_pf_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm,grid_pf_npts_local,int_array, &
                         PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm,grid_pf_npts_local, &
                         grid_pf_cell_ids_nindex,PETSC_COPY_VALUES,is_from, &
                         ierr);CHKERRQ(ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids,is_from,surf_ids_loc,is_to,vec_scat, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_from,ierr);CHKERRQ(ierr)
    call ISDestroy(is_to,ierr);CHKERRQ(ierr)

    call VecScatterBegin(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                         SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    count = 0
    do iconn = 1, grid_pf_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_pf_cell_ids_nindex(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc,ierr);CHKERRQ(ierr)

    !
    ! Step-2: Recompute 'map%s2d_iscr'
    !
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_nwts,surf_ids_loc, &
                      ierr);CHKERRQ(ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm,map%s2d_nwts,int_array, &
                         PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)


    do iconn = 1, map%s2d_nwts
      int_array(iconn) = map%s2d_icsr(iconn)
    enddo
    call ISCreateGeneral(option%mycomm,map%s2d_nwts,int_array, &
                         PETSC_COPY_VALUES,is_from,ierr);CHKERRQ(ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids,is_from,surf_ids_loc,is_to,vec_scat, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_from,ierr);CHKERRQ(ierr)
    call ISDestroy(is_to,ierr);CHKERRQ(ierr)

    call VecScatterBegin(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                         SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        map%s2d_icsr(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc,ierr);CHKERRQ(ierr)

    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied.'
      call PrintErrMsg(option)
    endif
    call VecDestroy(surf_ids,ierr);CHKERRQ(ierr)

    !
    ! Step-3: Find surface cells-ids of ELM subsurface domain
    !
    allocate(v_loc(grid_elm_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF,grid_elm_npts_local,surf_ids_loc, &
                      ierr);CHKERRQ(ierr)
    call VecCreateMPI(option%mycomm,elm_pf_idata%nlelm_sub,PETSC_DECIDE, &
                      surf_ids,ierr);CHKERRQ(ierr)
    call VecSet(surf_ids,-1.d0,ierr);CHKERRQ(ierr)

    ! Set 1.0 to all cells that make up surface of ELM subsurface domain
    call VecSetValues(surf_ids,grid_elm_npts_local, &
                      grid_elm_cell_ids_nindex_copy,v_loc,INSERT_VALUES, &
                      ierr);CHKERRQ(ierr)

    deallocate(v_loc)
    call VecAssemblyBegin(surf_ids,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(surf_ids,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids,v_loc,ierr);CHKERRQ(ierr)
    count = 0
    do local_id=1,elm_pf_idata%nlelm_sub
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count,istart,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                    option%mycomm,ierr);CHKERRQ(ierr)

    count = 0
    do local_id=1,elm_pf_idata%nlelm_sub
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(surf_ids,v_loc,ierr);CHKERRQ(ierr)

    !
    allocate(int_array(grid_elm_npts_local))
    do iconn = 1, grid_elm_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm,grid_elm_npts_local,int_array, &
                         PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm,grid_elm_npts_local, &
                         grid_elm_cell_ids_nindex_copy,PETSC_COPY_VALUES, &
                         is_from,ierr);CHKERRQ(ierr)


    ! create scatter context
    call VecScatterCreate(surf_ids,is_from,surf_ids_loc,is_to,vec_scat, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_from,ierr);CHKERRQ(ierr)
    call ISDestroy(is_to,ierr);CHKERRQ(ierr)

    call VecScatterBegin(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                         SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    count = 0
    do iconn = 1, grid_elm_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_elm_cell_ids_nindex_copy(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc,ierr);CHKERRQ(ierr)


    !
    ! Step-4: Recompute 'map%s2d_jscr'
    !
    call VecCreateSeq(PETSC_COMM_SELF,map%s2d_nwts,surf_ids_loc, &
                      ierr);CHKERRQ(ierr)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm,map%s2d_nwts,int_array, &
                         PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)


    do iconn = 1, map%s2d_nwts
      int_array(iconn) = map%s2d_jcsr(iconn)
    enddo
    call ISCreateGeneral(option%mycomm,map%s2d_nwts,int_array, &
                         PETSC_COPY_VALUES,is_from,ierr);CHKERRQ(ierr)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(surf_ids,is_from,surf_ids_loc,is_to,vec_scat, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_from,ierr);CHKERRQ(ierr)
    call ISDestroy(is_to,ierr);CHKERRQ(ierr)

    call VecScatterBegin(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                         SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(vec_scat,surf_ids,surf_ids_loc,INSERT_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterDestroy(vec_scat,ierr);CHKERRQ(ierr)

    call VecGetArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        map%s2d_jcsr(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(surf_ids_loc,v_loc,ierr);CHKERRQ(ierr)
    call VecDestroy(surf_ids_loc,ierr);CHKERRQ(ierr)

    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of surface cells in mapping dataset does not ' // &
        'match surface cells on which BC is applied.'
      call PrintErrMsg(option)
    endif
    call VecDestroy(surf_ids,ierr);CHKERRQ(ierr)

    select case(source_mesh_id)
      case(ELM_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, grid_elm_npts_local, &
                                         grid_elm_cell_ids_nindex_copy)
        call MappingSetDestinationMeshCellIds(map, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SUB_MESH)
        call MappingSetSourceMeshCellIds(map, grid_pf_npts_local, &
                                        grid_pf_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, grid_elm_npts_local, &
                                              grid_elm_npts_ghost, &
                                              grid_elm_cell_ids_nindex_copy, &
                                              grid_elm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to ' // &
          'pflotranModelInitMappingSurf3D'
        call PrintErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_elm_local_nindex)

    call MappingDecompose(map, option%mycomm)
    call MappingFindDistinctSourceMeshCellIds(map)
    call MappingCreateWeightMatrix(map, option%myrank)
    call MappingCreateScatterOfSourceMesh(map, option%mycomm)

  end subroutine pflotranModelInitMapSrfTo2DSub

! ************************************************************************** !

  subroutine pflotranModelStepperRunTillPauseTime(model, pause_time)
  !
  ! It performs the model integration
  ! till the specified pause_time.
  ! NOTE: It is assumed 'pause_time' is in seconds.
  ! NOTE(bja, 2013-07) the strange waypoint insertion of t+30min /
  ! deletion of t is to ensure that we always have a valid waypoint in
  ! front of us, but pflotran does not delete them, so we don't want
  ! to accumulate too large of a list.
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !


    implicit none

    type(pflotran_model_type), pointer :: model
    PetscReal, intent(in) :: pause_time

    PetscReal :: pause_time1

    if (OptionIsIORank(model%option)) then
       write(model%option%fid_out, *), '>>>> Inserting waypoint at pause_time = ', pause_time
    endif
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][pflotran_model_mod::pflotranModelStepperRunTillPauseTime] pause_time = ', pause_time
#endif
    pause_time1 = pause_time + 1800.0d0
    call pflotranModelInsertWaypoint(model, pause_time1)
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][pflotran_model_mod::pflotranModelStepperRunTillPauseTime] pause_time1 = ', pause_time1
#endif
    call model%simulation%RunToTime(pause_time)

    call pflotranModelDeleteWaypoint(model, pause_time)

    ! TODO(GB): Use XXXUpdateMassBalancePatch() to ensure mass balance
    ! betweent ELM calls

  end subroutine pflotranModelStepperRunTillPauseTime

! ************************************************************************** !

  subroutine pflotranModelInsertWaypoint(model, waypoint_time)
  !
  ! Inserts a waypoint within the waypoint list
  ! so that the model integration can be paused when that waypoint is
  ! reached
  ! NOTE: It is assumed the 'waypoint_time' is in seconds
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Subsurface_class, only : realization_subsurface_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointInsertInList
    use Units_module, only : UnitsConvertToInternal
    use Option_module

    implicit none

    type(pflotran_model_type), pointer :: model
    type(waypoint_type), pointer       :: waypoint
    type(waypoint_type), pointer       :: waypoint2
    PetscReal                          :: waypoint_time
    character(len=MAXWORDLENGTH)       :: word

    class(realization_subsurface_type), pointer    :: realization

    word = 's'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time !* UnitsConvertToInternal(word,'sec',model%option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 3153600.d0
    waypoint2 => WaypointCreate(waypoint)

    select type (simulation => model%simulation)
      class is (simulation_subsurface_type)
         call WaypointInsertInList(waypoint, simulation%waypoint_list_subsurface)
      class default
         nullify(realization)
         model%option%io_buffer = "pflotranModelInsertWaypoint only " // &
              "works on subsurface simulations."
         call PrintErrMsg(model%option)
    end select

  end subroutine pflotranModelInsertWaypoint

! ************************************************************************** !

  subroutine pflotranModelDeleteWaypoint(model, waypoint_time)

    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    use Realization_Subsurface_class, only : realization_subsurface_type

    use Waypoint_module, only : waypoint_type, WaypointCreate, WaypointDeleteFromList
    use Units_module, only : UnitsConvertToInternal
    use Option_module

    implicit none

    type(pflotran_model_type), pointer :: model
    type(waypoint_type), pointer       :: waypoint
    PetscReal                          :: waypoint_time
    character(len=MAXWORDLENGTH)       :: word
    character(len=MAXWORDLENGTH)       :: internal_units

    class(realization_subsurface_type), pointer    :: realization

    word = 's'
    waypoint => WaypointCreate()
    waypoint%time              = waypoint_time !* UnitsConvertToInternal(word,model%option)
    waypoint%update_conditions = PETSC_TRUE
    waypoint%dt_max            = 3153600

    select type (simulation => model%simulation)
      class is (simulation_subsurface_type)
         call WaypointDeleteFromList(waypoint, simulation%waypoint_list_subsurface)
      class default
         nullify(realization)
         model%option%io_buffer = "pflotranModelInsertWaypoint only " // &
              "works on subsurface simulations."
         call PrintErrMsg(model%option)
    end select

  end subroutine pflotranModelDeleteWaypoint

! ************************************************************************** !

  subroutine pflotranModelSetupRestart(model, restart_stamp)
  !
  ! pflotranModelSetupRestart()
  ! This checks to see if a restart file stamp was provided by the
  ! driver. If so, we set the restart flag and reconstruct the
  ! restart file name. The actual restart is handled by the standard
  ! pflotran mechanism in TimeStepperInitializeRun()
  ! NOTE: this must be called between pflotranModelCreate() and
  ! pflotranModelStepperRunInit()
  !

    use Option_module
    use String_module

    implicit none

    type(pflotran_model_type), pointer :: model
    character(len=MAXWORDLENGTH) :: restart_stamp

    call PrintWrnMsg(model%option)

    if (.not. StringNull(restart_stamp)) then
       model%option%restart_flag = PETSC_TRUE
       model%option%restart_filename = &
            trim(model%option%global_prefix) // &
            trim(model%option%group_prefix) // &
            '-' // trim(restart_stamp) // '.chk'
    end if

  end subroutine pflotranModelSetupRestart

! ************************************************************************** !

  subroutine pflotranModelUpdateSourceSink(pflotran_model)
  !
  ! Update the source/sink term
  !
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  !

    use elm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Grid_module
    use Global_Aux_module
    use Mapping_module
    use Option_module
    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

    class(realization_subsurface_type), pointer          :: subsurf_realization
    type(coupler_type), pointer               :: source_sink
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: qflx_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscInt                                  :: local_id
    PetscInt                                  :: ghosted_id
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_aux_vars(:)

    call MappingSourceToDestination(pflotran_model%map_elm_sub_to_pf_sub, &
                                    elm_pf_idata%qflx_elm, &
                                    elm_pf_idata%qflx_pf)

    ! Get pointer to subsurface-realization
    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         subsurf_realization => simulation%realization
      class default
         pflotran_model%option%io_buffer = " Unsupported simulation_type " // &
            " in pflotranModelUpdateSourceSink."
         call PrintErrMsg(pflotran_model%option)
    end select

    global_aux_vars  => subsurf_realization%patch%aux%Global%auxvars
    grid             => subsurf_realization%patch%grid

    ! Find value of pressure-dof depending on flow mode
    select case (pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        press_dof = RICHARDS_PRESSURE_DOF
      case (TH_MODE)
        press_dof = TH_PRESSURE_DOF
      case default
        pflotran_model%option%io_buffer = 'Unsupported Flow mode'
        call PrintErrMsg(pflotran_model%option)
    end select

    ! Update the 'elm_et_ss' source/sink term
    call VecGetArrayF90(elm_pf_idata%qflx_pf,qflx_pf_loc,ierr);CHKERRQ(ierr)
    found = PETSC_FALSE
    source_sink => subsurf_realization%patch%source_sink_list%first
    do
      if (.not.associated(source_sink)) exit

      cur_connection_set => source_sink%connection_set

      ! Find appropriate Source/Sink from the list of Source/Sinks
      if(StringCompare(source_sink%name,'elm_et_ss')) then

        found = PETSC_TRUE
        if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) then
          call PrintErrMsg(pflotran_model%option,'elm_et_ss is not of ' // &
                           'HET_MASS_RATE_SS')
        endif

        do iconn = 1, cur_connection_set%num_connections
          source_sink%flow_aux_real_var(press_dof,iconn) = qflx_pf_loc(iconn)

          if (pflotran_model%option%iflowmode == TH_MODE) then
            local_id = cur_connection_set%id_dn(iconn)
            ghosted_id = grid%nL2G(local_id)
            source_sink%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
              global_aux_vars(ghosted_id)%temp
          end if
        enddo
      endif

      source_sink => source_sink%next
    enddo
    call VecRestoreArrayF90(elm_pf_idata%qflx_pf,qflx_pf_loc, &
                            ierr);CHKERRQ(ierr)

    if(.not.found) &
      call PrintErrMsg(pflotran_model%option,'elm_et_ss not found in ' // &
                       'source-sink list of subsurface model.')

  end subroutine pflotranModelUpdateSourceSink

! ************************************************************************** !

  subroutine pflotranModelUpdateFlowConds(pflotran_model)
  !
  ! This routine Updates boundary and source/sink condtions for PFLOTRAN that
  ! are prescribed by ELM
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  !

    use elm_pflotran_interface_data
    use Mapping_module
    use Option_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

    call pflotranModelUpdateSourceSink(pflotran_model)

  end subroutine pflotranModelUpdateFlowConds

! ************************************************************************** !

  subroutine pflotranModelUpdateSubsurfTCond(pflotran_model)
  !
  ! This routine updates subsurface boundary condtions of PFLOTRAN related to
  ! energy equation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/08/2013
  !

    use elm_pflotran_interface_data
    use Connection_module
    use Coupler_module
    use Mapping_module
    use Option_module
    use Realization_Subsurface_class, only : realization_subsurface_type
    use Simulation_Base_class, only : simulation_base_type
    use String_module
    use Simulation_Subsurface_class, only : simulation_subsurface_type

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model

    class(realization_subsurface_type), pointer          :: subsurf_realization
    type(coupler_type), pointer               :: boundary_condition
    type(connection_set_type), pointer        :: cur_connection_set
    PetscScalar, pointer                      :: gflux_subsurf_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof

    ! Map ground-heat flux from ELM--to--PF grid
    call MappingSourceToDestination(pflotran_model%map_elm_srf_to_pf_2dsub, &
                                    elm_pf_idata%gflux_subsurf_elm, &
                                    elm_pf_idata%gflux_subsurf_pf)

    ! Get pointer to subsurface-realization
    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         subsurf_realization => simulation%realization
      class default
         pflotran_model%option%io_buffer = " Unsupported simulation_type " // &
            " in pflotranModelUpdateSubSurfTCond."
         call PrintErrMsg(pflotran_model%option)
    end select

    ! Update the 'elm_et_ss' source/sink term
    call VecGetArrayF90(elm_pf_idata%gflux_subsurf_pf,gflux_subsurf_pf_loc, &
                        ierr);CHKERRQ(ierr)
    found = PETSC_FALSE
    boundary_condition => subsurf_realization%patch%boundary_condition_list%first
    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      ! Find appropriate BC from the list of boundary conditions
      if(StringCompare(boundary_condition%name,'elm_gflux_bc')) then

        if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) &
            /= NEUMANN_BC) then
          call PrintErrMsg(pflotran_model%option,'elm_gflux_bc is not of ' // &
                           'NEUMANN_BC')
        endif
        found = PETSC_TRUE

        do iconn = 1, cur_connection_set%num_connections
          boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
            gflux_subsurf_pf_loc(iconn)*pflotran_model%option%scale
        enddo
      endif

      boundary_condition => boundary_condition%next
    enddo
    call VecRestoreArrayF90(elm_pf_idata%gflux_subsurf_pf, &
                            gflux_subsurf_pf_loc,ierr);CHKERRQ(ierr)

    if(.not.found) &
      call PrintErrMsg(pflotran_model%option,'elm_gflux_bc not found in ' // &
                       'source-sink list of subsurface model.')

  end subroutine pflotranModelUpdateSubsurfTCond

! ************************************************************************** !

  subroutine pflotranModelGetUpdatedData(pflotran_model)
  !
  ! This routine get updated states evoloved by PFLOTRAN.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 5/14/2013
  !

    use Option_module
    use Richards_module
    use Richards_Aux_module
    use TH_module
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use elm_pflotran_interface_data
    use Realization_Subsurface_class, only : realization_subsurface_type

    type(pflotran_model_type), pointer  :: pflotran_model
    class(realization_subsurface_type), pointer     :: realization

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select

    select case(pflotran_model%option%iflowmode)
      case (RICHARDS_MODE)
        call RichardsUpdateAuxVars(realization)
        call pflotranModelGetSaturation(pflotran_model)
      case (TH_MODE)
        call THUpdateAuxVars(realization)
        call pflotranModelGetSaturation(pflotran_model)
        call pflotranModelGetTemperature(pflotran_model)
        call pflotranModelGetEffThermCond(pflotran_model)
      case default
        pflotran_model%option%io_buffer='pflotranModelGetSaturation ' // &
          'not implmented for this mode.'
        call PrintErrMsg(pflotran_model%option)
    end select

  end subroutine pflotranModelGetUpdatedData

! ************************************************************************** !

  subroutine pflotranModelGetSaturation(pflotran_model)
  !
  ! Extract soil saturation values simulated by
  ! PFLOTRAN in a PETSc vector.
  !
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  !

    use Option_module
    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Material_Aux_module, only : material_auxvar_type
    use elm_pflotran_interface_data
    use Mapping_module
    use TH_Aux_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: sat_pf_p(:)
    PetscReal, pointer :: mass_pf_p(:)
    type(TH_auxvar_type),pointer :: TH_auxvars(:)
    type(option_type), pointer :: option
    class(material_auxvar_type), pointer :: material_auxvars(:)

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    global_aux_vars => patch%aux%Global%auxvars
    option          => realization%option
    material_auxvars=> patch%aux%Material%auxvars
! #ifdef DEBUG_ELMPFEH
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] before check mass_pf and mass_elm '
!      !stop
! #endif
    ! Save the saturation values
    call VecGetArrayF90(elm_pf_idata%sat_pf,sat_pf_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%mass_pf,mass_pf_p, ierr);CHKERRQ(ierr)
    do local_id=1, grid%nlmax
      ghosted_id=grid%nL2G(local_id)
      sat_pf_p(local_id)=global_aux_vars(ghosted_id)%sat(1)
      mass_pf_p(local_id)= &
      global_aux_vars(ghosted_id)%sat(1) * &
      global_aux_vars(ghosted_id)%den_kg(1) * &
      material_auxvars(ghosted_id)%volume* &
      material_auxvars(ghosted_id)%porosity
! #ifdef DEBUG_ELMPFEH
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] check mass_pf and mass_elm '
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] local_id = ', local_id
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] ghosted_id = ', ghosted_id
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] sat_pf_p = ', sat_pf_p(local_id)
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] mass_pf_p = ', mass_pf_p(local_id)
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation]     |- global_aux_vars(ghosted_id)%sat(1) = ', global_aux_vars(ghosted_id)%sat(1)
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation]     |- global_aux_vars(ghosted_id)%den_kg(1) = ', global_aux_vars(ghosted_id)%den_kg(1)
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation]     |- material_auxvars(ghosted_id)%volume = ', material_auxvars(ghosted_id)%volume
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation]     |- material_auxvars(ghosted_id)%porosity = ', material_auxvars(ghosted_id)%porosity
!      !stop
! #endif
    enddo
    call VecRestoreArrayF90(elm_pf_idata%sat_pf,sat_pf_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%mass_pf,mass_pf_p,ierr);CHKERRQ(ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%sat_pf, &
                                    elm_pf_idata%sat_elm)
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                     elm_pf_idata%mass_pf, &
                                     elm_pf_idata%mass_elm)

    if (pflotran_model%option%iflowmode == TH_MODE .and. &
        option%flow%th_freezing) then

      TH_auxvars => patch%aux%TH%auxvars

      call VecGetArrayF90(elm_pf_idata%sat_ice_pf,sat_pf_p, &
                          ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        sat_pf_p(local_id) = TH_auxvars(ghosted_id)%ice%sat_ice
      enddo
      call VecGetArrayF90(elm_pf_idata%sat_ice_pf,sat_pf_p, &
                          ierr);CHKERRQ(ierr)

      call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                      elm_pf_idata%sat_ice_pf, &
                                      elm_pf_idata%sat_ice_elm)
    endif

  end subroutine pflotranModelGetSaturation

  ! ************************************************************************** !

  subroutine pflotranModelGetInternalflow(pflotran_model)
  !
  ! Extract internal flow fluxes simulated by
  ! PFLOTRAN in a PETSc vector.
  !
  ! Author: Yi Xiao
  ! Date: 07/11/2024
  !

    use Option_module
    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Material_Aux_module, only : material_auxvar_type
    use elm_pflotran_interface_data
    use Mapping_module
    use TH_Aux_module
    use Connection_module
    use Output_HDF5_module
    use Output_Aux_module
    use Output_Common_module
    use hdf5
    use HDF5_module
    use HDF5_Aux_module
    use String_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscInt           :: temp_int, iconn, sum_connection
    PetscInt           :: skip_conn_type, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
    ! PetscReal, pointer ::
    type(TH_auxvar_type),pointer :: TH_auxvars(:)
    type(option_type), pointer :: option
    class(material_auxvar_type), pointer :: material_auxvars(:)
    type(output_option_type), pointer :: output_option
    type(connection_set_list_type), pointer :: connection_set_list
    type(connection_set_type), pointer :: cur_connection_set
    PetscReal, pointer :: temp_vertical_influx_p(:), temp_lateral_influx_p(:), temp_vertical_efflux_p(:), temp_lateral_efflux_p(:)
    PetscReal, parameter :: eps = 1.0e-5
    PetscReal :: temp_dirz
    PetscViewer :: viewer
    character(len=MAXSTRINGLENGTH) :: string
    character(len=MAXSTRINGLENGTH) :: filename
    character(len=MAXWORDLENGTH) :: word
    PetscInt :: var_list_type
    PetscBool :: hdf5_first
    integer(HID_T) :: file_id
    integer(HID_T) :: grp_id
    type(output_variable_type), pointer :: cur_variable
    Vec :: global_vec

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
         !write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] simulation%waypoint_list_subsurface%last%time = ', simulation%waypoint_list_subsurface%last%time
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    global_aux_vars => patch%aux%Global%auxvars
    option          => realization%option
    material_auxvars=> patch%aux%Material%auxvars
    output_option   => realization%output_option

    ! dev: try output connection set here; then see how to classify into lateral and vertical flow; then see how to map to ELM
#ifdef DEBUG_ELMPFEH
    temp_int = ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
    write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] associated(patch%internal_flow_fluxes) = ' , associated(patch%internal_flow_fluxes)
    if (associated(patch%internal_flow_fluxes)) then
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] patch%internal_flow_fluxes(1,:)=', patch%internal_flow_fluxes(1,:)
    endif
    write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_int=', temp_int
    !write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] patch%imat=', patch%imat

    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0
    do
      if (.not.associated(cur_connection_set)) exit

      call VecGetArrayF90(elm_pf_idata%vertical_influx_pf,temp_vertical_influx_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(elm_pf_idata%lateral_influx_pf,temp_lateral_influx_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(elm_pf_idata%vertical_efflux_pf,temp_vertical_efflux_p,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(elm_pf_idata%lateral_efflux_pf,temp_lateral_efflux_p,ierr);CHKERRQ(ierr)

      temp_vertical_influx_p = 0.0
      temp_lateral_influx_p = 0.0
      temp_vertical_efflux_p = 0.0
      temp_lateral_efflux_p = 0.0

      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] initial temp_vertical_influx_p = ', temp_vertical_influx_p
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] initial temp_lateral_influx_p = ', temp_lateral_influx_p
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] initial temp_vertical_efflux_p = ', temp_vertical_efflux_p
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] initial temp_lateral_efflux_p = ', temp_lateral_efflux_p


      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        ghosted_id_up = cur_connection_set%id_up(iconn)
        ghosted_id_dn = cur_connection_set%id_dn(iconn)

        local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

        if (patch%imat(ghosted_id_up) <= 0 .or.  &
            patch%imat(ghosted_id_dn) <= 0) cycle

        ! if (.not.(skip_conn_type == NO_CONN)) then
        !   if (skip_conn(cur_connection_set%dist(1:3,iconn), skip_conn_type)) cycle
        ! endif

        if (associated(patch%internal_flow_fluxes)) then
          temp_dirz = abs(cur_connection_set%dist(3,iconn))
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_dirz = ', temp_dirz
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] local_id_up = ', local_id_up
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] local_id_dn = ', local_id_dn
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] patch%internal_flow_fluxes(1,sum_connection) = ', patch%internal_flow_fluxes(1,sum_connection)
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_vertical_efflux_p(local_id_up) = ', temp_vertical_efflux_p(local_id_up)
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_vertical_influx_p(local_id_dn) = ', temp_vertical_influx_p(local_id_dn)
          if (abs(temp_dirz - 1.0)<eps) then
            ! vertical flow
            if (local_id_up>0) then
              temp_vertical_efflux_p(local_id_up) = temp_vertical_efflux_p(local_id_up) + patch%internal_flow_fluxes(1,sum_connection)
            endif
            if (local_id_dn>0) then
              temp_vertical_influx_p(local_id_dn) = temp_vertical_influx_p(local_id_dn) + patch%internal_flow_fluxes(1,sum_connection)
            endif
          else
            ! "lateral" flow cross column
            if (local_id_up>0) then
              temp_lateral_efflux_p(local_id_up) = temp_lateral_efflux_p(local_id_up) + patch%internal_flow_fluxes(1,sum_connection)
            endif
            if (local_id_dn>0) then
              temp_lateral_influx_p(local_id_dn) = temp_lateral_influx_p(local_id_dn) + patch%internal_flow_fluxes(1,sum_connection)
            endif
          endif
        endif

        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] iconn = ', iconn
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] sum_connection = ', sum_connection
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- ghosted_id_up = ', ghosted_id_up
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- ghosted_id_dn = ', ghosted_id_dn
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- local_id_up = ', local_id_up
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- local_id_dn = ', local_id_dn
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- patch%internal_flow_fluxes(1,sum_connection) = ', patch%internal_flow_fluxes(1,sum_connection)
        ! !write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- patch%internal_velocities(1,sum_connection) = ', patch%internal_velocities(1,sum_connection)
        ! !write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- cur_connection_set%itype = ', cur_connection_set%itype
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- cur_connection_set%dist(1:3,iconn) = ', cur_connection_set%dist(1:3,iconn)
        ! write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- abs(temp_dirz - 1.0)<eps = ', abs(temp_dirz - 1.0)<eps
        ! !write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- istart up = ', (local_id_up-1)*option%nflowdof + 1
        ! !write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow]   |- istart dn = ', (local_id_dn-1)*option%nflowdof + 1

      end do

      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_vertical_influx_p = ', temp_vertical_influx_p
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_lateral_influx_p = ', temp_lateral_influx_p
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_vertical_efflux_p = ', temp_vertical_efflux_p
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] temp_lateral_efflux_p = ', temp_lateral_efflux_p

      call VecRestoreArrayF90(elm_pf_idata%vertical_influx_pf,temp_vertical_influx_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(elm_pf_idata%lateral_influx_pf,temp_lateral_influx_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(elm_pf_idata%vertical_efflux_pf,temp_vertical_efflux_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(elm_pf_idata%lateral_efflux_pf,temp_lateral_efflux_p,ierr);CHKERRQ(ierr)

#ifdef PRINT_INTERNALFLOW
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] print elm_pf_idata%vertial/lateral_influx/efflux_pf to pf_internalflow.out'
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] option%time = ', option%time
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] realization%output_option%tconv = ', realization%output_option%tconv

      ! modified based on output_hdf9.F90, subroutine OutputHDF5
      !call OutputHDF5OpenFile(option, output_option, var_list_type, file_id, first)
      if (abs(option%time - 1800.0) < eps) then
        hdf5_first = PETSC_TRUE
      else
        hdf5_first = PETSC_FALSE
      endif
      filename = 'pf_internalflow.h5'
      if (.not.hdf5_first) then
        call HDF5FileTryOpen(filename,file_id,hdf5_first,option%comm)
      endif
      if (hdf5_first) then
        call HDF5FileOpen(filename,file_id,PETSC_TRUE,option)
      endif

      ! create a group for the data set
      write(string,'(''Time:'',es13.5,x,a1)') &
            option%time/output_option%tconv,output_option%tunit
      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetInternalflow] string = ', string
      call HDF5GroupOpenOrCreate(file_id,string,grp_id,option)

      ! write group attributes
      call OutputHDF5WriteSnapShotAtts(grp_id,option)

      ! write elm_pf_idata%{vertical_influx_pf, lateral_influx_pf, vertical_efflux_pf, lateral_efflux_pf} to file
      ! write elm_pf_idata%{mass_elm, area_top_face_elm, qflx_elm} to file
      string = 'vertical_influx_pf'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%vertical_influx_pf,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'lateral_influx_pf'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%lateral_influx_pf,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'vertical_efflux_pf'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%vertical_efflux_pf,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'lateral_efflux_pf'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%lateral_efflux_pf,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'area_top_face_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%area_top_face_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'area_top_face_pf'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%area_top_face_pf,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mflx_infl_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mflx_infl_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mflx_et_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mflx_et_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mflx_dew_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mflx_dew_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mflx_sub_snow_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mflx_sub_snow_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mflx_snowlyr_disp_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mflx_snowlyr_disp_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mflx_drain_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mflx_drain_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      string = 'mass_elm'
      call HDF5WriteDataSetFromVec(string,option,elm_pf_idata%mass_elm,grp_id, &
                                    H5T_NATIVE_DOUBLE)
      call HDF5GroupClose(grp_id,option)

      call OutputHDF5CloseFile(option, file_id)

      hdf5_first = PETSC_FALSE

#endif

      cur_connection_set => cur_connection_set%next
    end do
    !stop
#endif

!     do local_id=1, grid%nlmax
!       ghosted_id=grid%nL2G(local_id)
!       sat_pf_p(local_id)=global_aux_vars(ghosted_id)%sat(1)
!       mass_pf_p(local_id)= &
!       global_aux_vars(ghosted_id)%sat(1) * &
!       global_aux_vars(ghosted_id)%den_kg(1) * &
!       material_auxvars(ghosted_id)%volume* &
!       material_auxvars(ghosted_id)%porosity
! ! #ifdef DEBUG_ELMPFEH
! !      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] check mass_pf and mass_elm '
! !      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] local_id = ', local_id
! !      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] ghosted_id = ', ghosted_id
! !      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetSaturation] sat_pf_p = ', sat_pf_p(local_id)
! !      !stop
! ! #endif
!     enddo
!     call VecRestoreArrayF90(elm_pf_idata%internalflow_pf,internalflow_pf_p,ierr);CHKERRQ(ierr)

!     call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
!                                     elm_pf_idata%internalflow_pf, &
!                                     elm_pf_idata%internalflow_elm)

!     ! if (pflotran_model%option%iflowmode == TH_MODE .and. &
!     !     option%flow%th_freezing) then
!     ! ! special treatment for TH_MODE?
!     ! endif

  end subroutine pflotranModelGetInternalflow
! ************************************************************************** !

subroutine OutputHDF5WriteSnapShotAtts(parent_id,option)
  !
  ! Writes attributes associated with a snapshot time in the output file.
  !
  ! Author: Glenn Hammond
  ! Date: 07/31/19
  !
  use hdf5
  use Option_module

  implicit none

  integer(HID_T) :: parent_id
  type(option_type) :: option

  integer(HID_T) :: attribute_id
  integer(HID_T) :: dataspace_id
  character(len=MAXSTRINGLENGTH) :: string
  integer(HSIZE_T) :: dims(1)
  PetscMPIInt :: hdf5_err
  PetscMPIInt, parameter :: ON=1, OFF=0

  dims = 1
  call h5screate_simple_f(1,dims,dataspace_id,hdf5_err)
  string = 'Time (s)'
  call h5eset_auto_f(OFF,hdf5_err)
  call h5aopen_f(parent_id, string, attribute_id, hdf5_err)
  if (hdf5_err /= 0) then
    call h5acreate_f(parent_id,string,H5T_NATIVE_DOUBLE,dataspace_id, &
                     attribute_id,hdf5_err)
  endif
  call h5eset_auto_f(ON,hdf5_err)
  call h5awrite_f(attribute_id,H5T_NATIVE_DOUBLE,option%time,dims,hdf5_err)
  call h5aclose_f(attribute_id, hdf5_err)
  call h5sclose_f(dataspace_id, hdf5_err)

end subroutine OutputHDF5WriteSnapShotAtts
! ************************************************************************** !

  subroutine pflotranModelGetTemperature(pflotran_model)
  !
  ! This routine get updated states evoloved by PFLOTRAN.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 5/14/2013
  !

    use Option_module
    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use elm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer           :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(global_auxvar_type), pointer         :: global_aux_vars(:)
    type(th_auxvar_type), pointer             :: th_aux_vars(:)

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscScalar, pointer :: temp_pf_p(:)
    PetscReal, pointer :: sat_ice_pf_p(:)

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    global_aux_vars => patch%aux%Global%auxvars
    th_aux_vars     => patch%aux%TH%auxvars

    call VecGetArrayF90(elm_pf_idata%temp_pf,temp_pf_p,ierr);CHKERRQ(ierr)
    do ghosted_id=1,grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id>0) then
        !call VecSetValues(elm_pf_idata%sat_ice_pf,1,local_id-1, &
        !                 global_aux_vars(ghosted_id)%sat_ice(1),INSERT_VALUES,ierr)
        temp_pf_p(local_id) = global_aux_vars(ghosted_id)%temp
      endif
    enddo
    call VecRestoreArrayF90(elm_pf_idata%temp_pf,temp_pf_p, &
                            ierr);CHKERRQ(ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%temp_pf, &
                                    elm_pf_idata%temp_elm)

!    call VecAssemblyBegin(elm_pf_idata%sat_ice_pf,ierr)
!    call VecAssemblyEnd(elm_pf_idata%sat_ice_pf,ierr)
!    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
!                                    pflotran_model%option, &
!                                    elm_pf_idata%sat_ice_pf, &
!                                    elm_pf_idata%sat_ice_elm)

  end subroutine pflotranModelGetTemperature

! ************************************************************************** !

  subroutine pflotranModelGetEffThermCond(pflotran_model)
  !
  ! This routine get updated effective soil thermal conductivity value.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 9/18/2014
  !

    use Option_module
    use Realization_Subsurface_class
    use Patch_module
    use Grid_module
    use Global_Aux_module
    use TH_Aux_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use elm_pflotran_interface_data
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer        :: pflotran_model
    class(realization_subsurface_type), pointer          :: realization
    type(patch_type), pointer                 :: patch
    type(grid_type), pointer                  :: grid
    type(th_auxvar_type), pointer             :: th_aux_vars(:)

    PetscErrorCode       :: ierr
    PetscInt             :: local_id, ghosted_id
    PetscScalar, pointer :: eff_tc_p(:)

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    th_aux_vars     => patch%aux%TH%auxvars

    call VecGetArrayF90(elm_pf_idata%eff_therm_cond_pf,eff_tc_p, &
                        ierr);CHKERRQ(ierr)
    do ghosted_id=1,grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      if (local_id > 0) then
        ! Converting MJ/m/K to J/m/K
        eff_tc_p(local_id) = th_aux_vars(ghosted_id)%Dk_eff/pflotran_model%option%scale
      endif
    enddo
    call VecRestoreArrayF90(elm_pf_idata%eff_therm_cond_pf,eff_tc_p, &
                            ierr);CHKERRQ(ierr)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%eff_therm_cond_pf, &
                                    elm_pf_idata%eff_therm_cond_elm)

  end subroutine pflotranModelGetEffThermCond

! ************************************************************************** !

  subroutine pflotranModelStepperRunFinalize(model)
  !
  ! It performs the same execution of commands
  ! that are carried out in StepperRun() once the model integration is
  ! finished
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    implicit none

    type(pflotran_model_type), pointer :: model

    call model%simulation%FinalizeRun()

  end subroutine pflotranModelStepperRunFinalize

! ************************************************************************** !

  function pflotranModelNSurfCells3DDomain(pflotran_model)
  !
  ! This function returns the number of control volumes forming surface of
  ! the sub-surface domain. The subroutines assumes the following boundary
  ! condition is specified in inputdeck:
  ! - 'elm_gflux_bc': when running subsurface only simulation.
  ! - 'from_surface_bc': when running surface-subsurface simulation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 6/03/2013
  !

    use Option_module
    use Coupler_module
    use String_module
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model

    class(realization_subsurface_type), pointer :: realization
    type(coupler_list_type), pointer :: coupler_list
    type(coupler_type), pointer :: coupler
    type(simulation_base_type), pointer :: simulation
    character(len=MAXWORDLENGTH) :: condition_name
    PetscInt :: pflotranModelNSurfCells3DDomain
    PetscBool :: found

    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select

    ! Determine the BC coupler name to search from list of BCs depending on
    ! subsurface or surface-subsurface simulation.
    condition_name = 'elm_gflux_bc'

    coupler_list => realization%patch%boundary_condition_list
    coupler => coupler_list%first
    found = PETSC_FALSE

    do
      if (.not.associated(coupler)) exit
      if (StringCompare(coupler%name,trim(condition_name))) then
        pflotranModelNSurfCells3DDomain=coupler%connection_set%num_connections
        found = PETSC_TRUE
      endif
      coupler => coupler%next
    enddo

    if(.not.found)  &
      call PrintErrMsg(pflotran_model%option, &
            'Missing from the input deck a BC named elm_gflux_bc')

  end function pflotranModelNSurfCells3DDomain

! ************************************************************************** !

  subroutine pflotranModelGetTopFaceArea(pflotran_model)
  !
  ! This subroutine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 6/10/2013
  !

    use Option_module
    use Patch_module
    use Discretization_module
    use Grid_Unstructured_Aux_module
    use Grid_Unstructured_Cell_module
    use Grid_Unstructured_module
    use Grid_module
    use elm_pflotran_interface_data
    use Utility_module, only : DotProduct, CrossProduct
    use Simulation_Base_class, only : simulation_base_type
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type
    use Mapping_module

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model

    type(option_type), pointer :: option
    class(realization_subsurface_type), pointer :: realization
    type(discretization_type), pointer :: discretization
    type(patch_type), pointer :: patch
    type(grid_type), pointer :: grid
    !type(point_type) :: point1, point2, point3, point4

    PetscInt :: local_id
    PetscInt :: ghosted_id
    PetscInt :: iface
    PetscInt :: face_id
    PetscInt :: cell_type
    PetscInt :: vertex_ids(4)

    PetscReal :: area1

    PetscScalar, pointer :: area_p(:)
    PetscErrorCode :: ierr

    option => pflotran_model%option
    select type (simulation => pflotran_model%simulation)
      type is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         nullify(realization)
         pflotran_model%option%io_buffer = "ERROR: XXX only works on subsurface simulations."
         call PrintErrMsg(pflotran_model%option)
    end select
    discretization => realization%discretization
    patch => realization%patch
    grid => patch%grid
! #ifdef DEBUG_ELMPFEH
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetTopFaceArea] grid%itype = ', grid%itype
!      !stop
! #endif
    call VecGetArrayF90(elm_pf_idata%area_top_face_pf,area_p, &
                        ierr);CHKERRQ(ierr)

    if (grid%itype == STRUCTURED_GRID) then
      ! Structured grid
      do ghosted_id=1,grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if(local_id>0) then
          area1 = grid%structured_grid%dx(ghosted_id)* &
                  grid%structured_grid%dy(ghosted_id)
          area_p(local_id) = area1
        endif
      enddo
    else if (grid%itype == UNSTRUCTURED_GRID .or. grid%itype == IMPLICIT_UNSTRUCTURED_GRID .or. grid%itype == EXPLICIT_UNSTRUCTURED_GRID) then
      ! Unstructured grid
      do local_id = 1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        cell_type = grid%unstructured_grid%cell_type(local_id)

        ! Find iface
        if (cell_type == HEX_TYPE) then
          iface = 6
        else if (cell_type == WEDGE_TYPE) then
          iface = 5
        else
          call PrintErrMsg(pflotran_model%option, &
            'Only hex and wedge cell_type supported in ELM-PFLOTRAN')
        endif
! #ifdef DEBUG_ELMPFEH
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetTopFaceArea] cell_type = ', cell_type
!      write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetTopFaceArea] iface = ', iface
!      !stop
! #endif
        ! Get face-id
        face_id = grid%unstructured_grid%cell_to_face_ghosted(iface, ghosted_id)

        ! Save face area
        area_p(local_id) = grid%unstructured_grid%face_area(face_id)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetTopFaceArea] face_id = ', face_id
     write(*,*) '[YX DEBUG][pflotran_model::pflotranModelGetTopFaceArea] area_p(local_id) = ', area_p(local_id)
     !stop
#endif
      enddo
    endif
    call VecRestoreArrayF90(elm_pf_idata%area_top_face_pf,area_p, &
                            ierr);CHKERRQ(ierr)
    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_elm_sub, &
                                    elm_pf_idata%area_top_face_pf, &
                                    elm_pf_idata%area_top_face_elm)
  end subroutine pflotranModelGetTopFaceArea

! ************************************************************************** !

  subroutine pflotranModelDestroy(model)
  !
  ! Deallocates the pflotranModel object
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  !

    use Factory_PFLOTRAN_module, only : FactoryPFLOTRANFinalize
    use Mapping_module, only : MappingDestroy
    use Communicator_Aux_module
    use Factory_Forward_module
    use Driver_class

    implicit none

    type(pflotran_model_type), pointer :: model
    PetscInt :: iflag
    type(comm_type), pointer :: comm
    class(driver_type), pointer :: driver

    ! FIXME(bja, 2013-07) none of the mapping information appears to
    ! be cleaned up, so we are leaking memory....

    driver => model%simulation%driver
    call model%simulation%FinalizeRun()
    call SimulationBaseDestroy(model%simulation)

    if (associated(model%map_elm_sub_to_pf_sub)) then
      call MappingDestroy(model%map_elm_sub_to_pf_sub)
      nullify(model%map_elm_sub_to_pf_sub)
    endif
    if (associated(model%map_elm_sub_to_pf_extended_sub)) then
      call MappingDestroy(model%map_elm_sub_to_pf_extended_sub)
      nullify(model%map_elm_sub_to_pf_extended_sub)
    endif
    if (associated(model%map_elm_srf_to_pf_2dsub)) then
      call MappingDestroy(model%map_elm_srf_to_pf_2dsub)
      nullify(model%map_elm_srf_to_pf_2dsub)
    endif
    if (associated(model%map_elm_srf_to_pf_srf)) then
      call MappingDestroy(model%map_elm_srf_to_pf_srf)
      nullify(model%map_elm_srf_to_pf_srf)
    endif
    if (associated(model%map_pf_sub_to_elm_sub)) then
      call MappingDestroy(model%map_pf_sub_to_elm_sub)
      nullify(model%map_pf_sub_to_elm_sub)
    endif
    if (associated(model%map_pf_srf_to_elm_srf)) then
      call MappingDestroy(model%map_pf_srf_to_elm_srf)
      nullify(model%map_pf_srf_to_elm_srf)
    endif

    call FactoryPFLOTRANFinalize(driver)
    iflag = driver%exit_code
    call DriverDestroy(driver)
    call exit(iflag)

  end subroutine pflotranModelDestroy

end module pflotran_model_module

