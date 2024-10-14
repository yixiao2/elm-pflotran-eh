module elm_pflotran_interface_data
!
! NOTES
!     (1) For data transfer from ELM to PFLOTRAN <- create mpi vec for ELM + seq vec for PFLOTRAN
!         For data transfer from PFLOTRAN to ELM <- create mpi vec for PFLOTRAN + seq vec for ELM
!         '*_elmp': mpi vec for ELM variables; '*_pfp': mpi vec for PFLOTRAN variables
!         '*_elms': seq vec for ELM variables; '*_pfs': seq vec for PFLOTRAN variables
!     (2) '*_subsurf_':
!     (3) '*_srf_':


#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
  use petscsys
  use petscvec

  implicit none

  private

  type, public :: elm_pflotran_idata_type

  !******************* A few global constants *******************************
  PetscInt :: nzelm_mapped ! num of ELM soil layers that are mapped

  !******************* Mesh property *******************************
  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlelm_sub ! num of local elm cells
  PetscInt :: ngelm_sub ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub  ! num of local pflotran cells
  PetscInt :: ngpf_sub  ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the surface of the 3D subsurface domain
  PetscInt :: nlelm_2dsub  ! num of local elm cells
  PetscInt :: ngelm_2dsub  ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dsub   ! num of local pflotran cells
  PetscInt :: ngpf_2dsub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlelm_srf  ! num of local elm cells
  PetscInt :: ngelm_srf  ! num of ghosted elm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf   ! num of local pflotran cells
  PetscInt :: ngpf_srf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Area of top face of PF 3D domain
  ! Vec :: area_top_face_elmp ! mpi vec
  ! Vec :: area_top_face_pfs  ! seq vec
  Vec :: area_top_face_pfp  ! mpi vec
  Vec :: area_top_face_elms ! seq vec

  ! Area of top face of ELM domain (projected area of PFLOTRAN surface)
  ! for elm flux convert to pflotran source sink
  ! Vec :: area_proj_top_face_elmp ! mpi vec
  ! Vec :: area_proj_top_face_pfs  ! seq vec
  ! Vec :: area_proj_top_face_pfp  ! mpi vec
  Vec :: area_proj_top_face_elms ! seq vec, to replace area_top_face_elms used in ExternalModelPFLOTRAN.F90

  ! cell IDs
  Vec :: pfgrid_nG2A_pfs ! mpi vec

  ! Time invariant data:

  ! (i) Soil properties -
  ! Local for ELM  - mpi vectors
  Vec :: hksat_x_elm
  Vec :: hksat_y_elm
  Vec :: hksat_z_elm
  Vec :: sucsat_elm
  Vec :: watsat_elm
  Vec :: bsw_elm
  Vec :: hksat_x2_elm
  Vec :: hksat_y2_elm
  Vec :: hksat_z2_elm
  Vec :: sucsat2_elm
  Vec :: watsat2_elm
  Vec :: bsw2_elm
  Vec :: thetares2_elm
  Vec :: press_elm

  ! Local for PFLOTRAN - seq. vec
  Vec :: hksat_x_pf
  Vec :: hksat_y_pf
  Vec :: hksat_z_pf
  Vec :: sucsat_pf
  Vec :: watsat_pf
  Vec :: bsw_pf
  Vec :: hksat_x2_pf
  Vec :: hksat_y2_pf
  Vec :: hksat_z2_pf
  Vec :: sucsat2_pf
  Vec :: watsat2_pf
  Vec :: bsw2_pf
  Vec :: thetares2_pf
  Vec :: press_pf

  ! to output vertical and lateral internal fluxes - seq vec
  ! classify patch%internal_flow_fluxes in PFLOTRAN, in pflotran_model_module::pflotranModelGetInternalflow
  Vec :: internal_flow_flux_vertical_pfs
  Vec :: internal_flow_flux_lateral_pfs

  ! (ii) Mesh property



  ! Time variant data

  ! (i) Sink/Source of water for PFLOTRAN's 3D subsurface domain
  Vec :: qflx_elm   ! mpi vec
  Vec :: qflx_pf    ! seq vec

  ! (ii) Source of water and temperature of rain for PFLOTRAN's 2D surface domain
  Vec :: rain_elm   ! mpi vec
  Vec :: rain_pf    ! seq vec
  Vec :: rain_temp_elm ! mpi vec
  Vec :: rain_temp_pf  ! seq vec

  ! (iii) Ground heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top surface of the subsurface domain
  Vec :: gflux_subsurf_elm  ! mpi vec
  Vec :: gflux_subsurf_pf   ! seq vec
  !       When running PFLOTRAN surface-subsurface simulation, ground heat flux
  !       is a SS for PFLOTRAN's surface domain.
  !
  !       Note: ELM decomposes the domain across processors in a horizontal.
  !       Thus, nlelm_2dsub = nlelm_srf across all processors. Thus, there is
  !       no need for 'gflux_surf_elm'

  ! (iv) Saturation and mass
  Vec :: sat_elms    ! seq vec
  Vec :: sat_pfp     ! mpi vec
  Vec :: mass_elms    ! seq vec
  Vec :: mass_pfp     ! mpi vec

  ! (v) Subsurface temperature
  Vec :: temp_elm   ! seq vec
  Vec :: temp_pf    ! mpi vec

  ! (vi) Ice saturation
  Vec :: sat_ice_elm ! seq vec
  Vec :: sat_ice_pf  ! mpi vec

  ! (vii) Stand water head
  Vec :: h2osfc_elm ! seq vec
  Vec :: h2osfc_pf  ! mpi vec

  Vec :: eff_therm_cond_elm ! seq vec
  Vec :: eff_therm_cond_pf  ! mpi vec


#ifdef PRINT_INTERNALFLOW
  ! to output component of qflx_elm; including mflx_infl/dew/et/sub_snow/snowlry_disp/drain
  Vec :: mflx_infl_elms
  Vec :: mflx_et_elms
  Vec :: mflx_dew_elms
  Vec :: mflx_sub_snow_elms
  Vec :: mflx_snowlyr_disp_elms
  Vec :: mflx_drain_elms
#endif

  end type elm_pflotran_idata_type

  type(elm_pflotran_idata_type) , public, target , save :: elm_pf_idata

  public :: ELMPFLOTRANIDataInit, &
            ELMPFLOTRANIDataCreateVec, &
            ELMPFLOTRANIDataDestroy

contains

! ************************************************************************** !

  subroutine ELMPFLOTRANIDataInit()
  !
  ! This routine initialized the data transfer type.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  !

    implicit none

    elm_pf_idata%nlelm_sub = 0
    elm_pf_idata%ngelm_sub = 0
    elm_pf_idata%nlpf_sub = 0
    elm_pf_idata%ngpf_sub = 0

    elm_pf_idata%nlelm_2dsub = 0
    elm_pf_idata%ngelm_2dsub = 0
    elm_pf_idata%nlpf_2dsub = 0
    elm_pf_idata%ngpf_2dsub = 0

    elm_pf_idata%nlelm_srf = 0
    elm_pf_idata%ngelm_srf = 0
    elm_pf_idata%nlpf_srf = 0
    elm_pf_idata%ngpf_srf = 0

    elm_pf_idata%hksat_x_elm = PETSC_NULL_VEC
    elm_pf_idata%hksat_y_elm = PETSC_NULL_VEC
    elm_pf_idata%hksat_z_elm = PETSC_NULL_VEC
    elm_pf_idata%sucsat_elm = PETSC_NULL_VEC
    elm_pf_idata%watsat_elm = PETSC_NULL_VEC
    elm_pf_idata%bsw_elm = PETSC_NULL_VEC
    elm_pf_idata%hksat_x2_elm = PETSC_NULL_VEC
    elm_pf_idata%hksat_y2_elm = PETSC_NULL_VEC
    elm_pf_idata%hksat_z2_elm = PETSC_NULL_VEC
    elm_pf_idata%sucsat2_elm = PETSC_NULL_VEC
    elm_pf_idata%watsat2_elm = PETSC_NULL_VEC
    elm_pf_idata%bsw2_elm = PETSC_NULL_VEC
    elm_pf_idata%thetares2_elm = PETSC_NULL_VEC
    elm_pf_idata%press_elm = PETSC_NULL_VEC

    elm_pf_idata%hksat_x_pf = PETSC_NULL_VEC
    elm_pf_idata%hksat_y_pf = PETSC_NULL_VEC
    elm_pf_idata%hksat_z_pf = PETSC_NULL_VEC
    elm_pf_idata%sucsat_pf = PETSC_NULL_VEC
    elm_pf_idata%watsat_pf = PETSC_NULL_VEC
    elm_pf_idata%bsw_pf = PETSC_NULL_VEC
    elm_pf_idata%hksat_x2_pf = PETSC_NULL_VEC
    elm_pf_idata%hksat_y2_pf = PETSC_NULL_VEC
    elm_pf_idata%hksat_z2_pf = PETSC_NULL_VEC
    elm_pf_idata%sucsat2_pf = PETSC_NULL_VEC
    elm_pf_idata%watsat2_pf = PETSC_NULL_VEC
    elm_pf_idata%bsw2_pf = PETSC_NULL_VEC
    elm_pf_idata%thetares2_pf = PETSC_NULL_VEC
    elm_pf_idata%press_pf = PETSC_NULL_VEC

    elm_pf_idata%internal_flow_flux_vertical_pfs = PETSC_NULL_VEC
    elm_pf_idata%internal_flow_flux_lateral_pfs = PETSC_NULL_VEC

    elm_pf_idata%qflx_elm = PETSC_NULL_VEC
    elm_pf_idata%qflx_pf = PETSC_NULL_VEC

    elm_pf_idata%rain_elm = PETSC_NULL_VEC
    elm_pf_idata%rain_pf = PETSC_NULL_VEC
    elm_pf_idata%rain_temp_elm = PETSC_NULL_VEC
    elm_pf_idata%rain_temp_pf = PETSC_NULL_VEC

    elm_pf_idata%gflux_subsurf_elm = PETSC_NULL_VEC
    elm_pf_idata%gflux_subsurf_pf = PETSC_NULL_VEC

    elm_pf_idata%sat_elms = PETSC_NULL_VEC
    elm_pf_idata%sat_pfp = PETSC_NULL_VEC
    elm_pf_idata%mass_elms = PETSC_NULL_VEC
    elm_pf_idata%mass_pfp = PETSC_NULL_VEC

    elm_pf_idata%temp_elm = PETSC_NULL_VEC
    elm_pf_idata%temp_pf = PETSC_NULL_VEC

    elm_pf_idata%sat_ice_elm = PETSC_NULL_VEC
    elm_pf_idata%sat_ice_pf = PETSC_NULL_VEC

    elm_pf_idata%h2osfc_elm = PETSC_NULL_VEC
    elm_pf_idata%h2osfc_pf = PETSC_NULL_VEC

    elm_pf_idata%eff_therm_cond_elm = PETSC_NULL_VEC
    elm_pf_idata%eff_therm_cond_pf = PETSC_NULL_VEC

    elm_pf_idata%nzelm_mapped = 0

    elm_pf_idata%area_top_face_pfp  = PETSC_NULL_VEC
    elm_pf_idata%area_top_face_elms = PETSC_NULL_VEC
    elm_pf_idata%area_proj_top_face_elms = PETSC_NULL_VEC

    elm_pf_idata%pfgrid_nG2A_pfs = PETSC_NULL_VEC

#ifdef PRINT_INTERNALFLOW
    elm_pf_idata%mflx_infl_elms = PETSC_NULL_VEC
    elm_pf_idata%mflx_et_elms = PETSC_NULL_VEC
    elm_pf_idata%mflx_dew_elms = PETSC_NULL_VEC
    elm_pf_idata%mflx_sub_snow_elms = PETSC_NULL_VEC
    elm_pf_idata%mflx_snowlyr_disp_elms = PETSC_NULL_VEC
    elm_pf_idata%mflx_drain_elms = PETSC_NULL_VEC
#endif

  end subroutine ELMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine ELMPFLOTRANIDataCreateVec(mycomm)
  !
  ! This routine creates PETSc vectors required for data transfer between
  ! ELM and PFLOTRAN.
  !
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  !

    implicit none

    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank
    PetscReal      :: zero = 0.0d0
    Vec :: vec_test

    call MPI_Comm_rank(mycomm,rank,ierr);CHKERRQ(ierr)

    !
    ! For data transfer from ELM to PFLOTRAN
    !

    ! Create MPI Vectors for ELM
    call VecCreateMPI(mycomm,elm_pf_idata%nlelm_sub,PETSC_DECIDE, &
                      elm_pf_idata%hksat_x_elm,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%hksat_x_elm,0.d0,ierr);CHKERRQ(ierr)

    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%hksat_y_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%hksat_z_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%sucsat_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%watsat_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%bsw_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%press_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_elm,elm_pf_idata%qflx_elm, &
                      ierr);CHKERRQ(ierr)

    call VecCreateMPI(mycomm,elm_pf_idata%nlelm_2dsub,PETSC_DECIDE, &
                      elm_pf_idata%gflux_subsurf_elm,ierr);CHKERRQ(ierr)
    call VecCreateMPI(mycomm,elm_pf_idata%nlelm_srf,PETSC_DECIDE, &
                      elm_pf_idata%rain_elm,ierr);CHKERRQ(ierr)
    call VecCreateMPI(mycomm,elm_pf_idata%nlelm_srf,PETSC_DECIDE, &
                      elm_pf_idata%rain_temp_elm,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%gflux_subsurf_elm,0.d0,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%rain_elm,0.d0,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%rain_temp_elm,0.d0,ierr);CHKERRQ(ierr)

    ! Create Seq. Vectors for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngpf_sub, &
                      elm_pf_idata%hksat_x_pf,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%hksat_x_pf,0.d0,ierr);CHKERRQ(ierr)

    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%hksat_y_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%hksat_z_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%sucsat_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%watsat_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%bsw_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%press_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%hksat_x_pf,elm_pf_idata%qflx_pf, &
                      ierr);CHKERRQ(ierr)

    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngpf_2dsub, &
                      elm_pf_idata%gflux_subsurf_pf,ierr);CHKERRQ(ierr)
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngpf_srf, &
                      elm_pf_idata%rain_pf,ierr);CHKERRQ(ierr)
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngpf_srf, &
                      elm_pf_idata%rain_temp_pf,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%gflux_subsurf_pf,0.d0,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%rain_pf,0.d0,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%rain_temp_pf,0.d0,ierr);CHKERRQ(ierr)

    !
    ! For data transfer from PFLOTRAN to ELM
    !

    ! Create MPI Vectors for PFLOTRAN
    ! 3D Subsurface PFLOTRAN ---to--- 3D Subsurface ELM
    call VecCreateMPI(mycomm,elm_pf_idata%nlpf_sub,PETSC_DECIDE, &
                      elm_pf_idata%sat_pfp,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%sat_pfp,0.d0,ierr);CHKERRQ(ierr)

    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%mass_pfp, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%temp_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%sat_ice_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%area_top_face_pfp, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%eff_therm_cond_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%hksat_x2_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%hksat_y2_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%hksat_z2_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%sucsat2_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%watsat2_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%bsw2_pf, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_pfp,elm_pf_idata%thetares2_pf, &
                      ierr);CHKERRQ(ierr)

    call VecCreateMPI(PETSC_COMM_SELF,elm_pf_idata%nlpf_sub,PETSC_DECIDE, &
                      elm_pf_idata%pfgrid_nG2A_pfs,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%pfgrid_nG2A_pfs,0.d0,ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%pfgrid_nG2A_pfs,elm_pf_idata%internal_flow_flux_vertical_pfs, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%pfgrid_nG2A_pfs,elm_pf_idata%internal_flow_flux_lateral_pfs, &
                      ierr);CHKERRQ(ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface ELM
    call VecCreateMPI(mycomm,elm_pf_idata%nlpf_srf,PETSC_DECIDE,elm_pf_idata%h2osfc_pf, &
    ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%h2osfc_pf,0.d0,ierr);CHKERRQ(ierr)

    ! Create Seq. Vectors for ELM
    ! 3D Subsurface PFLOTRAN ---to--- 3D Subsurface ELM
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%ngelm_sub, &
                      elm_pf_idata%sat_elms,ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%sat_elms,0.d0,ierr);CHKERRQ(ierr)

    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mass_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%temp_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%sat_ice_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%area_top_face_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%eff_therm_cond_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%hksat_x2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%hksat_y2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%hksat_z2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%sucsat2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%watsat2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%bsw2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%thetares2_elm, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%area_proj_top_face_elms, &
                      ierr);CHKERRQ(ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface ELM
    call VecCreateSeq(PETSC_COMM_SELF,elm_pf_idata%nlelm_2dsub,elm_pf_idata%h2osfc_elm, &
    ierr);CHKERRQ(ierr)
    call VecSet(elm_pf_idata%h2osfc_elm,0.d0,ierr);CHKERRQ(ierr)

#ifdef PRINT_INTERNALFLOW
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mflx_et_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mflx_drain_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mflx_infl_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mflx_dew_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mflx_sub_snow_elms, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(elm_pf_idata%sat_elms,elm_pf_idata%mflx_snowlyr_disp_elms, &
                      ierr);CHKERRQ(ierr)
#endif

  end subroutine ELMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine ELMPFLOTRANIDataDestroy()
  !
  ! This routine destroys PETSc vectors that were created for data transfer.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  !

    implicit none

    PetscErrorCode :: ierr

    if(elm_pf_idata%hksat_x_elm       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_x_elm,ierr)
    if(elm_pf_idata%hksat_y_elm       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_y_elm,ierr)
    if(elm_pf_idata%hksat_z_elm       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_z_elm,ierr)
    if(elm_pf_idata%sucsat_elm        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sucsat_elm,ierr)
    if(elm_pf_idata%watsat_elm        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%watsat_elm,ierr)
    if(elm_pf_idata%bsw_elm           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%bsw_elm,ierr)
    if(elm_pf_idata%hksat_x2_elm      /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_x2_elm,ierr)
    if(elm_pf_idata%hksat_y2_elm      /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_y2_elm,ierr)
    if(elm_pf_idata%hksat_z2_elm      /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_z2_elm,ierr)
    if(elm_pf_idata%sucsat2_elm       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sucsat2_elm,ierr)
    if(elm_pf_idata%watsat2_elm       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%watsat2_elm,ierr)
    if(elm_pf_idata%bsw2_elm          /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%bsw2_elm,ierr)
    if(elm_pf_idata%thetares2_elm     /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%thetares2_elm,ierr)
    if(elm_pf_idata%press_elm         /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%press_elm,ierr)

    if(elm_pf_idata%hksat_x_pf        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_x_pf,ierr)
    if(elm_pf_idata%hksat_y_pf        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_y_pf,ierr)
    if(elm_pf_idata%hksat_z_pf        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_z_pf,ierr)
    if(elm_pf_idata%sucsat_pf         /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sucsat_pf,ierr)
    if(elm_pf_idata%watsat_pf         /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%watsat_pf,ierr)
    if(elm_pf_idata%bsw_pf            /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%bsw_pf,ierr)
    if(elm_pf_idata%hksat_x2_pf       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_x2_pf,ierr)
    if(elm_pf_idata%hksat_y2_pf       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_y2_pf,ierr)
    if(elm_pf_idata%hksat_z2_pf       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%hksat_z2_pf,ierr)
    if(elm_pf_idata%sucsat2_pf        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sucsat2_pf,ierr)
    if(elm_pf_idata%watsat2_pf        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%watsat2_pf,ierr)
    if(elm_pf_idata%bsw2_pf           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%bsw2_pf,ierr)
    if(elm_pf_idata%thetares2_pf      /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%thetares2_pf,ierr)
    if(elm_pf_idata%press_pf          /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%press_pf,ierr)

    if(elm_pf_idata%internal_flow_flux_vertical_pfs /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%internal_flow_flux_vertical_pfs,ierr)
    if(elm_pf_idata%internal_flow_flux_lateral_pfs  /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%internal_flow_flux_lateral_pfs,ierr)

    if(elm_pf_idata%qflx_elm          /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%qflx_elm,ierr)
    if(elm_pf_idata%qflx_pf           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%qflx_pf,ierr)

    if(elm_pf_idata%rain_elm          /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%rain_elm,ierr)
    if(elm_pf_idata%rain_pf           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%rain_pf,ierr)
    if(elm_pf_idata%rain_temp_elm     /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%rain_temp_elm,ierr)
    if(elm_pf_idata%rain_temp_pf      /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%rain_temp_pf,ierr)

    if(elm_pf_idata%gflux_subsurf_elm /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%gflux_subsurf_elm,ierr)
    if(elm_pf_idata%gflux_subsurf_pf  /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%gflux_subsurf_pf,ierr)

    if(elm_pf_idata%sat_elms           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sat_elms,ierr)
    if(elm_pf_idata%sat_pfp            /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sat_pfp,ierr)
    if(elm_pf_idata%mass_elms          /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%mass_elms,ierr)
    if(elm_pf_idata%mass_pfp           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%mass_pfp,ierr)

    if(elm_pf_idata%temp_elm          /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%temp_elm,ierr)
    if(elm_pf_idata%temp_pf           /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%temp_pf,ierr)

    if(elm_pf_idata%sat_ice_elm       /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sat_ice_elm,ierr)
    if(elm_pf_idata%sat_ice_pf        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%sat_ice_pf,ierr)

    if(elm_pf_idata%h2osfc_elm        /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%h2osfc_elm,ierr)
    if(elm_pf_idata%h2osfc_pf         /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%h2osfc_pf,ierr)

    if(elm_pf_idata%area_top_face_pfp   /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_top_face_pfp,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%area_top_face_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_top_face_elms,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%area_proj_top_face_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%area_proj_top_face_elms,ierr);CHKERRQ(ierr)

    if(elm_pf_idata%pfgrid_nG2A_pfs  /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%pfgrid_nG2A_pfs,ierr)

    if(elm_pf_idata%eff_therm_cond_elm  /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%eff_therm_cond_elm,ierr)
    if(elm_pf_idata%eff_therm_cond_pf  /= PETSC_NULL_VEC) call VecDestroy(elm_pf_idata%eff_therm_cond_pf,ierr)

#ifdef PRINT_INTERNALFLOW
    if(elm_pf_idata%mflx_infl_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%mflx_infl_elms ,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%mflx_et_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%mflx_et_elms ,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%mflx_dew_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%mflx_dew_elms ,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%mflx_sub_snow_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%mflx_sub_snow_elms ,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%mflx_snowlyr_disp_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%mflx_snowlyr_disp_elms ,ierr);CHKERRQ(ierr)
    if(elm_pf_idata%mflx_drain_elms  /= PETSC_NULL_VEC) &
      call VecDestroy(elm_pf_idata%mflx_drain_elms ,ierr);CHKERRQ(ierr)
#endif

  end subroutine ELMPFLOTRANIDataDestroy

end module elm_pflotran_interface_data
