module dist_fn
!  use init_g, only: ginit
  use redistribute, only: redist_type
  implicit none
  public :: init_dist_fn, init_dist_fn_minimal
  public :: timeadv, exb_shear
  public :: getfieldeq, getmoms
  public :: getmoms_r
  public :: getmoms_r_temp
  public :: flux
!  public :: ginit, get_epar, get_heat
  public :: get_epar, get_heat
  public :: t0, omega0, gamma0, source0
  public :: reset_init, write_f, reset_physics
  public :: write_vp !Velocity-space slices for each Fourier mode (GGH 23JAN08)
  public :: par_spectrum, get_init_field, write_vspec
  public :: Hankel_transform
  public :: get_init_field_h
  public :: get_dens_vel, get_jext !GGH
  public :: get_verr, pass_nwrite, par_transform !MAB

  ! RN>
  public :: init_eq
  public :: gamtot, gamtot1, gamtot2
  public :: mom_coeff, ncnt_mom_coeff
  public :: mom_coeff_npara, mom_coeff_nperp
  public :: mom_coeff_tpara, mom_coeff_tperp
  public :: mom_shift_para, mom_shift_perp
  ! <RN

  public :: g_adjust
  public :: is_g_gnew

  private

  ! knobs
  real, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  integer, dimension (:), allocatable :: bd_exp ! nspec
  real :: poisfac, cvdrift, gbdrift, alfy
  real :: t0, omega0, gamma0, thetas, source0
  real :: phi_ext
  real :: aky_star, akx_star
  real :: g_exb, g_pvg, g_exb_start_time
  integer :: g_exb_start_timestep
  logical :: exit_if_g_exb_fails

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_zero = 2, &
       source_option_test1 = 3, source_option_hm_force = 4

  ! TT>
  logical :: nonad_zero ! outgoing boundary condition for h
  integer :: boundary_option_switch
  integer, parameter :: boundary_option_zero = 1, &
       boundary_option_self_periodic = 2
  ! <TT

  logical :: mult_imp, test, def_parity, even, test_bes
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.
  
  ! must be initialized appropriately in init_g
  logical :: is_g_gnew = .true.

!! k_parallel filter items
!  real, dimension(:), allocatable :: work, tablekp
!  real :: scale
!  integer :: nwork, ntablekp

  ! internal arrays

  real, dimension (:), allocatable :: wdrift
  ! (-g-layout-)

  real, dimension (:,:,:), allocatable :: wstar
  ! (naky,negrid,nspec) replicated

  ! fieldeq
  real, dimension (:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  ! (nakx,naky) replicated

  complex, dimension (:), allocatable :: a, b, r, ainv
  ! (-g-layout-)

  complex, dimension (:,:,:), allocatable :: g0, g_h, gnl_1, gnl_2, gnl_3
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  complex, dimension (:,:), allocatable :: gk0
  
  ! momentum conservation
  real, dimension (:,:,:), allocatable :: sq

  logical :: initialized = .false.
  logical :: initializing = .true.

  real, allocatable :: mom_coeff(:,:,:,:)
  real, allocatable :: mom_coeff_npara(:,:,:), mom_coeff_nperp(:,:,:)
  real, allocatable :: mom_coeff_tpara(:,:,:), mom_coeff_tperp(:,:,:)
  real, allocatable :: mom_shift_para(:,:,:), mom_shift_perp(:,:,:)
  integer, parameter :: ncnt_mom_coeff=8

  logical :: debug = .false.

  integer :: nwrite = 1, navg = 0  ! MAB

contains

  subroutine init_dist_fn_minimal
    ! <doc>
    !  minimal set of initialization to use distribution function arrays
    ! </doc>
    use mp, only: proc0
    use kgrids, only: init_kgrids, naky, nakx, akx, aky
    use theta_grid, only: init_theta_grid, ntgrid
    use le_grids, only: init_le_grids, nlambda, negrid
    use species, only: init_species, nspec
    use agk_layouts, only: init_dist_fn_layouts, init_agk_layouts
    use run_parameters, only: init_run_parameters
    implicit none

    call init_agk_layouts
    call init_species
    call init_theta_grid
    call init_kgrids
    call init_le_grids (accelerated_x, accelerated_v)
    call read_parameters

    if (test) then
       if (proc0) then
          write (*,*) 'nspecies = ', nspec
          write (*,*) 'nlambda = ', nlambda
          write (*,*) 'negrid = ', negrid
          write (*,*) 'nakx = ', nakx
          write (*,*) 'naky = ', naky
       end if
!!$       call finish_mp
!!$       stop
    end if

    call init_run_parameters
    call init_dist_fn_layouts (ntgrid, naky, nakx, nlambda, negrid, nspec)
    if (.not. test_bes) then
       call allocate_arrays
    end if

    call init_eq
  end subroutine init_dist_fn_minimal

  subroutine init_dist_fn
    use mp, only: finish_mp
    use collisions, only: init_collisions
    use hyper, only: init_hyper
    use nonlinear_terms, only: init_nonlinear_terms
    implicit none

    if (initialized) return
    initialized = .true.

    if (.not. test_bes) then
       call init_nonlinear_terms 
    end if

    call init_dist_fn_minimal

    call init_vpar

    call init_bessel
    call init_wdrift
    call init_wstar

    if (test_bes) then
       call bes_out
       call finish_mp
       stop
    end if

    call init_par_filter
    call init_collisions 
    call init_invert_rhs
    call init_fieldeq
    call init_hyper

    call init_mom_coeff
  end subroutine init_dist_fn

  subroutine bes_out

    use agk_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx, idx_local, proc_id
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use species, only: spec
    use le_grids, only: e
    use dist_fn_arrays, only: aj0, aj1vp2
    use mp, only: proc0, barrier, send, receive

    real :: b0, b1
    integer :: it, ie, il, is, ik, ig, iglo, unit

    if (proc0) then
       call get_unused_unit (unit)
       call open_output_file (unit, ".bes")
    end if


    do iglo=g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo, iglo) 
       it = it_idx(g_lo, iglo) 
       is = is_idx(g_lo, iglo) 
       ie = ie_idx(g_lo, iglo) 
       il = il_idx(g_lo, iglo)

       if (idx_local (g_lo, ik, it, il, ie, is)) then
          if (proc0) then 
             b0 = spec(is)%z*aj0(iglo)/spec(is)%temp
             b1 = aj1vp2(iglo)
          else
             call send (spec(is)%z*aj0(iglo)/spec(is)%temp, 0)
             call send (aj1vp2(iglo), 0)
          end if
       else if (proc0) then
          call receive (b0, proc_id(g_lo, iglo))
          call receive (b1, proc_id(g_lo, iglo))
       end if
       
       if (proc0) write (unit, "(i8,1x,3(1x,e13.6))") iglo, b0, b1, e(ie, is)
       
       call barrier
    end do

    if (proc0) call close_output_file (unit)

  end subroutine bes_out


  subroutine read_parameters
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use theta_grid, only: nperiod
    use text_options, only: text_option, get_option_value
    use species, only: nspec
    use mp, only: proc0, broadcast, mp_abort
    use agk_mem, only: alloc8
    use agk_mem, only: alloc8
    implicit none
    type (text_option), dimension (5), parameter :: sourceopts = &
         (/ text_option('default', source_option_full), &
            text_option('full', source_option_full), &
            text_option('zero', source_option_zero), &
            text_option('test1', source_option_test1), &
            text_option('hm', source_option_hm_force) /)
    character(20) :: source_option

    ! TT>
    type (text_option), dimension (4), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_self_periodic), &
            text_option('zero', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic) /)
    character(20) :: boundary_option
    ! <TT

    type (text_option), dimension (8), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
! eventually add in iphi00 = 0 option:
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg), &! /)
            text_option('zero', adiabatic_option_zero) /)
    character(30) :: adiabatic_option

    ! TT>
    ! namelist /dist_fn_knobs/ poisfac, adiabatic_option, mult_imp, test, def_parity, even, &
    namelist /dist_fn_knobs/ boundary_option, nonad_zero, &
         poisfac, adiabatic_option, mult_imp, test, def_parity, even, &
         ! <TT
         g_exb, test_bes, g_pvg, exit_if_g_exb_fails, &
         g_exb_start_time, g_exb_start_timestep
    
    namelist /source_knobs/ t0, omega0, gamma0, source0, &
           phi_ext, source_option, aky_star, akx_star, cvdrift, gbdrift, &
           alfy
    integer :: ierr, is, in_file
    logical :: exist
    real :: bd
    logical :: done = .false.
    integer :: ireaderr1=0, ireaderr2=0
    logical :: abort_readerr=.true.

    if (done) return
    done = .true.

    if (proc0) then
       ! TT>
       boundary_option = 'default'
       nonad_zero = .false.
       ! <TT
       adiabatic_option = 'default'
       poisfac = 0.0
       t0 = 100.0
       source0 = 1.0
       omega0 = 0.0
       gamma0 = 0.0
       aky_star = 0.0
       akx_star = 0.0
       phi_ext = 0.0
       g_exb = 0.0
       g_pvg=0.0
       g_exb_start_timestep = -1
       g_exb_start_time = -1
       exit_if_g_exb_fails = .false.
       cvdrift = 0.0
       gbdrift = 0.0
       alfy = 0.0
       mult_imp = .false.
       test = .false.
       test_bes = .false.
       def_parity = .false.
       even = .true.
       source_option = 'default'

       in_file = input_unit_exist("dist_fn_knobs", exist)
       if (exist) read (unit=in_file, nml=dist_fn_knobs, iostat=ireaderr1)
       in_file = input_unit_exist("source_knobs", exist)
       if (exist) read (unit=in_file, nml=source_knobs, iostat=ireaderr2)

       ierr = error_unit()
       ! TT>
       call get_option_value &
            (boundary_option, boundaryopts, boundary_option_switch, &
            ierr, "boundary_option in dist_fn_knobs")
       ! <TT
       call get_option_value &
            (source_option, sourceopts, source_option_switch, &
            ierr, "source_option in source_knobs")
       call get_option_value &
            (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in dist_fn_knobs")

    end if
    call broadcast (ireaderr1)
    call broadcast (ireaderr2)
    if (abort_readerr) then
       if (ireaderr1 > 0) then
          call mp_abort ('Read error at dist_fn_knobs')
       else if (ireaderr1 < 0) then
          call mp_abort ('End of file/record occurred at dist_fn_knobs')
       endif
       if (ireaderr2 > 0) then
          call mp_abort ('Read error at source_knobs')
       else if (ireaderr2 < 0) then
          call mp_abort ('End of file/record occurred at source_knobs')
       endif
    endif

    if (.not.allocated(fexp)) then
       allocate (fexp(nspec)); call alloc8(r1=fexp,v="fexp")
       allocate (bkdiff(nspec)); call alloc8(r1=bkdiff,v="bkdiff")
       allocate (bd_exp(nspec)); call alloc8(i1=bd_exp,v="bd_exp")
    end if

    if (proc0) call read_species_knobs

    ! TT>
    call broadcast (boundary_option_switch)
    call broadcast (nonad_zero)
    ! <TT
    call broadcast (adiabatic_option_switch)
    call broadcast (poisfac)
    call broadcast (t0)
    call broadcast (source0)
    call broadcast (omega0)
    call broadcast (gamma0)
    call broadcast (aky_star)
    call broadcast (akx_star)
    call broadcast (phi_ext)
    call broadcast (g_exb)
    call broadcast (g_pvg)
    call broadcast (exit_if_g_exb_fails)
    call broadcast (g_exb_start_timestep)
    call broadcast (g_exb_start_time)
    call broadcast (cvdrift)
    call broadcast (gbdrift)
    call broadcast (alfy)
    call broadcast (source_option_switch)
    call broadcast (fexp)
    call broadcast (bkdiff)
    call broadcast (bd_exp)
    call broadcast (mult_imp)
    call broadcast (test)
    call broadcast (test_bes)
    call broadcast (def_parity)
    call broadcast (even)

    if (mult_imp) then
       ! nothing -- fine for linear runs, but not implemented nonlinearly
    else
! consistency check for bkdiff
       bd = bkdiff(1)
       do is = 1, nspec
          if (bkdiff(is) /= bd) then
             if (proc0) write(*,*) 'Forcing bkdiff for species ',is,' equal to ',bd
             if (proc0) write(*,*) 'If this is a linear run, and you want unequal bkdiff'
             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
             bkdiff(is) = bd
          endif
       end do
    end if

  end subroutine read_parameters 

  subroutine read_species_knobs
    use species, only: nspec
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer :: is, unit
    do is = 1, nspec
       fexp(is) = 0.4
       bkdiff(is) = 0.0
       bd_exp(is) = 0
       
       call get_indexed_namelist_unit (unit, "dist_fn_species_knobs", is)
       call fill_species_knobs (unit, fexp(is), bkdiff(is), bd_exp(is))
       close (unit=unit)
    end do
  end subroutine read_species_knobs

  subroutine fill_species_knobs (unit, fexp_out, bakdif_out, bd_exp_out)
    implicit none
    integer, intent (in) :: unit
    real, intent (in out) :: fexp_out
    real, intent (in out) :: bakdif_out
    integer, intent (in out) :: bd_exp_out
    integer :: bd_exp
    real :: fexp, bakdif
    namelist /dist_fn_species_knobs/ fexp, bakdif, bd_exp

    fexp = fexp_out
    bakdif = bakdif_out
    bd_exp = bd_exp_out
    read (unit=unit, nml=dist_fn_species_knobs)
    fexp_out = fexp
    bd_exp_out = bd_exp
    bakdif_out = bakdif
  end subroutine fill_species_knobs


  subroutine init_vpar
    use dist_fn_arrays, only: vpa, vpar, vpac, vperp2
    use species, only: spec
    use theta_grid, only: ntgrid, delthet, gradpar
    use le_grids, only: e, al
    use agk_time, only: dtime
    use agk_layouts, only: g_lo, il_idx, ie_idx, is_idx
    use agk_mem, only: alloc8
    implicit none
    integer :: iglo, is
    real :: al1, e1

    if (.not.allocated(vpa)) then
       allocate (vpa  (2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r2=vpa , v="vpa")
       allocate (vpac (2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r2=vpac, v="vpac")
       allocate (vpar (2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r2=vpar, v="vpar")
       allocate (vperp2 ( g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r1=vperp2, v="vperp2")
    endif
    vpa = 0. ; vpac = 0. ; vperp2 = 0. ; vpar = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       al1 = al(il_idx(g_lo,iglo))
       e1 = e(ie_idx(g_lo,iglo),is)

       vpa(1,iglo) = sqrt(e1*max(0.0, 1.0 - al1))
       vpa(2,iglo) = - vpa(1,iglo)
       vperp2(iglo) = al1*e1

       ! <doc> INFO (RN):
       !  Following coding may be unnecessary because al1 (lambda) never
       !  be unity. ("1-al1" gets close to the machine epsilon x 100 
       !  if Nlambda is greater than ~1/sqrt(eps*100).
       !  Thus, there's no distinction between vpa and vpac, and we can
       !  remove vpac.
       ! </doc>
       if (1.0 - al1 < 100.0*epsilon(0.0)) then
          vpa(1,iglo) = 0.0
          vpa(2,iglo) = 0.0
       end if

       if  (1.0 - al1 < 0.0) then
          vpac(1,iglo) = 1.0
          vpac(2,iglo) = -1.0
       else
          vpac(1,iglo) = vpa(1,iglo)
          vpac(2,iglo) = vpa(2,iglo)
       end if

!       vpac(:,iglo) = 0.0

! TT: variable step AB3
!       vpar(1,iglo) = spec(is)%zstm*dtime * abs(gradpar)/delthet * vpac(1,iglo)
!       vpar(2,iglo) = spec(is)%zstm*dtime * abs(gradpar)/delthet * vpac(2,iglo)
       vpar(1,iglo) = spec(is)%zstm*dtime(1) * abs(gradpar)/delthet * vpac(1,iglo)
       vpar(2,iglo) = spec(is)%zstm*dtime(1) * abs(gradpar)/delthet * vpac(2,iglo)
! <TT
! ISSUE?  vpar(ntgrid) was set to zero.  Now independent of theta.  Problem anywhere?

    end do

  end subroutine init_vpar

  subroutine init_wdrift
    use species, only: nspec
    use theta_grid, only: ntgrid
    use kgrids, only: naky, nakx, aky
    use le_grids, only: negrid, nlambda, al, e
    use agk_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use agk_mem, only: alloc8
    use species, only: spec
    use dist_fn_arrays, only: aj0
    use agk_time, only: dtime
    implicit none

    integer :: iglo
    logical :: alloc = .true.
    integer :: il, ie, ik, is
    if (alloc) then
       allocate (wdrift(g_lo%llim_proc:g_lo%ulim_alloc)); call alloc8(r1=wdrift,v="wdrift")
    end if
    wdrift = 0.

!!$    do iglo = g_lo%llim_proc, g_lo%ulim_proc
!!$       wdrift(iglo) &
!!$            = wdrift_func(il_idx(g_lo,iglo), ie_idx(g_lo,iglo), &
!!$            ik_idx(g_lo,iglo), is_idx(g_lo,iglo))
!!$    end do
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il=il_idx(g_lo,iglo)
       ie=ie_idx(g_lo,iglo)
       ik=ik_idx(g_lo,iglo)
       is=is_idx(g_lo,iglo)
       wdrift(iglo) = - aky(ik) * dtime(1) * ( &
            cvdrift * e(ie,is) * (1.0 - al(il)) & ! curvature drift
            + 0.5 * gbdrift * e(ie,is) * al(il) & ! grad-B drift
            - aj0(iglo) * alfy * spec(is)%zstm * & ! Alfven wave
            sqrt(e(ie,is) * (1.0 - al(il))) )
    end do

    alloc = .false.

  contains

    function wdrift_func (il, ie, ik, is)
      use kgrids, only: aky
      use le_grids, only: e, al
      use agk_time, only: dtime
      implicit none
      real :: wdrift_func
      integer, intent (in) :: ik, il, ie, is

      wdrift_func = - aky(ik) * ( cvdrift * e(ie,is) * (1.0 - al(il)) &
! TT: variable step AB3
!           + 0.5 * gbdrift * e(ie,is) * al(il) ) * dtime
           + 0.5 * gbdrift * e(ie,is) * al(il) ) * dtime(1)
! <TT

    end function wdrift_func

  end subroutine init_wdrift

  subroutine init_wstar
    use species, only: nspec, spec
    use kgrids, only: naky, aky
    use le_grids, only: negrid, e
    use agk_time, only: dtime
    use agk_mem, only: alloc8
    implicit none
    integer :: ik, ie, is

    if(.not.allocated(wstar)) then
       allocate (wstar(naky,negrid,nspec)); call alloc8(r3=wstar,v="wstar")
    end if

    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             wstar(ik,ie,is) = -dtime(1)*aky(ik)/2. &
                  *(spec(is)%fprim+spec(is)%tprim*(e(ie,is)-1.5))
          end do
       end do
    end do
  end subroutine init_wstar

  subroutine init_bessel
    use dist_fn_arrays, only: aj0, aj1vp2, vperp2, aj1
    use species, only: spec
    use theta_grid, only: ntgrid
    use kgrids, only: kperp2, single
    use le_grids, only: e, al
    use agk_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use agk_mem, only: alloc8
    use spfunc, only: j0, j1
    implicit none
    integer :: ig, ik, it, il, ie, is
    integer :: iglo
    real :: arg

    logical :: done = .false.

    if (.not. done) then
      allocate (aj0 (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r1=aj0, v="aj0")
      allocate (aj1 (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r1=aj1, v="aj1")
      allocate (aj1vp2 (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (r1=aj1vp2, v="aj1vp2")
    end if 

    ! If single and g_exb is non zero, need to recalc the bessels every time
    ! init_bessel is called EGH

    if (done .and. .not. (single .and. g_exb .ne. 0.0) ) return
    done = .true.
    !write (*,*) "kperp2 in  bessel function is ", kperp2(1,1)

    aj0 = 0. ; aj1vp2 = 0. ; aj1 = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       arg = spec(is)%smz*sqrt(e(ie,is)*al(il)*kperp2(it,ik))  ! |B| = 1 assumed
       aj0(iglo) = j0(arg)
       aj1(iglo) = j1(arg)*arg
       aj1vp2(iglo) = j1(arg)*vperp2(iglo)*2.0
    end do

  end subroutine init_bessel

  subroutine init_par_filter
    use theta_grid, only: ntgrid, nperiod
    use agk_transforms, only: init_zf

    call init_zf (ntgrid, nperiod)

  end subroutine init_par_filter

  subroutine par_spectrum(an, an2)

    use agk_transforms, only: kz_spectrum
    use theta_grid, only: ntgrid
    use kgrids, only: naky, nakx

    complex, dimension(:,:,:) :: an, an2    
    integer :: it, ik
    real :: scale

    call kz_spectrum (an, an2, nakx, naky)
    scale = 1./real(4*ntgrid**2)  
    an2 = an2*scale

  end subroutine par_spectrum

  subroutine par_transform (an_in, an_out)

    use agk_transforms, only: kz_transform
    use theta_grid, only: ntgrid
    use agk_layouts, only: g_lo

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: an_in
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (out) :: an_out
    real :: scale

    call kz_transform (an_in, an_out)
    scale = 1./real(2*ntgrid)
    an_out = an_out*scale

  end subroutine par_transform

  subroutine init_invert_rhs
    use mp, only: proc0
    use dist_fn_arrays, only: vpa, vpar, vpac
    use species, only: spec
    use theta_grid, only: ntgrid, theta
    use kgrids, only: naky, nakx, aky
    use le_grids, only: negrid, nlambda
    use constants
    use agk_layouts, only: g_lo, il_idx, is_idx
    use agk_mem, only: alloc8
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is
    real :: vp, bd, wd

    if (.not.allocated(a)) then
       allocate (a (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c1=a, v="a")
       allocate (b (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c1=b, v="b")
       allocate (r (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c1=r, v="r")
       allocate (ainv (g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c1=ainv, v="ainv")
    endif
    a = 0. ; b = 0. ; r = 0. ; ainv = 0.
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       vp = vpar(1,iglo)
       bd = bkdiff(is)
       wd = wdrift(iglo)

! NOTE: vp was zero for theta=ntgrid in gs2 coding
       ainv(iglo) &
            = 1.0/(1.0 + bd &
            + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0+bd) + 2.0*vp))
       r(iglo) &
            = (1.0 - bd &
            + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0-bd) - 2.0*vp)) &
            *ainv(iglo)
       a(iglo) &
            = 1.0 + bd &
            + fexp(is)*spec(is)%tz*(-zi*wd*(1.0+bd) - 2.0*vp)
       b(iglo) &
            = 1.0 - bd &
            + fexp(is)*spec(is)%tz*(-zi*wd*(1.0-bd) + 2.0*vp)
    end do

    initializing = .false.

  end subroutine init_invert_rhs

  subroutine allocate_arrays 
    use kgrids, only: naky, nakx
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: g, gnew !, g_eq
    use agk_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    use agk_mem, only: alloc8

    implicit none
    logical :: alloc = .true.

    if (alloc) then
       allocate (g   (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=g, v="g")
       allocate (gnew(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=gnew, v="gnew")
       allocate (g0  (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=g0, v="g0")
!       if (store_eq) then
!          allocate (g_eq(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=g_eq, v="g_eq")
!       end if

       allocate (gk0 (2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c2=gk0, v="gk0")

       if (nonlin) then
          allocate (gnl_1 (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=gnl_1, v="gnl_1")
          allocate (gnl_2 (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=gnl_2, v="gnl_2")
          allocate (gnl_3 (-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=gnl_3, v="gnl_3")
          gnl_1 = 0. ; gnl_2 = 0. ; gnl_3 = 0.
       else
          allocate (gnl_1(1,2,1)); call alloc8(c3=gnl_1,v="gnl_1")
          allocate (gnl_2(1,2,1)); call alloc8(c3=gnl_2,v="gnl_2")
          allocate (gnl_3(1,2,1)); call alloc8(c3=gnl_3,v="gnl_3")
       end if
    endif

    g = 0. ; gnew = 0. ; g0 = 0. ; gk0 = 0.

!    if(store_eq) g_eq=0.
    alloc = .false.
  end subroutine allocate_arrays

  subroutine timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, mode)

    use theta_grid, only: ntgrid
    use collisions, only: solfp1!, collisions_test => test
    use dist_fn_arrays, only: gnew, g
    use nonlinear_terms, only: add_nonlinear_terms
    use hyper, only: hyper_diff
    use run_parameters, only: fphi, fbpar
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, optional, intent (in) :: mode
    integer :: modep
    integer :: diagnostics = 1

    modep = 0
    if (present(mode)) modep = mode

    if (modep <= 0) then
       call add_nonlinear_terms (gnl_1, gnl_2, gnl_3, &
            phi, apar, bpar, istep, bkdiff(1), fexp(1))
! TT: replace invert_rhs for collision test mode
!       if (collisions_test) then
!          gnew = g
!      else
! <TT
       call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
! TT>
!       end if
! <TT
       call hyper_diff (gnew, phinew, bparnew)

! MAB>
! don't want to calculate collisional heating unless we are going
! to use it for some diagnostic
       !!! RN
       call g_adjust (gnew, phinew, bparnew, fphi, fbpar, is_g_gnew)
       call g_adjust (g, phi, bpar, fphi, fbpar)
       !!!
       if (nwrite-mod(istep+nwrite-1,nwrite) <= navg) then
          call solfp1 (gnew, g, g0, phi, apar, bpar, phinew, aparnew, bparnew, diagnostics)
       else
          call solfp1 (gnew, g, g0, phi, apar, bpar, phinew, aparnew, bparnew)
       end if
       call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar, is_g_gnew)
       call g_adjust (g, phi, bpar, -fphi, -fbpar)
! <MAB
       
       if (def_parity) then
          if (even) then
             gnew(-ntgrid:-1, 1,:) = gnew( ntgrid: 1:-1,2,:)
             gnew( 1: ntgrid, 1,:) = gnew(-1:-ntgrid:-1,2,:)
          else
             gnew( 1: ntgrid, 1,:) = -gnew(-1:-ntgrid:-1,2,:)
             gnew(-ntgrid:-1, 1,:) = -gnew( ntgrid: 1:-1,2,:)
          end if
       end if
       
    end if

  end subroutine timeadv

  subroutine exb_shear (g0, phi, apar, bpar, istep)

    use mp, only: proc0
    use agk_layouts, only: g_lo, is_kx_local, idx_local, idx
    use file_utils, only: error_unit
    use theta_grid, only: ntgrid
    use kgrids, only: akx, aky, naky, ikx, nakx, box, single, calculate_kgrids
    use le_grids, only: negrid, nlambda
    use species, only: nspec
    use run_parameters, only: use_Phi, use_Apar, use_Bpar
    use dist_fn_arrays, only: kx_shift
    use agk_time, only: dtime, time
    use agk_mem, only: alloc8
    use mp, only: mp_abort
    use agk_mem, only: alloc8

    integer, intent (in) :: istep
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    integer, dimension(:), allocatable, save :: jump, ikx_indexed
    integer, dimension(1), save :: itmin
    integer :: ik, it, ie, is, il, ig, isgn, to_iglo, from_iglo, ierr, ik_lower
    real :: dkx, gdt
    logical, save :: exb_first = .true.
    logical, save :: kx_local
    real, save :: wait_time = 0.

    if (abs(g_exb)<epsilon(0.0)) return

    
    if (single) then
      if (g_exb_start_timestep > 0) then 
        if (istep < g_exb_start_timestep) then 
          return
        else if (g_exb_start_timestep == istep) then
          wait_time = time
        end if 
      end if 
      if (g_exb_start_time >= 0) then 
        if (time < g_exb_start_time) then 
          return
        else if (wait_time .eq. 0.0) then
          wait_time = time
        end if 
      end if

      call calculate_kgrids(g_exb, time - wait_time)
      call init_bessel
      call init_fieldeq
      return
    end if 


! If not in box configuration, return
    if (.not.box) then
      if (exit_if_g_exb_fails) then
        call mp_abort("ExB shear only works in box configuration")
      else
        write (error_unit(), *) "ExB shear only works in box configuration"
        return
      end if
    end if


       

! If kx data is not local, no ExB shear will be included.

! Initialize kx_shift, jump, idx_indexed

    if (exb_first) then
       call is_kx_local (negrid, nspec, nlambda, naky, nakx, kx_local)
       if (.not. kx_local) then
          g_exb = 0.
          if (proc0) then
             ierr = error_unit()
             write (ierr, fmt="('Non-zero g_ExB not implemented for this layout.  g_ExB set to zero.')")
          end if
          if (exit_if_g_exb_fails) then 
            call mp_abort("kx_local is false: g_exb will not work. Aborting.")
          end if
          return
       end if
       exb_first = .false.
       allocate (kx_shift(naky)); call alloc8(r1=kx_shift,v="kx_shift")
       allocate (jump(naky)); call alloc8(i1=jump,v="jump")
       kx_shift = 0.
       jump = 0

       allocate (ikx_indexed(nakx)); call alloc8(i1=ikx_indexed,v="ikx_indexed")
       itmin = minloc (ikx)

       do it=itmin(1), nakx
          ikx_indexed (it+1-itmin(1)) = it
       end do

       do it=1,itmin(1)-1
          ikx_indexed (nakx - itmin(1) + 1 + it)= it
       end do
    end if

    dkx = akx(2)
    
! BD: To do: Put the right timestep in here.
!
! For now, approximate Greg's dt == 1/2 (t_(n+1) - t_(n-1))
! with dtime.  
!
! Note: at first time step, there is a difference of a factor of 2.
!

! TT: variable step AB3
!    gdt = dtime
    gdt = dtime(1)
! <TT

    if (g_exb_start_timestep > istep) return
    if (g_exb_start_time >= 0 .and. time < g_exb_start_time) return


! kx_shift is a function of time.   Update it here:


    do ik=1, naky
       kx_shift(ik) = kx_shift(ik) - aky(ik)*g_exb*gdt
       jump(ik) = nint(kx_shift(ik)/dkx)
       kx_shift(ik) = kx_shift(ik) - jump(ik)*dkx
    end do

! BD: To do: Should save kx_shift array in restart file

    ! EGH: Don't bother to do ik == 1 unless aky(1) is finite 
    if (aky(1) .ne. 0.0) then
      ik_lower = 1
    else
      ik_lower = 2
    end if 

    do ik = naky, ik_lower, -1
       if (jump(ik) == 0) exit

       if (jump(ik) < 0) then

          if (use_Phi) then
             do it = 1, nakx + jump(ik)
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = phi (ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = nakx + jump(ik) + 1, nakx
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (use_Apar) then
             do it = 1, nakx + jump(ik)
                do ig=-ntgrid,ntgrid
                   apar(ig,ikx_indexed(it),ik) = apar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = nakx + jump(ik) + 1, nakx
                do ig=-ntgrid,ntgrid
                   apar (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (use_Bpar) then
             do it = 1, nakx + jump(ik)
                do ig=-ntgrid,ntgrid
                   bpar(ig,ikx_indexed(it),ik) = bpar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = nakx + jump(ik) + 1, nakx
                do ig=-ntgrid,ntgrid
                   bpar (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For some data layouts, there is no communication required.  Try such a case first.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (kx_local) then
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda
                      do it = 1, nakx + jump(ik)
                         
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)
                         
                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = g0(ig,isgn,from_iglo)
                               end do
                            end do
                         end if
                      end do

                      do it = nakx + jump(ik) + 1, nakx
                         
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         
                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = 0.
                               end do
                            end do
                         end if
                      end do
                   end do
                end do
             end do
          end if

       else  ! case for jump(ik) > 0

          if (use_Phi) then
             do it = nakx, 1+jump(ik), -1
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = phi (ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = jump(ik), 1, -1
                do ig=-ntgrid,ntgrid
                   phi (ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (use_Apar) then
             do it = nakx, 1+jump(ik), -1
                do ig=-ntgrid,ntgrid
                   apar(ig,ikx_indexed(it),ik) = apar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = jump(ik), 1, -1
                do ig=-ntgrid,ntgrid
                   apar(ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

          if (use_Bpar) then
             do it = nakx, 1+jump(ik), -1
                do ig=-ntgrid,ntgrid
                   bpar(ig,ikx_indexed(it),ik) = bpar(ig,ikx_indexed(it-jump(ik)),ik)
                end do
             end do
             do it = jump(ik), 1, -1
                do ig=-ntgrid,ntgrid
                   bpar(ig,ikx_indexed(it),ik) = 0.
                end do
             end do
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For some data layouts, there is no communication required.  Try such a case first.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (kx_local) then
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda
                      do it = nakx, 1+jump(ik), -1

                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)

                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = g0(ig,isgn,from_iglo)
                               end do
                            end do
                         end if
                      end do

                      do it = jump(ik), 1, -1

                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)

                         if (idx_local (g_lo, to_iglo)) then
                            do isgn=1,2
                               do ig=-ntgrid,ntgrid
                                  g0(ig,isgn,to_iglo) = 0.
                               end do
                            end do
                         end if
                      end do
                   end do
                end do
             end do
          end if

       end if

    end do
       
  end subroutine exb_shear

  subroutine g_adjust (g, phi, bpar, facphi, facbpar, is_g)
    ! g -> h: facphi= 1, facbpar= 1
    ! h -> g: facphi=-1, facbpar=-1
    use species, only: spec
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0, aj1vp2
    use agk_layouts, only: g_lo, ik_idx, it_idx, is_idx
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (in) :: facphi, facbpar
    logical, intent(inout), optional :: is_g
    integer :: iglo, ig, ik, it, is
    complex :: adj

    ! invert is_g [if is_g, g is g, otherwise, g is h]
    if (present(is_g)) is_g = .not. is_g

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do ig = -ntgrid, ntgrid
          adj = aj1vp2(iglo)*bpar(ig,it,ik)*facbpar &
               + spec(is)%z*phi(ig,it,ik)*aj0(iglo) &
               /spec(is)%temp*facphi
          g(ig,1,iglo) = g(ig,1,iglo) + adj
          g(ig,2,iglo) = g(ig,2,iglo) + adj
       end do
    end do
  end subroutine g_adjust

  subroutine get_source_term &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        isgn, iglo, sourcefac, source)
    use dist_fn_arrays, only: aj0, aj1vp2, vpar, vpac, g
    use theta_grid, only: ntgrid, theta
    use kgrids, only: aky, akx
    use le_grids, only: nlambda, e, negrid
    use species, only: spec, nspec
    use run_parameters, only: fphi, fapar, fbpar, aant
    use agk_time, only: dtime
    use agk_layouts, only: g_lo, ik_idx, it_idx, ie_idx, is_idx
    use nonlinear_terms, only: nonlin
    use hyper, only: D_res
    use antenna, only: drive_upar, drive_upar_mode
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: isgn, iglo
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:), intent (out) :: source
    real :: timep
    complex, dimension (-ntgrid:ntgrid) :: upar_ext

    integer :: ig, ik, it, ie, is
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg!, bpargavg !GGH added bpargavg

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    phigavg  = (fexp(is)*phi(:,it,ik) + (1.0-fexp(is))*phinew(:,it,ik)) &
                *aj0(iglo)*fphi &
             + (fexp(is)*bpar(:,it,ik) + (1.0-fexp(is))*bparnew(:,it,ik))&
                *aj1vp2(iglo)*fbpar*spec(is)%tz
    apargavg = (fexp(is)*apar(:,it,ik) + (1.0-fexp(is))*aparnew(:,it,ik)) &
                *aj0(iglo)*fapar

! source term in finite difference equations
    select case (source_option_switch)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Default choice: solve self-consistent equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_full)
       call set_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms + include external i omega_d Phi * F_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_hm_force)
       call set_source
       if (istep > 0 .and. aky(int(aky_star)) == aky(ik) .and. akx(int(akx_star)) == akx(it)) &
            source(:ntgrid-1) = source(:ntgrid-1) - zi*phi_ext*sourcefac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Include no source term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_zero)
       source = 0.0

    end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Add upar antenna
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (aant) then
       if (any(drive_upar_mode)) then
          call drive_upar(upar_ext,is,it,ik)
          source = source + upar_ext*vpar(isgn,iglo)
       end if
    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do matrix multiplications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (isgn == 1) then
       do ig = -ntgrid, ntgrid-1
          source(ig) = source(ig) + b(iglo)*g(ig,1,iglo) + a(iglo)*g(ig+1,1,iglo)
       end do
    else
       do ig = -ntgrid, ntgrid-1
          source(ig) = source(ig) + a(iglo)*g(ig,2,iglo) + b(iglo)*g(ig+1,2,iglo)
       end do
    end if

    source(ntgrid) = source(-ntgrid)

  contains

    subroutine set_source

      use agk_mem, only: alloc8
      complex :: apar_p, apar_m, phi_p, phi_m!, bpar_p !GGH added bpar_p
      real, dimension(:), allocatable, save :: ufac
      real :: bd, bdfac_p, bdfac_m
      integer :: i_s
      logical :: first = .true.
! TT: variable step AB3
      real :: deltmp, c1, c2, c3
! <TT

      if (first) then
         first = .false.
         allocate (ufac(nspec)); call alloc8(r1=ufac,v="ufac")
         ufac = 2.0*spec%uprim
      endif

! try fixing bkdiff dependence
      bd = bkdiff(1)

      bdfac_p = 1.+bd*(3.-2.*real(isgn))
      bdfac_m = 1.-bd*(3.-2.*real(isgn))

      do ig = -ntgrid, ntgrid-1
         phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
         phi_m = phigavg(ig+1)-phigavg(ig)
         apar_p = bdfac_p*apargavg(ig+1)+bdfac_m*apargavg(ig)
         apar_m = aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
              -apar(ig+1,it,ik)-apar(ig,it,ik)
 
         source(ig) = -2.0 * vpar(isgn,iglo) * phi_m &
              - spec(is)%zstm * vpac(isgn,iglo) &
              * (aj0(iglo)*apar_m + D_res(it,ik)*apar_p) &
              -zi*wdrift(iglo)*phi_p & ! J0 factor included in phi_p
              ! TT: Do we need bpar term?
              + zi*(wstar(ik,ie,is) &
! TT: variable step AB3
!              + vpac(isgn,iglo)*dtime*aky(ik)/2.0*ufac(is)) &
              + vpac(isgn,iglo)*dtime(1)*aky(ik)/2.0*ufac(is) &
              -2.0*vpac(isgn,iglo)*dtime(1)*aky(ik)/2.0*g_pvg) &
! <TT
              *(phi_p - apar_p*spec(is)%stm*vpac(isgn,iglo)) 
      end do
        
      if (nonlin) then         
         select case (istep)
         case (0)
            ! nothing
         case (1) ! Euler
            do ig = -ntgrid, ntgrid-1
! TT: variable step AB3
!               source(ig) = source(ig) - 0.5*dtime*gnl_1(ig,isgn,iglo)
               source(ig) = source(ig) - 0.5*dtime(1)*gnl_1(ig,isgn,iglo)
! <TT
            end do
         case (2) ! AB2
!!$ TT: variable step AB2
!!$            do ig = -ntgrid, ntgrid-1
!!$               source(ig) = source(ig) - 0.5*dtime*( &
!!$                    1.5*gnl_1(ig,isgn,iglo) - 0.5*gnl_2(ig,isgn,iglo))
!!$            end do
            c1 = dtime(1) * (dtime(1) + 2.0*dtime(2)) / (2.0*dtime(2))
            c2 = - dtime(1)**2 / (2.0*dtime(2))
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) - 0.5*( &
                    c1*gnl_1(ig,isgn,iglo) + c2*gnl_2(ig,isgn,iglo) )
            end do
!!$ <TT
         case default ! AB3

!!$ TT: variable step AB3
!!$            do ig = -ntgrid, ntgrid-1
!!$               source(ig) = source(ig) - 0.5*dtime*( &
!!$                    (23./12.) * gnl_1(ig,isgn,iglo) &
!!$                    - (4./3.) * gnl_2(ig,isgn,iglo) &
!!$                    + (5./12.) * gnl_3(ig,isgn,iglo))
!!$            end do
            deltmp = dtime(1) + dtime(2)
            c1 = ( deltmp**2 * (deltmp/3. + dtime(3)/2.) / dtime(2) &
                 & - dtime(2)*(dtime(2)/3.+dtime(3)/2.) ) / (dtime(2)+dtime(3))
            deltmp = dtime(2) + dtime(3)
            c2 = - dtime(1)**2 * (dtime(1)/3. + deltmp/2.) / dtime(2) / dtime(3)
            c3 = dtime(1)**2 * (dtime(1)/3. + dtime(2)/2.) / dtime(3) / deltmp
            ! The above calculation may be avoided to execute for each iglo
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) - 0.5*( &
                    c1 * gnl_1(ig,isgn,iglo) &
                    + c2 * gnl_2(ig,isgn,iglo) &
                    + c3 * gnl_3(ig,isgn,iglo) )
            end do
!!$ <TT
         end select
      end if

    end subroutine set_source

  end subroutine get_source_term

  subroutine invert_rhs_1 &
       (phi, apar, bpar, phinew, aparnew, bparnew, istep, &
        iglo, sourcefac)
    ! TT>
    !use dist_fn_arrays, only: gnew
    use dist_fn_arrays, only: gnew, aj0, aj1vp2
    use run_parameters, only: fphi, fbpar
    ! <TT
    use run_parameters, only: ieqzip
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda
    use kgrids, only: aky, nakx
    ! TT>
    !use agk_layouts, only: g_lo, ik_idx, it_idx
    use agk_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use species, only: spec
    ! <TT
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, intent (in) :: iglo
    complex, intent (in) :: sourcefac

    ! TT>
    !integer :: ig, isgn
    integer :: ig, isgn, it, ik, is
    complex :: adjleft, adjright
    ! <TT
    complex :: ftmp
    complex, dimension (-ntgrid:ntgrid,2) :: source, g1

    ! TT>
    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    !if(ieqzip(it_idx(g_lo,iglo),ik_idx(g_lo,iglo))==0) return
    if (ieqzip(it,ik)==0) return
    ! <TT

    do isgn = 1, 2
       call get_source_term (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, isgn, iglo, sourcefac, source(:,isgn))
    end do

    gnew(:,:,iglo) = 0.0              ! gnew is the inhomogeneous solution

    ! TT>
    if (boundary_option_switch == boundary_option_self_periodic) then
       g1 = 0.0                          ! g1 is the homogeneous solution 

       g1(-ntgrid,1) = 1.0               ! initial conditions for g1
       g1( ntgrid,2) = 1.0               ! initial conditions for g1
    else if (nonad_zero) then
       ! TT (Jan 4, 2013): ported CMR's implementation of outgoing boundary
       ! condition for h from GS2
       is = is_idx(g_lo,iglo)
       adjleft = aj1vp2(iglo) * bparnew(-ntgrid,it,ik) * fbpar &
            + spec(is)%z * phinew(-ntgrid,it,ik) * aj0(iglo) &
            /spec(is)%temp * fphi
       gnew(-ntgrid,1,iglo) = gnew(-ntgrid,1,iglo) - adjleft
       adjright = aj1vp2(iglo) * bparnew(ntgrid,it,ik) * fbpar &
            + spec(is)%z * phinew(ntgrid,it,ik) * aj0(iglo) &
            /spec(is)%temp * fphi
       gnew(ntgrid,2,iglo) = gnew(ntgrid,2,iglo) - adjright
    end if
    ! <TT

    !RN> forward elimination
                                      ! time advance vpar < 0 inhomogeneous part:
    do ig = ntgrid-1, -ntgrid, -1
       gnew(ig,2,iglo) = -gnew(ig+1,2,iglo)*r(iglo) + ainv(iglo)*source(ig,2)
    end do

                                      ! time advance vpar > 0 inhomogeneous part:
    do ig = -ntgrid, ntgrid-1
       gnew(ig+1,1,iglo) = -gnew(ig,1,iglo)*r(iglo) + ainv(iglo)*source(ig,1)
    end do

    !RN<

    !RN> taking care of boundary component
    !    (equivalent to backward substitution)

    ! TT: put homogeneous part and periodic correction into if statement
    if (boundary_option_switch == boundary_option_self_periodic) then
       do ig = ntgrid-1, -ntgrid, -1
          g1(ig,2) = -g1(ig+1,2)*r(iglo) ! time advance vpar < 0 homogeneous part:
       end do

       do ig = -ntgrid, ntgrid-1
          g1(ig+1,1) = -g1(ig,1)*r(iglo) ! time advance vpar > 0 homogeneous part:
       end do
    
       call self_periodic                ! add correct amount of homogeneous solution now
    end if
    ! <TT

    !RN<

    if (def_parity) then
       if (even) then
          gnew(-ntgrid:-1,1,iglo) = gnew( ntgrid:1:-1,2,iglo)
          gnew(1:ntgrid, 1,iglo) = gnew(-1:-ntgrid:-1,2,iglo)
       else
          gnew(1:ntgrid, 1,iglo) = -gnew(-1:-ntgrid:-1,2,iglo)
          gnew(-ntgrid:-1,1,iglo) = -gnew( ntgrid:1:-1,2,iglo)
       end if
    end if

  contains

    subroutine self_periodic

      ! g1(ntgrid,1)=1 only if Vpar=0
      if (g1(ntgrid,1) /= 1.) then
         ! gnew(-ntgrid,...)=0
         ftmp = (gnew(ntgrid,1,iglo) - gnew(-ntgrid,1,iglo))/(1.0 - g1(ntgrid,1))
         gnew(:,1,iglo) = gnew(:,1,iglo) + ftmp*g1(:,1)
      end if

      if (g1(-ntgrid,2) /= 1.) then
         ftmp = (gnew(-ntgrid,2,iglo) - gnew(ntgrid,2,iglo))/(1.0 - g1(-ntgrid,2))
         gnew(:,2,iglo) = gnew(:,2,iglo) + ftmp*g1(:,2)
      end if

    end subroutine self_periodic

  end subroutine invert_rhs_1

  subroutine invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    ! <doc>
    !  Gyrokinetic Solver: see Sec. 3.3 of AstroGK paper
    !                      solve (38) for g_{i+1}^{n+1} 
    ! </doc>
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use agk_layouts, only: g_lo
    use agk_time, only: time
    use constants
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep

    integer :: iglo

    real :: timep
    complex :: sourcefac

    if (time > t0) then
       sourcefac = source0*exp(-zi*omega0*time+gamma0*time)
    else
       sourcefac = (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time+gamma0*time)
    end if


! This loop is completely parallelizable over the iglo index.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       call invert_rhs_1 (phi, apar, bpar, phinew, aparnew, bparnew, &
            istep, iglo, sourcefac)
    end do

  end subroutine invert_rhs

  subroutine getan (antot, antota, antotp)
    use dist_fn_arrays, only: vpa, aj0, aj1vp2, gnew
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species
    use run_parameters, only: beta, use_Phi, use_Apar, use_Bpar
    use agk_layouts, only: g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    real, dimension (nspec) :: wgt

    integer :: isgn, iglo, ig

    if (use_Phi) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do

       wgt = spec%z*spec%dens
       call integrate_species (g0, wgt, antot)

    end if

    if (use_Apar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = aj0(iglo)*vpa(isgn,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       
       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g0, wgt, antota)

    end if

    if (use_Bpar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = 0.5*aj1vp2(iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do

       wgt = spec%temp*spec%dens
       call integrate_species (g0, wgt, antotp)

    end if
  end subroutine getan

  subroutine getan2 (antot, antota, antotp, adjust_in)
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use run_parameters, only: beta, use_Phi, use_Apar, use_Bpar
    use dist_fn_arrays, only: gnew
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew
    use mp, only: iproc
    use agk_mem, only: alloc8, dealloc8
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (out) :: antot, antota, antotp
    logical, optional :: adjust_in
    logical :: adjust_loc
    integer :: is

    integer :: nakx,naky
    complex, allocatable :: dens(:,:,:,:),upar(:,:,:,:),ppp(:,:,:,:)

    adjust_loc=.true.
    if (present(adjust_in)) adjust_loc=adjust_in

    nakx=size(antot,2)
    naky=size(antot,3)
    allocate(dens(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=dens,v="dens")
    allocate(upar(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=upar,v="upar")
    allocate(ppp(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=ppp,v="ppp")
    dens(:,:,:,:)=cmplx(0.,0.)
    upar(:,:,:,:)=cmplx(0.,0.)
    ppp(:,:,:,:)=cmplx(0.,0.)
    
    ! adjust=F/T moments of g/h
    call getmoms_r (dens=dens,uz=upar,ppp=ppp, adjust_in=adjust_loc)
    
    antot(:,:,:)=cmplx(0.,0.)
    antota(:,:,:)=cmplx(0.,0.)
    antotp(:,:,:)=cmplx(0.,0.)

    if (use_Phi) then
       do is=1,nspec
          antot(:,:,:)=antot(:,:,:)+dens(:,:,:,is)*spec(is)%z
       end do
    endif

    if (use_Apar) then
       do is=1,nspec
          antota(:,:,:)=antota(:,:,:) &
               & + 2.*beta*upar(:,:,:,is)*spec(is)%z*spec(is)%dens
       end do
    endif

    if (use_Bpar) then
       do is=1,nspec
          antotp(:,:,:)=antotp(:,:,:) + .5 * ppp(:,:,:,is)
       end do
    endif

    call dealloc8(c4=dens,v="dens")
    call dealloc8(c4=upar,v="upar")
    call dealloc8(c4=ppp,v="ppp")
  end subroutine getan2

  subroutine getmoms (ntot, density, upar, uperp, tpar, tperp)
    use dist_fn_arrays, only: vpa, vperp2, aj0, gnew
    use agk_layouts, only: g_lo, is_idx, ik_idx, it_idx
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    use run_parameters, only: use_Phi, use_Bpar
    use dist_fn_arrays, only: aj1vp2
    use fields_arrays, only: phinew, bparnew
    implicit none
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density, &
         upar, uperp, tpar, tperp, ntot

    integer :: ik, it, isgn, is, iglo, ig

! returns moment integrals to PE 0

! total density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1,2
          g0(:,isgn,iglo) = aj0(iglo) * gnew(:,isgn,iglo)
       end do
    end do

    if (use_Phi) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) = g0(:,isgn,iglo) + &
                  & (aj0(iglo)**2-1.0) * spec(is)%zt * phinew(:,it,ik)
          end do
       end do
    endif

    if (use_Bpar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          do isgn=1, 2
             g0(:,isgn,iglo) = g0(:,isgn,iglo) + &
                  & aj1vp2(iglo)*aj0(iglo) * bparnew(:,it,ik)
          end do
       end do
    end if
    call integrate_moment (g0, ntot)

! guiding center density
    call integrate_moment (gnew, density)

! guiding center uperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = sqrt(vperp2(iglo))*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, uperp)

! guiding center upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = vpa(isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, upar)

! guiding center tpar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*vpa(isgn,iglo)*g0(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g0, tpar)
    tpar = tpar - density

! guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g0(ig,isgn,iglo) = vperp2(iglo)*gnew(ig,isgn,iglo)
          end do
       end do
    end do

    call integrate_moment (g0, tperp)
    tperp = tperp - density

    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       density(:,:,:,is)=density(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       uperp(:,:,:,is)=uperp(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
    end do

  end subroutine getmoms

  subroutine getmoms_r (dens, ux,uy,uz, pxx,pyy,pzz,pxy,pyz,pzx, ppp, qzzz, qzpp,&
       & adjust_in)
    ! <doc>
    !  Take moments up to 2nd order at particle position (small r).
    !  This is a just wrapper of getmoms_r_* routines.
    !  This routine uses gnew, phinew, aparnew, bparnew.
    !  If adjust_loc=T, gnew is translated to h, then moments are taken.
    !  Otherwise, moments of gnew itself are taken.
    ! </doc>
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    use fields_arrays, only: phinew, bparnew
    use run_parameters, only: fphi, fbpar
    implicit none
    complex, intent (out), optional :: dens(-ntgrid:,:,:,:)
    complex, intent (out), optional :: &
         & ux(-ntgrid:,:,:,:), uy(-ntgrid:,:,:,:),  uz(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ppp(-ntgrid:,:,:,:), &
         & pxx(-ntgrid:,:,:,:), pyy(-ntgrid:,:,:,:), pzz(-ntgrid:,:,:,:), &
         & pxy(-ntgrid:,:,:,:), pyz(-ntgrid:,:,:,:), pzx(-ntgrid:,:,:,:)
    complex, intent (out), optional :: qzzz(-ntgrid:,:,:,:), qzpp(-ntgrid:,:,:,:)

    logical, intent (in), optional :: adjust_in
    logical :: adjust_loc = .true.

    adjust_loc=.true.
    if(present(adjust_in)) adjust_loc = adjust_in

    ! <doc>
    !  g -> h
    ! </doc>
    if (adjust_loc) &
         & call g_adjust(gnew, phinew, bparnew, fphi, fbpar, is_g_gnew)

    if (present(dens)) call getmoms_r_0 (dens)
    if (present(uz))   call getmoms_r_1para (uz)
    if (present(pzz))  call getmoms_r_2parapara (pzz)
    if (present(ux)) then
       if (present(uy)) then
          call getmoms_r_1perp (ux=ux,uy=uy)
       else
          call getmoms_r_1perp (ux=ux)
       end if
    else
       if (present(uy)) then
          call getmoms_r_1perp (uy=uy)
       end if
    end if
    if (present(ppp)) then
       call getmoms_r_2perptheta (ppp=ppp)
    end if
    if (present(pxx)) then
       if (present(pyy)) then
          if (present(pxy)) then
             call getmoms_r_2perptheta (pxx=pxx,pyy=pyy,pxy=pxy)
          else
             call getmoms_r_2perptheta (pxx=pxx,pyy=pyy)
          end if
       else
          if (present(pxy)) then
             call getmoms_r_2perptheta (pxx=pxx,pxy=pxy)
          else
             call getmoms_r_2perptheta (pxx=pxx)
          end if
       end if
    else
       if (present(pyy)) then
          if (present(pxy)) then
             call getmoms_r_2perptheta (pyy=pyy,pxy=pxy)
          else
             call getmoms_r_2perptheta (pyy=pyy)
          end if
       else
          if (present(pxy)) then
             call getmoms_r_2perptheta (pxy=pxy)
          end if
       end if
    end if
    if (present(pyz)) then
       if (present(pzx)) then
          call getmoms_r_2thetapara (pyz=pyz,pzx=pzx)
       else
          call getmoms_r_2thetapara (pyz=pyz)
       end if
    else
       if (present(pzx)) then
          call getmoms_r_2thetapara (pzx=pzx)
       end if
    end if
    if (present(qzzz)) call getmoms_r_3(qzzz=qzzz)
    if (present(qzpp)) call getmoms_r_3(qzpp=qzpp)
    ! <doc>
    !  h -> g
    ! </doc>
    if (adjust_loc) &
         & call g_adjust(gnew, phinew, bparnew, -fphi, -fbpar, is_g_gnew)

  end subroutine getmoms_r

  subroutine getmoms_r_0 (dens)
    ! <doc>
    !  zeroth order moment at particle position (small r)
    !    n = n0 * int ... J0 dv^3
    ! </doc>
    use dist_fn_arrays, only: aj0, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    complex, intent (out) :: dens(-ntgrid:,:,:,:)
    integer :: isgn, iglo, is
    g0 = 0. ! g0 is common in this module
    dens = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(iglo)*gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, dens, 1)
    do is=1,nspec
       dens(:,:,:,is)=dens(:,:,:,is)*spec(is)%dens
    end do
  end subroutine getmoms_r_0

  subroutine getmoms_r_1para (uz)
    ! <doc>
    !  first order parallel moment at particle position (small r)
    !    uz = vth0 * int ... J0 v_para dv^3
    ! </doc>
    use dist_fn_arrays, only: aj0, vpa, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    complex, intent (out) :: uz(-ntgrid:,:,:,:)
    integer :: isgn, iglo, is
    g0 = 0. ! g0 is common in this module
    uz = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj0(iglo)*vpa(isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, uz, 1)
    do is=1,nspec
       uz(:,:,:,is)=uz(:,:,:,is)*spec(is)%stm
    end do
  end subroutine getmoms_r_1para

  subroutine getmoms_r_1perp (ux,uy)
    ! <doc>
    !  first order perpendicular moment at particle position (small r)
    !  and decomposition into Cartesian components
    !    u_{perp} = vth0 * int ... imag J1 v_perp dv^3
    !    ux = - ky/|k| * u_{perp}
    !    uy =   kx/|k| * u_{perp}
    ! </doc>
    use dist_fn_arrays, only: aj1, vperp2, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    use kgrids, only: nakx, naky, aky, akx, kinv
    use constants, only: zi
    use agk_mem, only: alloc8, dealloc8
    complex, intent (out), optional :: ux(-ntgrid:,:,:,:)
    complex, intent (out), optional :: uy(-ntgrid:,:,:,:)
    complex, allocatable :: wk1(:,:,:,:)
    integer :: isgn, iglo, is
    integer :: it, ik

    if (.not. (present(ux) .or. present(uy))) return

    g0 = 0. ! g0 is common in this module

    ! we know the size of the array
    if (.not.allocated(wk1)) then
       allocate (wk1(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=wk1,v="wk1")
    end if
    wk1 = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = zi*aj1(iglo)*sqrt(vperp2(iglo)) &
               & *gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, wk1, 1)

    if (present(ux)) then
       ux = 0.
       do it=1,nakx
          do ik=1,naky
             ux(:,it,ik,:)=-aky(ik) * kinv(it,ik) * wk1(:,it,ik,:)
          end do
       end do
       ! clear meaningless zero-zero mode only in nonlinear case
       if (.not. (nakx .eq. 1 .and. naky .eq. 1)) then
          ux(:,1,1,:) = 0.
       end if
       do is=1,nspec
          ux(:,:,:,is)=ux(:,:,:,is)*spec(is)%stm
       end do
    end if

    if (present(uy)) then
       uy = 0.
       do it=1,nakx
          do ik=1,naky
             uy(:,it,ik,:)= akx(it) * kinv(it,ik) * wk1(:,it,ik,:)
          end do
       end do
       ! clear meaningless zero-zero mode
       uy(:,1,1,:) = 0.
       do is=1,nspec
          uy(:,:,:,is)=uy(:,:,:,is)*spec(is)%stm
       end do
    end if

    if (allocated(wk1)) call dealloc8(c4=wk1,v="wk1")
    
  end subroutine getmoms_r_1perp

  subroutine getmoms_r_2parapara (pzz)
    ! <doc>
    !  second order para-para moment at particle position (small r)
    !    Pzz = (n0*T0) * int ... 2 J0 v_para^2 dv^3
    ! </doc>
    use dist_fn_arrays, only: aj0, vpa, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    complex, intent (out) :: pzz(-ntgrid:,:,:,:)
    integer :: isgn, iglo, is
    g0 = 0. ! g0 is common in this module
    pzz = 0.
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*aj0(iglo)*vpa(isgn,iglo)**2*gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, pzz, 1)
    do is=1,nspec
       pzz(:,:,:,is)= pzz(:,:,:,is)*spec(is)%dens*spec(is)%temp
    end do
  end subroutine getmoms_r_2parapara

  subroutine getmoms_r_2perptheta (pxx,pyy,pxy,ppp,ptt)
    ! <doc>
    !  second order perpendicular moment at particle position (small r)
    !  and decomposition into Cartesian components
    !    P_{perp,perp}   = (n0*T0) * int ... 2 J1/alpha v_perp^2 dv^3
    !    P_{theta,theta} = (n0*T0) * int ... 2 (J0 - J1/alpha) v_perp^2 dv^3
    !    Pxx = (kx/|k|)^2 * P_{perp,perp} + (ky/|k|)^2 * P_{theta,theta}
    !    Pyy = (ky/|k|)^2 * P_{perp,perp} + (kx/|k|)^2 * P_{theta,theta}
    !    Pxy = kx*ky/|k|^2 * (P_{perp,perp} - P_{theta,theta})
    ! </doc>
    use dist_fn_arrays, only: aj0, aj1vp2, vperp2, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    use kgrids, only: nakx, naky, aky, akx, kinv
    use agk_mem, only: alloc8, dealloc8
    complex, intent (out), optional :: pxx(-ntgrid:,:,:,:)
    complex, intent (out), optional :: pyy(-ntgrid:,:,:,:)
    complex, intent (out), optional :: pxy(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ppp(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ptt(-ntgrid:,:,:,:)
    complex, allocatable :: wk1(:,:,:,:), wk2(:,:,:,:)
    integer :: isgn, iglo, is
    integer :: it, ik

    if (.not. (present(pxx) .or. present(pyy) .or. present(pxy) &
         .or. present(ppp) .or. present(ptt))) return

    g0 = 0. ! g0 is common in this module

    ! we know the size of the array
    if (.not.allocated(wk1)) then
       allocate (wk1(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=wk1,v="wk1")
    end if
    if (.not.allocated(wk2)) then
       allocate (wk2(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=wk2,v="wk2")
    end if
    wk1 = 0. ; wk2 = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = aj1vp2(iglo)*gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, wk1, 1) ! perp-perp component
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = (2.*aj0(iglo)*vperp2(iglo)-aj1vp2(iglo)) & 
               & *gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, wk2, 1) ! theta-theta component

    if (present(pxx)) then
       pxx = 0.
       do it=1,nakx
          do ik=1,naky
             pxx(:,it,ik,:)=(wk1(:,it,ik,:)*akx(it)**2 + &
                  & wk2(:,it,ik,:)*aky(ik)**2)*kinv(it,ik)**2
          end do
       end do
       ! clear meaningless zero-zero mode
       pxx(:,1,1,:) = 0.
       do is=1,nspec
          pxx(:,:,:,is)= pxx(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if

    if (present(pyy)) then
       pyy = 0.
       do it=1,nakx
          do ik=1,naky
             pyy(:,it,ik,:)=(wk1(:,it,ik,:)*aky(ik)**2 + &
                  & wk2(:,it,ik,:)*akx(it)**2)*kinv(it,ik)**2
          end do
       end do
       ! clear meaningless zero-zero mode
       pyy(:,1,1,:) = 0.
       do is=1,nspec
          pyy(:,:,:,is)= pyy(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if

    if (present(pxy)) then
       pxy = 0.
       do it=1,nakx
          do ik=1,naky
             pxy(:,it,ik,:)=(wk1(:,it,ik,:)-wk2(:,it,ik,:)) &
                  & *akx(it)*aky(ik)*kinv(it,ik)**2
          end do
       end do
       ! clear meaningless zero-zero mode
       pxy(:,1,1,:) = 0.
       do is=1,nspec
          pxy(:,:,:,is)= pxy(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if

    if (present(ppp)) then
       ppp = 0.
       do is=1,nspec
          ppp(:,:,:,is)=wk1(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if

    if (present(ptt)) then
       ptt = 0.
       do is=1,nspec
          ptt(:,:,:,is)=wk2(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if

    if (allocated(wk1)) call dealloc8(c4=wk1,v="wk1")
    if (allocated(wk2)) call dealloc8(c4=wk2,v="wk2")

  end subroutine getmoms_r_2perptheta

  subroutine getmoms_r_2thetapara (pyz, pzx)
    ! <doc>
    !  second order theta-para moment at particle position (small r)
    !  and decomposition into Cartesian components
    !    P_{theta,para} = (n0*T0) * int ... imag 2 J1 v_par dv^3
    !    Pyz =   kx/|k| * P_{theta,para}
    !    Pzx = - ky/|k| * P_{theta,para}
    ! </doc>
    use dist_fn_arrays, only: aj1, vpa, vperp2, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    use kgrids, only: nakx, naky, aky, akx, kinv
    use constants, only: zi
    use agk_mem, only: alloc8, dealloc8
    complex, intent (out), optional :: pyz(-ntgrid:,:,:,:)
    complex, intent (out), optional :: pzx(-ntgrid:,:,:,:)
    complex, allocatable :: wk1(:,:,:,:)
    integer :: isgn, iglo, is
    integer :: it, ik
    
    if (.not. (present(pyz) .or. present(pzx))) return

    g0 = 0. ! g0 is common in this module

    ! we know the size of the array
    if (.not.allocated(wk1)) then
       allocate (wk1(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=wk1,v="wk1")
    end if
    wk1 = 0.

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g0(:,isgn,iglo) = 2.*zi*aj1(iglo) &
               & *sqrt(vperp2(iglo))*vpa(isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do
    call integrate_moment (g0, wk1, 1)

    if (present(pyz)) then
       pyz = 0.
       do it=1,nakx
          do ik=1,naky
             pyz(:,it,ik,:)= akx(it)*kinv(it,ik)*wk1(:,it,ik,:)
          end do
       end do
       ! clear meaningless zero-zero mode
       pyz(:,1,1,:) = 0.
       do is=1,nspec
          pyz(:,:,:,is)= pyz(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if
    
    if (present(pzx)) then
       pzx = 0.
       do it=1,nakx
          do ik=1,naky
             pzx(:,it,ik,:)=-aky(ik)*kinv(it,ik)*wk1(:,it,ik,:)
          end do
       end do
       ! clear meaningless zero-zero mode
       pzx(:,1,1,:) = 0.
       do is=1,nspec
          pzx(:,:,:,is)= pzx(:,:,:,is)*spec(is)%dens*spec(is)%temp
       end do
    end if

    if (allocated(wk1)) call dealloc8(c4=wk1,v="wk1")
    
  end subroutine getmoms_r_2thetapara

  subroutine getmoms_r_3 (qzzz, qzpp)
    ! <doc>
    !  third order moments at particle position (small r)
    !    Qzzz = (n0*T0*vth0) * int ... 2 J0 v_para^3 dv^3
    ! </doc>
    use dist_fn_arrays, only: aj0, vpa, vperp2, gnew
    use agk_layouts, only: g_lo
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    complex, intent (out), optional :: qzzz(-ntgrid:,:,:,:), qzpp(-ntgrid:,:,:,:)
    integer :: isgn, iglo, is
    
    if (.not. (present(qzzz) .or. present(qzpp))) return
    
    if (present(qzzz)) then
       qzzz = 0.
       g0 = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(iglo)*vpa(isgn,iglo)**3*gnew(:,isgn,iglo)
          end do
       end do
       call integrate_moment (g0, qzzz, 1)
       do is=1,nspec
          qzzz(:,:,:,is)= qzzz(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%stm
       end do
    endif

    if (present(qzpp)) then
       qzpp = 0.
       g0 = 0.
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(iglo)*vpa(isgn,iglo)*vperp2(iglo)*gnew(:,isgn,iglo)
          end do
       end do
       call integrate_moment (g0, qzpp, 1)
       do is=1,nspec
          qzpp(:,:,:,is)= qzpp(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%stm
       end do
    endif
  end subroutine getmoms_r_3


  subroutine getmoms_r_temp (tpara,tperp,adjust_in)
    ! <doc>
    !  Get temperatures T_{perp,perp} and T_{para,para} by solving
    !
    !   P_{para,para}/(n0*T0) = 
    !     a11 * n/n0 + a12 * T_{para,para}/T0 + a13 * T_{perp,perp}/T0
    !   P_{perp,perp}/(n0*T0) = 
    !     a21 * n/n0 + a22 * T_{para,para}/T0 + a23 * T_{perp,perp}/T0
    !
    !  where coefficients are
    !   a11 = mom_coeff_npara = 1
    !   a12 = 1
    !   a13 = mom_coeff_tperp = 0
    !   a21 = mom_coeff_npara = 1
    !   a22 = mom_coeff_tpara = 0
    !   a23 = 1
    !  Note that mom_coeff_... are moments of a Maxwellian whose analytical
    !  values are given at right most equality. See [[init_mom_coeff]].
    ! </doc>
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use kgrids, only: nakx, naky
    use agk_mem, only: alloc8, dealloc8

    complex, intent (out) :: tperp(-ntgrid:,:,:,:), tpara(-ntgrid:,:,:,:)

    complex, allocatable :: dens(:,:,:,:), ppp(:,:,:,:), pzz(:,:,:,:) ! local

    real :: a, b, tpara2, tperp2
    integer :: ig, it, ik, is

    logical, intent (in), optional :: adjust_in
    logical :: adjust_loc = .true.

    adjust_loc=.true.
    if(present(adjust_in)) adjust_loc = adjust_in

    ! we know the size of the array
    if (.not.allocated(dens)) then
       allocate (dens(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=dens,v="dens")
    end if
    dens = 0.
    if (.not.allocated(ppp)) then
       allocate (ppp(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=ppp,v="ppp")
    end if
    ppp = 0.
    if (.not.allocated(pzz)) then
       allocate (pzz(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=pzz,v="pzz")
    end if
    pzz = 0.

    call getmoms_r (dens=dens,pzz=pzz,ppp=ppp, adjust_in=adjust_loc)

    do ig=-ntgrid,ntgrid
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                tpara2=pzz(ig,it,ik,is)/(spec(is)%dens*spec(is)%temp) &
                     & - dens(ig,it,ik,is)/spec(is)%dens &
                     &   * mom_coeff_npara(it,ik,is)
                tperp2=ppp(ig,it,ik,is)/(spec(is)%dens*spec(is)%temp) &
                     & - dens(ig,it,ik,is)/spec(is)%dens &
                     &   * mom_coeff_nperp(it,ik,is)

                a=mom_coeff_tperp(it,ik,is)
                b=mom_coeff_tpara(it,ik,is)

                tpara(ig,it,ik,is)=(   tpara2-a*tperp2)/(1.-a*b)
                tperp(ig,it,ik,is)=(-b*tpara2+  tperp2)/(1.-a*b)
             end do
          end do
       end do
    end do

    do is=1,nspec
       tpara(:,:,:,is)=tpara(:,:,:,is)*spec(is)%temp
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
    end do

    if (allocated(dens)) call dealloc8(c4=dens,v="dens")
    if (allocated(ppp)) call dealloc8(c4=ppp,v="ppp")
    if (allocated(pzz)) call dealloc8(c4=pzz,v="pzz")

  end subroutine getmoms_r_temp
    
  subroutine init_fieldeq
    use dist_fn_arrays, only: aj0, aj1vp2
    use species, only: nspec, spec, has_electron_species
    use theta_grid, only: ntgrid
    use kgrids, only: naky, nakx, aky, kperp2
    use le_grids, only: integrate_species
    use agk_layouts, only: g_lo
    use run_parameters, only: tite
    use agk_mem, only: alloc8
    implicit none
    integer :: iglo, isgn
    integer :: ik, it
    complex, dimension (nakx,naky) :: tot
    real, dimension (nspec) :: wgt
    logical :: done = .false.


    if (.not. done) then 
      allocate (gamtot(nakx,naky)); call alloc8(r2=gamtot,v="gamtot")
      allocate (gamtot1(nakx,naky)); call alloc8(r2=gamtot1,v="gamtot1")
      allocate (gamtot2(nakx,naky)); call alloc8(r2=gamtot2,v="gamtot2")
      if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         allocate (gamtot3(nakx,naky)); call alloc8(r2=gamtot3,v="gamtot3")
      end if
    end if 
    done = .true.
    !               Z^2 n_{0s}
    ! gamtot = sum ------------ [1 - Gamma_0(alpha_s)]
    !           s     T_{0s}
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       gk0(:,iglo) = 1.0 - aj0(iglo)**2
    end do
    wgt = spec%z*spec%z * spec%dens / spec%temp
    call integrate_species (gk0, wgt, tot)
    do ik =1, naky
       do it = 1, nakx
          gamtot(it,ik) = real(tot(it,ik)) + kperp2(it,ik)*poisfac
       end do
    end do
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       gk0(:,iglo) = aj0(iglo)*aj1vp2(iglo)
    end do
    wgt = spec%z*spec%dens
    call integrate_species (gk0, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       gk0(:,iglo) = 0.5*aj1vp2(iglo)**2
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (gk0, wgt, tot)
    gamtot2 = real(tot)

! adiabatic electrons 
    if (.not. has_electron_species(spec)) then
       select case (adiabatic_option_switch)
       case (adiabatic_option_yavg)
          do ik = 1, naky
             if (aky(ik) > epsilon(0.0)) gamtot(:,ik) = gamtot(:,ik) + tite
          end do
       case (adiabatic_option_fieldlineavg)
          gamtot  = gamtot + tite
          gamtot3 = (gamtot-tite) / gamtot
          where (gamtot3 < 2.*epsilon(0.0)) gamtot3 = 1.0
       case (adiabatic_option_zero)
          ! nothing to add
       case default
          gamtot = gamtot + tite
       end select
    endif
  end subroutine init_fieldeq

  subroutine getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid, delthet, jacob
    use kgrids, only: naky, nakx, aky, kperp2
    use run_parameters, only: use_Phi, use_Apar, use_Bpar
    use run_parameters, only: beta, tite
    use species, only: spec, has_electron_species
    use agk_mem, only: alloc8
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    real, allocatable, dimension(:,:), save :: fl_avg, awgt
    integer :: ik, it
    logical :: first = .true.
    
    if (first) then
       allocate (fl_avg(nakx, naky)); call alloc8(r2=fl_avg,v="fl_avg")
    end if
    fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          
          if (first) then 
             allocate (awgt(nakx, naky)); call alloc8(r2=awgt,v="awgt")
             awgt = 0.
             do ik = 1, naky
                if (aky(ik) > epsilon(0.0)) cycle
                do it = 1, nakx
                   awgt(it,ik) = tite/(gamtot(it,ik)*gamtot3(it,ik))
                end do
             end do
          endif
           
          do ik = 1, naky
             do it = 1, nakx
                fl_avg(it,ik) = sum(antot(:,it,ik))*awgt(it,ik)
             end do
          end do

       end if
    end if

    if (use_Phi) then
       do ik = 1, naky
          do it = 1, nakx
             fieldeq(:,it,ik) = antot(:,it,ik) &
                  + bpar(:,it,ik)*gamtot1(it,ik) - gamtot(it,ik)*phi(:,it,ik) 
          end do
       end do
             
       if (.not. has_electron_species(spec)) then
          do ik = 1, naky
             do it = 1, nakx
                fieldeq(:,it,ik) = fieldeq(:,it,ik) + fl_avg(it,ik)
             end do
          end do
       end if
    end if

    if (use_Apar) then
       do ik = 1, naky
          do it = 1, nakx
             fieldeqa(:,it,ik) = antota(:,it,ik) - kperp2(it,ik)*apar(:,it,ik)
          end do
       end do
    end if

    if (use_Bpar) then
       do ik =1, naky
          do it = 1, nakx
             fieldeqp(:,it,ik) = (antotp(:,it,ik) &
                  + bpar(:,it,ik)*gamtot2(it,ik)+0.5*phi(:,it,ik)*gamtot1(it,ik))*beta
          end do
       end do

       do ik = 1, naky
          do it = 1, nakx
             fieldeqp(:,it,ik) = fieldeqp(:,it,ik) + bpar(:,it,ik)
          end do
       end do
    end if

    first = .false.

  end subroutine getfieldeq1

  subroutine getfieldeq (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp)
    use theta_grid, only: ntgrid
    use kgrids, only: naky, nakx
    use agk_mem, only: alloc8, dealloc8
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (out) ::fieldeq,fieldeqa,fieldeqp
    complex, dimension (:,:,:), allocatable :: antot, antota, antotp

    allocate (antot (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=antot,  v="antot")
    allocate (antota(-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=antota, v="antota")
    allocate (antotp(-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=antotp, v="antotp")

    call getan (antot, antota, antotp)
!!!    call getan2 (antot, antota, antotp) !RN
    call getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp)

    call dealloc8 (c3=antot,  v="antot")
    call dealloc8 (c3=antota, v="antota")
    call dealloc8 (c3=antotp, v="antotp")

  end subroutine getfieldeq

! TT: Given initial distribution function this obtains consistent fields
  subroutine get_init_field (phi, apar, bpar)
    ! inverts the field equations:
    !   gamtot * phi - gamtot1 * bpar = antot
    !   kperp2 * apar = antota
    !   beta/2 * gamtot1 * phi + (beta * gamtot2 + 1) * bpar = - beta * antotp
    ! I haven't made any check for use_Bpar=T case.
    use run_parameters, only: beta, use_Phi, use_Apar, use_Bpar
    use theta_grid, only: ntgrid
    use kgrids, only: nakx, naky, kperp2
    use species, only: nspec
    use dist_fn_arrays, only: aj0, vpa
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
    real, dimension (nspec) :: wgt
    real, dimension (-ntgrid:ntgrid,nakx,naky) :: denominator
    complex, dimension (-ntgrid:ntgrid,nakx,naky) :: antot, antota, antotp
    complex, dimension (-ntgrid:ntgrid,nakx,naky) :: numerator

    antot=0.0 ; antota=0.0 ; antotp=0.0
    call getan (antot, antota, antotp)
!!!    call getan2 (antot, antota, antotp, adjust_in=.false.) !RN

    ! get phi
    if (use_Phi) then
       numerator = spread(beta * gamtot2 + 1.0, 1, ntgrid*2+1) * antot &
            & - spread(beta * gamtot1, 1, ntgrid*2+1) * antotp
       denominator = spread( (beta * gamtot2 + 1.0) * gamtot &
            & + (beta/2.0) * gamtot1 * gamtot1, 1, ntgrid*2+1)
       where (abs(denominator) < epsilon(0.0)) ! it == ik == 1 only
          phi = 0.0
       elsewhere
          phi = numerator / denominator
       end where
    end if

    ! get apar
    if (use_Apar) then
       denominator = spread(kperp2,1,ntgrid*2+1)
       where (abs(denominator) < epsilon(0.0)) ! it == ik == 1 only
          apar = 0.0
       elsewhere
          apar = antota / denominator
       end where
    end if

    ! get bpar
    if (use_Bpar) then
       numerator = - spread(beta * gamtot, 1, ntgrid*2+1) * antotp &
            & - spread((beta/2.0) * gamtot1, 1, ntgrid*2+1) * antot
       ! following line is actually same with the denom for phi
       denominator = spread(gamtot * (beta * gamtot2 + 1.0) &
            & + (beta/2.0) * gamtot1 * gamtot1, 1, ntgrid*2+1)
       where (abs(denominator) < epsilon(0.0)) ! it == ik == 1 only
          bpar = 0.0
       elsewhere
          bpar = numerator / denominator
       end where
    end if
    
  end subroutine get_init_field

  subroutine get_init_field_h (phi, apar, bpar)
    ! calculate fields from the moments of h (not g)
    ! antot, antota, antotp are calcuated for h (not g)
    ! this routine should only be used in recon3 in init_g
    ! inverts the field equations:
    !   phi = antot / sum_s(q_{s}^2 n_{0s}/T_{0s})
    !   apar = antota / kperp2
    !   bpar = - beta * antotp
    use run_parameters, only: beta, use_Phi, use_Apar, use_Bpar
    use species, only: nspec,spec
    use theta_grid, only: ntgrid
    use kgrids, only: kperp2, nakx, naky
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
    complex, dimension (-ntgrid:ntgrid,nakx,naky) :: antot, antota, antotp

    antot=0.0 ; antota=0.0 ; antotp=0.0

    ! moments of h
    call getan2 (antot, antota, antotp, adjust_in=.true.)

    ! get phi
    if (use_Phi) then
       phi = antot / sum(spec(1:nspec)%z**2*spec(1:nspec)%dens/spec(1:nspec)%temp)
    end if

    ! get apar
    if (use_Apar) then
       where (abs(spread(kperp2,1,ntgrid*2+1)) < epsilon(0.0))
          apar = 0.0
       elsewhere
          apar = antota / spread(kperp2,1,ntgrid*2+1)
       end where
    end if

    ! get bpar
     if (use_Bpar) then
       bpar = - beta*antotp
    end if
 
  end subroutine get_init_field_h

  subroutine fields2g (phi,apar,bpar,ispec,g)
    ! <doc>
    !  determines g from the given fields assuming the fields are
    !  given by g of one species.
    ! </doc>
    use agk_layouts, only: g_lo, it_idx, ik_idx
    use run_parameters, only: beta
    use species, only: spec
    use theta_grid, only: ntgrid
    use kgrids, only: nakx, naky, kperp2
    use dist_fn_arrays, only: vpa, vperp2
    use agk_mem, only: alloc8, dealloc8
    complex, intent(in) :: &
         phi(-ntgrid:,:,:), apar(-ntgrid:,:,:), bpar(-ntgrid:,:,:)
    integer, intent(in) :: ispec
    complex, intent(out) :: g(:,:,:)

    complex, allocatable :: dens(:,:,:), upar(:,:,:), tperp(:,:,:)
    
    real :: b
    integer :: ig, iglo, isgn, it, ik

    if (.not.allocated(dens)) then
       allocate(dens(-ntgrid:ntgrid,nakx,naky)); call alloc8(c3=dens,v='dens')
    endif
    if (.not.allocated(upar)) then
       allocate(upar(-ntgrid:ntgrid,nakx,naky)); call alloc8(c3=upar,v='upar')
    endif
    if (.not.allocated(dens)) then
       allocate(tperp(-ntgrid:ntgrid,nakx,naky)); call alloc8(c3=tperp,v='tperp')
    endif

    do ig=-ntgrid, ntgrid
       upar(ig,:,:)=apar(ig,:,:)*.5*kperp2 &
            & /(beta*spec(ispec)%z*spec(ispec)%dens*spec(ispec)%stm)
    end do

    do isgn=1,2
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          g(:,isgn,iglo)=exp(b**2*.25)*(dens(:,it,ik)+2.*vpa(isgn,ig)*upar(:,it,ik)+tperp(:,it,ik)*(vperp2(iglo)-1.+b**2*.25))
       end do
    end do
    
    if (allocated(dens)) call dealloc8(c3=dens,v='dens')
    if (allocated(upar)) call dealloc8(c3=upar,v='upar')
    if (allocated(tperp)) call dealloc8(c3=tperp,v='tperp')

  end subroutine fields2g

  subroutine flux (phi, apar, bpar, &
        pflux,  qflux,  vflux, &
       pmflux, qmflux, vmflux, &
       pbflux, qbflux, vbflux)
    use species, only: spec
    use theta_grid, only: ntgrid, gradpar
    use kgrids, only: naky, nakx
    use le_grids, only: e
    use dist_fn_arrays, only: g, aj0, vpac, vpa, aj1vp2, vperp2
    use agk_layouts, only: g_lo, ie_idx, is_idx
    use mp, only: proc0
    use run_parameters, only: use_Phi, use_Apar, use_Bpar
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:,:), intent (out) :: pflux, pmflux, pbflux
    real, dimension (:,:,:), intent (out) :: vflux, vmflux, vbflux
    real, dimension (:,:,:,:), intent (out) :: qflux, qmflux, qbflux
    real :: anorm
    integer :: ig, it, ik, is, isgn
    integer :: iglo

    if (proc0) then
       pflux = 0.0;   qflux = 0.0;   vflux = 0.0
       pmflux = 0.0;  qmflux = 0.0;  vmflux = 0.0
       pbflux = 0.0;  qbflux = 0.0;  vbflux = 0.0
    end if
    
    anorm = sum(cabs(phi) + cabs(apar) + cabs(bpar))
    if (anorm < epsilon(0.0)) return


    if (use_Phi) then
       do isgn = 1, 2
          do ig=-ntgrid,ntgrid
             g0(ig,isgn,:) = g(ig,isgn,:)*aj0
          end do
       end do
       call get_flux (phi, pflux)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call get_flux (phi, qflux(:,:,:,1))

       do isgn = 1, 2
          do ig=-ntgrid,ntgrid
             g0(ig,isgn,:) = g(ig,isgn,:)*2.*vpa(isgn,:)**2*aj0
          end do
       end do
       call get_flux (phi, qflux(:,:,:,2))

       do isgn = 1, 2
          do ig=-ntgrid,ntgrid
             g0(ig,isgn,:) = g(ig,isgn,:)*vperp2*aj0
          end do
       end do
       call get_flux (phi, qflux(:,:,:,3))

       do isgn = 1, 2
          do ig=-ntgrid,ntgrid
             g0(ig,isgn,:) = g(ig,isgn,:)*aj0*vpac(isgn,:)
          end do
       end do
       call get_flux (phi, vflux)

    else
       pflux = 0.
       qflux = 0.
       vflux = 0.
    end if

    if (use_Apar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(iglo)*spec(is)%stm*vpa(isgn,iglo)
          end do
       end do
       call get_flux (apar, pmflux)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call get_flux (apar, qmflux(:,:,:,1))
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(iglo)*spec(is)%stm*vpa(isgn,iglo) &
                  *2.*vpa(isgn,iglo)**2
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,2))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(iglo)*spec(is)%stm*vpa(isgn,iglo) &
                  *vperp2(iglo)
          end do
       end do
       call get_flux (apar, qmflux(:,:,:,3))
       
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = -g(:,isgn,iglo)*aj0(iglo)*spec(is)%stm*vpa(isgn,iglo)*vpac(isgn,iglo)
          end do
       end do
       call get_flux (apar, vmflux)
    else
       pmflux = 0.
       qmflux = 0.
       vmflux = 0.
    end if

    if (use_Bpar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) = g(:,isgn,iglo)*aj1vp2(iglo)*spec(is)%tz
          end do
       end do
       call get_flux (bpar, pbflux)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g0(:,:,iglo) = g0(:,:,iglo)*e(ie_idx(g_lo,iglo), is_idx(g_lo,iglo))
       end do
       call get_flux (bpar, qbflux(:,:,:,1))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) &
                  = g(:,isgn,iglo)*aj1vp2(iglo)*spec(is)%tz*2.*vpa(isgn,iglo)**2
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,2))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) = g(:,isgn,iglo)*aj1vp2(iglo)*spec(is)%tz*vperp2(iglo)
          end do
       end do
       call get_flux (bpar, qbflux(:,:,:,3))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g0(:,isgn,iglo) = g(:,isgn,iglo)*aj1vp2(iglo)*spec(is)%tz*vpac(isgn,iglo)
          end do
       end do
       call get_flux (bpar, vbflux)
    else
       pbflux = 0.
       qbflux = 0.
       vbflux = 0.
    end if

  end subroutine flux

  subroutine get_flux (fld, flx)
    use theta_grid, only: ntgrid
    use kgrids, only: nakx, aky, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0
    use agk_mem, only: alloc8, dealloc8
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (in out) :: flx
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is, ig

    allocate (total (-ntgrid:ntgrid, nakx, naky, nspec)) ; call alloc8 (c4=total, v="total")
    call integrate_moment (g0, total)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, nakx
                flx(it,ik,is) = sum(aimag(total(:,it,ik,is)*conjg(fld(:,it,ik)))) &
                     *aky(ik)
             end do
          end do
       end do

       flx = flx*0.5

    end if

    call dealloc8 (c4=total,v="total")

  end subroutine get_flux
!>GGH
!-----------------------------------------------------------------------------
! Density: Calculate Density perturbations
!-----------------------------------------------------------------------------
  subroutine get_dens_vel (dv, dvk, phi, apar, bpar, phinew, aparnew, bparnew)
    use constants, only: pi
    use dist_fn_arrays, only: aj0, vpac, vperp2, g, gnew
    use agk_heating, only: dens_vel_diagnostics
    use agk_layouts, only: g_lo, ik_idx, it_idx, is_idx
    use kgrids, only: nakx, naky,aky
    use le_grids, only: integrate_moment
    use mp, only: proc0
    use nonlinear_terms, only: nonlin
    use run_parameters, only: fphi, fbpar
    use species, only: spec,nspec
    use theta_grid, only: jacob, delthet, ntgrid
    use agk_mem, only: alloc8, dealloc8 
    implicit none
    !Passed
    type (dens_vel_diagnostics) :: dv
    type (dens_vel_diagnostics), dimension(:,:) :: dvk
    complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar, phinew, aparnew, bparnew
    !Local 
    integer :: isgn, iglo, ig, is, ik, it            !Indices
    complex :: j0phiavg                              !J_0 q phi/T 
    complex :: havg                                  !Time and z centered h
    complex :: phiavg                               !Time and z centered phi
    
    complex, dimension(:,:,:,:), allocatable :: tot  !Integrated DF
    complex, dimension(:,:,:,:), allocatable :: tot2  !Integrated DF
    real :: fac2                                     !Factor
    real, dimension (:), allocatable :: wgt

    !Allocate tot variable for each calculation
    allocate (tot (-ntgrid:ntgrid, nakx, naky, nspec)) ; call alloc8 (c4=tot, v="tot")

    !Set up weighting factors for z-sums on proc0
    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid)); call alloc8(r1=wgt,v="wgt")
       wgt = 0.
       do ig=-ntgrid,ntgrid
          wgt(ig) = delthet*jacob
       end do
       wgt = wgt/sum(wgt)         
    endif

    !CONVERT FROM g TO h
    call g_adjust (g,    phi,    bpar,    fphi, fbpar)
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar, is_g_gnew)

    !Parallel velocity perturbation--------------------------------------------
    !Loop through and constuct time and z-centered integrands
!    do iglo=g_lo%llim_proc, g_lo%ulim_proc
!       is = is_idx(g_lo, iglo)
!       it = it_idx(g_lo, iglo)
!       ik = ik_idx(g_lo, iglo)
!       if (nonlin .and. it == 1 .and. ik == 1) cycle
!       do isgn=1,2
!          do ig=-ntgrid, ntgrid-1
!
!             !Time and z-centered h
!             havg = favg (g   (ig  ,isgn,iglo), &
!                          g   (ig+1,isgn,iglo), &
!                          gnew(ig  ,isgn,iglo), &
!                          gnew(ig+1,isgn,iglo))
!
!             !NOTE: Do I need vpac here?
!             g0(ig,isgn,iglo) = vpac(isgn,iglo)*havg
!
!          end do
!       end do
!    end do
!
!    !Integrate over velocity
!    call integrate_moment (g0, tot)
!    
!    !Sum over theta (z-direction) to get perpendicular Fourier components only
!    if (proc0) then
!       do ik = 1, naky
!          fac2 = 0.5
!          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!          do it = 1, nakx
!             !Skip mean value (k=0) components
!             if (nonlin .and. it == 1 .and. ik == 1) cycle
!             do is = 1, nspec
!                do ig = -ntgrid, ntgrid-1
!                   dvk(it,ik)%dvpar(is) = dvk(it,ik)%dvpar(is) &
!                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                end do
!                dv % dvpar(is) = dv%dvpar(is) + dvk(it,ik)%dvpar(is)
!             end do
!          end do
!       end do
!    end if
!
!    !Perpendicular velocity perturbation---------------------------------------
!    !Loop through and constuct time and z-centered integrands
!    !Non-Boltzmann part
!    do iglo=g_lo%llim_proc, g_lo%ulim_proc
!       is = is_idx(g_lo, iglo)
!       it = it_idx(g_lo, iglo)
!       ik = ik_idx(g_lo, iglo)
!       do isgn=1,2
!          do ig=-ntgrid, ntgrid-1
!
!             !Time and z-centered h
!             havg = favg (g   (ig  ,isgn,iglo), &
!                          g   (ig+1,isgn,iglo), &
!                          gnew(ig  ,isgn,iglo), &
!                          gnew(ig+1,isgn,iglo))
!
!             g0(ig,isgn,iglo) = sqrt(vperp2(iglo))*havg
!
!          end do
!       end do
!    end do
!
!    !Integrate over velocity
!    call integrate_moment (g0, tot)
!    
!    !Add in adiabatic part 
!    do iglo=g_lo%llim_proc, g_lo%ulim_proc
!       is = is_idx(g_lo, iglo)
!       it = it_idx(g_lo, iglo)
!       ik = ik_idx(g_lo, iglo)
!       do ig= -ntgrid, ntgrid-1
!          !Time and z-center J_0 Phi
!          j0phi= 0.25*spec(is)%zt* &
!               ( aj0(iglo) *(phi(ig  ,it,ik)+ phinew(ig,  it,ik) ) + &
!                 aj0(iglo) *(phi(ig+1,it,ik)+ phinew(ig+1,it,ik)))
!
!          tot(ig,it,ik,is)=tot(ig,it,ik,is)+ sqrt(pi/2.)*(1.-j0phi)               
!       enddo
!    enddo
!
!    !Sum over theta (z-direction) to get perpendicular Fourier components only
!    if (proc0) then
!       do ik = 1, naky
!          fac2 = 0.5
!          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!           do it = 1, nakx
!             !Skip mean value (k=0) components
!             if (nonlin .and. it == 1 .and. ik == 1) cycle
!             do is = 1, nspec
!                do ig = -ntgrid, ntgrid-1
!                   dvk(it,ik)%dvperp(is) = dvk(it,ik)%dvperp(is) &
!                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                end do
!                dv % dvperp(is) = dv%dvperp(is) + dvk(it,ik)%dvperp(is)
!             end do
!          end do
!       end do
!    end if
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
    !Density perturbation------------------------------------------------------
    !Loop through and constuct time and z-centered integrands
    !Non-Boltzmann part
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid, ntgrid

             havg = gnew(ig,isgn,iglo)

             !Adiabatic part
             phiavg = phinew(ig+1,it,ik) * fphi * spec(is)%zt

             g0(ig,isgn,iglo) = aj0(iglo) * havg - phiavg * spec(is)%dens

          end do
       end do
    end do

    !Integrate over velocity
    call integrate_moment (g0, tot)
!    tot2=0.

    !Calculate Boltzmann part 
!    if (proc0) then
!       do ik=1,naky
!          fac2 = 0.5
!          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!          do it = 1,nakx
!             if (nonlin .and. it == 1 .and. ik == 1) cycle
!             do is=1,nspec
!                iglo=iglo_idx
!                do ig=-ntgrid, ntgrid-1
!                   
!                   phi_avg = favg (phi   (ig  ,it,ik), &
!                        phi   (ig+1,it,ik), &
!                        phinew(ig  ,it,ik), &
!                        phinew(ig+1,it,ik))
!                   
!                   dvk(it,ik)%dvpar(is) = dvk(it,ik)%dvpar(is) &
!                        + spec(is)%zt*abs(phi_avg)*wgt(ig)*fac2
!
!                end do
!                dv % dvpar(is) = dv%dvpar(is) + dvk(it,ik)%dvpar(is)
!             enddo
!          end do
!       end do
!    endif

    !Sum over theta (z-direction) to get perpendicular Fourier components only
    if (proc0) then
!       dv%dvpar(:)=0.  !Initialize
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, nakx
             !Skip mean value (k=0) components
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
!                dvk(it,ik)%dvpar(is)=0. !Initialize
                do ig = -ntgrid, ntgrid
!                   dvk(it,ik)%dn(is) = dvk(it,ik)%dn(is) &
!                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
!                   dvk(it,ik)%dn(is) = dvk(it,ik)%dn(is) &
!                        + abs(tot(ig,it,ik,is))*wgt(ig)*fac2
                   dvk(it,ik)%dn(is) = dvk(it,ik)%dn(is) &
                        + abs(tot(ig,it,ik,is))*wgt(ig)*fac2
!                   dvk(it,ik)%dvpar(is) = dvk(it,ik)%dvpar(is) &
!                        + abs(tot2(ig,it,ik,is))*wgt(ig)*fac2
                end do
                dv % dn(is) = dv%dn(is) + dvk(it,ik)%dn(is)
!                dv % dvpar(is) = dv%dvpar(is) + dvk(it,ik)%dvpar(is)
             end do
          end do
       end do
    end if

    !Calculate Boltzmann part 
!    do iglo=g_lo%llim_proc, g_lo%ulim_proc
!       is = is_idx(g_lo, iglo)
!       it = it_idx(g_lo, iglo)
!       ik = ik_idx(g_lo, iglo)
!       if (nonlin .and. it == 1 .and. ik == 1) cycle
!       do isgn=1,2
!          do ig=-ntgrid, ntgrid-1
!             j0phiavg = favg (phi   (ig  ,it,ik), &
!                  phi   (ig+1,it,ik), &  
!                  phinew(ig  ,it,ik), &  
!                  phinew(ig+1,it,ik)) * aj0(iglo) * fphi * spec(is)%zt
!             
!             g0(ig,isgn,iglo) = -j0phiavg*spec(is)%dens
!
!          end do
!       end do
!    end do
!
!    call integrate_moment (g0, tot)
!
!    !Sum over theta (z-direction) to get perpendicular Fourier components only
!    if (proc0) then
!       do ik = 1, naky
!          fac2 = 0.5
!          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
!          do it = 1, nakx
!             !Skip mean value (k=0) components
!             if (nonlin .and. it == 1 .and. ik == 1) cycle
!             do is = 1, nspec
!                do ig = -ntgrid, ntgrid-1
!                   dvk(it,ik)%dvpar(is) = dvk(it,ik)%dvpar(is) &
!                        + abs(tot(ig,it,ik,is))*wgt(ig)*fac2
!                end do
!                dv % dvpar(is) = dv%dvpar(is) + dvk(it,ik)%dvpar(is)
!             end do
!          end do
!       end do
!    end if
!
    
    call dealloc8 (c4=tot,v="tot")

    !CONVERT FROM h BACK TO g
    call g_adjust (g,    phi,    bpar,    -fphi, -fbpar)
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar, is_g_gnew)

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  end subroutine get_dens_vel

!-----------------------------------------------------------------------------
! Calculate J_external
!-----------------------------------------------------------------------------
  subroutine get_jext(j_extx, j_exty, j_extz)
    use kgrids, only: nakx, naky, aky, kperp2
    use mp, only: proc0
    use theta_grid, only: jacob, delthet, ntgrid
    use antenna, only: antenna_apar
    use agk_mem, only: alloc8, dealloc8
    implicit none
    !Passed
    real, dimension(:,:) ::   j_extx, j_exty, j_extz
    !Local 
    complex, dimension(:,:,:), allocatable :: j_extxa, j_extya, j_extza
    integer :: ig,ik, it                             !Indices
    real :: fac2                                     !Factor
    real, dimension (:), allocatable :: wgt
    
    !Get j_ext at current time
    allocate (j_extza (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=j_extza, v="j_extza")
    allocate (j_extya (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=j_extya, v="j_extya")
    allocate (j_extxa (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=j_extxa, v="j_extxa")
    j_extxa = 0.; j_extya = 0.; j_extza = 0.

    call antenna_apar (j_extxa,j_extya,j_extza)            

    !Set weighting factor for z-averages

    allocate (wgt(-ntgrid:ntgrid)); call alloc8(r1=wgt,v="wgt")
    !GGH NOTE: Here wgt is 1/(2*ntgrid)
    wgt = 0.
    do ig=-ntgrid,ntgrid
       wgt(ig) = delthet*jacob
    end do
    wgt = wgt/sum(wgt)         
    
    !Take average over z (and weight k modes for later sum?)
    do ig=-ntgrid, ntgrid
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, nakx
             j_extx(it,ik)=j_extx(it,ik)+real(j_extxa(ig,it,ik))*wgt(ig)*fac2
             j_exty(it,ik)=j_exty(it,ik)+real(j_extya(ig,it,ik))*wgt(ig)*fac2
             j_extz(it,ik)=j_extz(it,ik)+real(j_extza(ig,it,ik))*wgt(ig)*fac2
          end do
       end do
    enddo

    call dealloc8 (c3=j_extxa, v="j_extxa")
    call dealloc8 (c3=j_extya, v="j_extya")
    call dealloc8 (c3=j_extza, v="j_extza")
!    deallocate (wgt)
    call dealloc8(r1=wgt,v="wgt")

  end subroutine get_jext
!<GGH
!-----------------------------------------------------------------------------
  subroutine get_heat (hk, phi, apar, bpar, phinew, aparnew, bparnew)
    use mp, only: proc0, iproc
    use constants, only: pi, zi
    use kgrids, only: nakx, naky, aky, akx, kperp2
    use dist_fn_arrays, only: vpa, vpac, aj0, aj1vp2, g, gnew
    use dist_fn_arrays, only: c_rate
    use agk_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
    use le_grids, only: integrate_moment, e
    use species, only: spec, nspec,has_electron_species
    use theta_grid, only: dz, ntgrid, delthet, gradpar
    use run_parameters, only: fphi, fapar, fbpar, beta, tite
    use agk_time, only: dtime
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_apar, a_ext_data
    use hyper, only: D_v, D_eta, nexp
    use agk_heating_index
    use agk_mem, only: alloc8, dealloc8
    use dist_fn_arrays, only: g_eq
    use fields_arrays, only: phi_eq, apar_eq, bpar_eq
    use run_parameters, only: store_eq
    implicit none
    !Passed
    real, dimension(:,:,:), intent (out) :: hk
    complex, dimension (-ntgrid:,:,:) :: phi, apar, bpar, phinew, aparnew, bparnew
    !Local
    complex, dimension(:,:,:,:), allocatable :: tot
    complex, dimension(:,:,:), allocatable :: epar, bpardot, apardot, phidot, a_ext_old, a_ext_new
    complex, dimension(:,:,:), allocatable :: j_extx, j_exty, j_extz
    complex :: fac, chi, havg, phi_m, apar_m, hdot
    complex :: chidot, j0phiavg, j1bparavg, j0aparavg
    complex :: de, denew, phi_avg
! TT: wgt is homogeneous in AGK
    real, dimension (:), allocatable :: wgt
!!$    real :: wgt
! <TT
    real :: dtinv, akperp4,  akperp2
    real :: bpar2, bpar2new, bperp2, bperp2new, fac3
    integer :: isgn, iglo, ig, is, ik, it, ie
    complex :: phidot1,apardot1,bpardot1
    complex, dimension(:,:,:,:), allocatable :: upar, ux, uy
    complex, dimension(:,:,:,:), allocatable :: dens, pxx, pyy, pzz
    real :: bperp_pt2, bpar_pt2
    complex :: phi_pt_avg

    hk=0.

    if(proc0) then
      if(debug) write(*,*) "dist_fn: starting get_heat." !DEBUG-KDN-110718
       !GGH ERROR?: I think this should be weighted over only -ntgrid:ntgrid-1,
       !            since everything is periodic [f(ntgrid)=f(-ntgrid)]   18 MAR 2009
       !GGH ERROR?: Should we be summing or averaging  along z (theta) ?  18 MAR 2009
       ! TT: trying to fix the problem
       if(.not.allocated(wgt)) then
          allocate (wgt(-ntgrid:ntgrid-1)); call alloc8(r1=wgt,v="wgt")
       end if
       ! Same coding with get_vol_average_all in agk_diagnostics.f90
       ! May be combined?
       wgt = dz ! = delthet*jacob ! 2pi * z0 / ntheta
       wgt = wgt/sum(wgt)  ! 1 / (Nz-1)
       ! TT: The following should give the same result for general nperiod
       ! and do not require array
!!$       wgt = 1.0 / (2*ntgrid) ! 1 / (Nz-1)
      if(debug) write(*,*) "dist_fn:get_heat: proc0 allocation finished?" !DEBUG-KDN-110718

    endif

!-----------------------------------------------------------------------------
! Ion/Electron heating------------------------------------------------------
!-----------------------------------------------------------------------------

    allocate ( phidot (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=phidot, v="phidot")
    allocate (apardot (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=apardot, v="apardot")
    allocate (bpardot (-ntgrid:ntgrid, nakx, naky)) ; call alloc8 (c3=bpardot, v="bpardot")

    call dot ( phi,  phinew,  phidot, fphi)
    call dot (apar, aparnew, apardot, fapar)
    call dot (bpar, bparnew, bpardot, fbpar)

    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: allocate/dot finished." !DEBUG-KDN-110718

! TT: Added lwg1
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       g0(:,:,iglo) = 0.5 * spec(is)%temp * spec(is)%dens &
            * conjg(gnew(:,:,iglo)) * gnew(:,:,iglo) * exp(-e(ie,is))
    end do

    allocate (tot (-ntgrid:ntgrid, nakx, naky, nspec)) ; call alloc8 (c4=tot, v="tot")
    call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          do it = 1, nakx
             do is = 1, nspec             
                hk(it,ik,lwg1) = hk(it,ik,lwg1) &
                     + sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)

             end do
          end do
       end do
    end if
! <TT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Next two calls make the variables g, gnew = h, hnew 
!!! until the end of this procedure!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call g_adjust (g,    phi,    bpar,    fphi, fbpar)
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar, is_g_gnew)
    if (store_eq) call g_adjust (g_eq, phi_eq, bpar_eq, fphi, fbpar)

    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: g_adjust finished." !DEBUG-KDN-110718

! TT: variable step AB3
!    dtinv = 1./dtime
    dtinv = 1./dtime(1)
! <TT
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          
          do ig=-ntgrid, ntgrid
             
             chidot = aj0(iglo)*(phidot(ig,it,ik) - vpac(isgn,iglo) * spec(is)%stm * apardot(ig,it,ik)) &
                  + aj1vp2(iglo)*bpardot(ig,it,ik)*spec(is)%tz

!GGH             phidot1 = fphi* ((1.-fexp(is))*phinew(ig,it,ik)  - fexp(is)*phi(ig,it,ik) )*2.*dtinv
!GGH             apardot1= fapar*((1.-fexp(is))*aparnew(ig,it,ik) - fexp(is)*apar(ig,it,ik))*2.*dtinv
!GGH             bpardot1= fbpar*((1.-fexp(is))*bparnew(ig,it,ik) - fexp(is)*bpar(ig,it,ik))*2.*dtinv         
!GGH             chidot = aj0(iglo)* (  phidot1 - vpac(isgn,iglo) * spec(is)%stm * apardot1) &
!GGH                  + aj1vp2(iglo)*bpardot1*spec(is)%tz
             
             hdot = (gnew(ig,isgn,iglo) - g(ig,isgn,iglo))*dtinv
             havg =  0.5* ( gnew(ig,isgn,iglo) + g(ig,isgn,iglo) )

!GGH             hdot = ((1.-fexp(is))*gnew(ig,isgn,iglo) - fexp(is)*g(ig,isgn,iglo))*2.*dtinv
!GGH             havg =  (1.-fexp(is))*gnew(ig,isgn,iglo) + fexp(is)*g(ig,isgn,iglo)

! First term on RHS and LHS of Eq B-10 of H1:

             g0(ig,isgn,iglo) = spec(is)%dens*conjg(havg) * (chidot * spec(is)%z - hdot * spec(is)%temp)

          end do
       end do
    end do

    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: 3-layer do-loop finished." !DEBUG-KDN-110718
    
    call dealloc8 (c3=phidot,  v="phidot"  )
    call dealloc8 (c3=apardot, v="apardot")
    call dealloc8 (c3=bpardot, v="bpardot")

! TT: moved up for lwg1
!    allocate (tot (-ntgrid:ntgrid, nakx, naky, nspec)) ; call alloc8 (c4=tot, v="tot")
! <TT

    call integrate_moment (g0, tot)

    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: integrate_moment finished." !DEBUG-KDN-110718

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, nakx
                hk(it,ik,lheating(is)) = &
                     & sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
             end do
          end do
       end do
      if(debug) write(*,*) "dist_fn:get_heat: hk totalizing done." !DEBUG-KDN-110718
    end if

!-----------------------------------------------------------------------------
! Antenna Power and B-field contribution to E and E_dot---------------------
!-----------------------------------------------------------------------------
    if (proc0) then
       if(debug) write(*,*) "dist_fn:get_heat: beginning antenna calls." !DEBUG-KDN-110718
       
       allocate (j_extx(-ntgrid:ntgrid, nakx, naky)); call alloc8(c3=j_extx,v="j_extx")
       allocate (j_exty(-ntgrid:ntgrid, nakx, naky)); call alloc8(c3=j_exty,v="j_exty")
       allocate (j_extz(-ntgrid:ntgrid, nakx, naky)); call alloc8(c3=j_extz,v="j_extz")
       j_extx = 0.; j_exty = 0.; j_extz = 0.
       call antenna_apar (j_extx, j_exty, j_extz)       

       !NOTE: We'll use time decentering from species 1 for the fields 
       !WARNING: THIS MAY BE INCORRECT! GGH
       is=1

       if (beta > epsilon(0.)) then
! TT: variable step AB3
!          dtinv = 1./dtime
          dtinv = 1./dtime(1)
! <TT
          do ik=1,naky
             do it = 1,nakx
                do ig=-ntgrid, ntgrid-1
                   
                   ! J_ext.E when driving antenna only includes A_parallel:
                   apar_m = (aparnew(ig,it,ik)-apar(ig,it,ik)) * dtinv * fapar
                   allocate (epar(-ntgrid:ntgrid, nakx, naky)); call alloc8(c3=epar,v="epar")
                   call get_epar (phi, apar, phinew, aparnew, epar)
                   !This is incorrect. It only includes the d A_z / dt contribution to E_z, but k_z phi is the same order.
                   hk(it,ik,lantenna) = hk(it, ik,lantenna) - real(conjg(j_extz(ig,it,ik))*epar(ig,it,ik)) * wgt(ig)
                   !hk(it,ik,lantenna) = hk(it, ik,lantenna) + real(conjg(j_extz(ig,it,ik))*apar_m) * wgt(ig)
                   !phi_m = (phinew(ig+1,it,ik)-phinew(ig,it,ik))/delthet*gradpar*fphi
                   !hk(it,ik,lantenna) = hk(it, ik,lantenna) + real(conjg(j_extz(ig,it,ik))*phi_m) * wgt(ig) !Other component
                   call dealloc8(c3=epar,v="epar")

                   !Account for energy input of slow mode (Bpar) antennae
                   bpardot1  = (bparnew(ig,it,ik) -  bpar(ig,it,ik)) * dtinv * fbpar !Recycled variable name
                   hk(it,ik,lantenna) = hk(it, ik,lantenna) - real(conjg(j_extx(ig,it,ik))*2. * zi * aky(ik) * bpardot1 / kperp2(it,ik))*wgt(ig)
                   hk(it,ik,lantenna) = hk(it, ik,lantenna) - real(conjg(j_exty(ig,it,ik))*2. * zi * akx(it) * bpardot1 / kperp2(it,ik))*wgt(ig)
!GGH                   apar_m = ((1.-fexp(is))*aparnew(ig,it,ik)-fexp(is)*apar(ig,it,ik)) * 2.*dtinv * fapar

                   bpar2  = cabs(bpar(ig,it,ik))**2 * fbpar
                   bperp2 = cabs(apar(ig,it,ik))**2 * fapar * kperp2(it,ik) * 0.25 ! 0.25 from normalizations

                   bpar2new  = cabs(bparnew(ig,it,ik))**2 * fbpar
                   bperp2new = cabs(aparnew(ig,it,ik))**2 * fapar * kperp2(it,ik) * 0.25 ! 0.25 from normalizations

                   fac3 = wgt(ig) * (2.0/beta)

                   hk(it,ik,lenergy_dot) = hk(it,ik,lenergy_dot) + 0.5 * (bpar2new + bperp2new - (bpar2 + bperp2)) * dtinv * fac3

!GGH                   hk(it,ik,lenergy_dot = hk(it,ik,lenergy_dot + 0.5 * ((1.-fexp(is))*(bpar2new + bperp2new) - fexp(is)*(bpar2 + bperp2)) * dtinv * fac3

                   hk(it,ik,lenergy) = hk(it,ik,lenergy) + 0.5 * (bpar2new + bperp2new) * fac3

                   !Eapar = int k_perp^2 A_par^2/(8 pi)                   
                   !RN> Eapar = int k_perp^2 A_par^2/(4beta0) [normalized unit]
                   hk(it,ik,leapar) = hk(it,ik,leapar) + 0.5 * bperp2new * fac3 

                   !Ebpar = int B_par^2/(8 pi)
                   !RN> Ebpar = int B_par^2/beta0 [normalize unit]
                   hk(it,ik,lebpar) = hk(it,ik,lebpar) + 0.5 * bpar2new * fac3 
                   
                end do
             end do
          end do
          if (store_eq) then
             do ik=1,naky
                do it = 1,nakx
                   do ig=-ntgrid, ntgrid
                      bperp_pt2 = cabs(aparnew(ig,it,ik)-apar_eq(ig,it,ik))**2 &
                           & * fapar * kperp2(it,ik) * 0.25
                      hk(it,ik,leapar_pt) = hk(it,ik,leapar_pt) &
                           & + 0.5 * bperp_pt2 * fac3
                      bpar_pt2 = cabs(bparnew(ig,it,ik)-bpar_eq(ig,it,ik))**2 &
                           & * fbpar
                      hk(it,ik,lebpar_pt) = hk(it,ik,lebpar_pt) &
                           & + 0.5 * bpar_pt2 * fac3 
                   end do
                end do
             end do
          end if
       else ! if not beta > 0
          hk(:,:,lantenna) = 0.
          hk(:,:,lenergy_dot) = 0.
          hk(:,:,lenergy) = 0.
          hk(:,:,lebpar) = 0.

          ! put the 2D ES invariant: chosen ig=0 to output
          hk(1:nakx,1:naky,leapar) = gamtot(1:nakx,1:naky) * &
               cabs(phinew(0,1:nakx,1:naky))**2 / 2.0
          ! TT: yavg electron response is taken care of by gamtot
       end if
!       deallocate (j_ext)
       call dealloc8(c3=j_extx,v="j_extx")
       call dealloc8(c3=j_exty,v="j_exty")
       call dealloc8(c3=j_extz,v="j_extz")
      if(debug) write(*,*) "dist_fn:get_heat: antenna calls finished." !DEBUG-KDN-110718

    end if

!-----------------------------------------------------------------------------
! Finish E_dot--------------------------------------------------------------
!-----------------------------------------------------------------------------

!GGH Include response of Boltzmann species for single-species runs

! TT> need to check if this is correct (10/13/08)
!    if (.not. has_electron_species(spec)) then
    if ( (.not. has_electron_species(spec)) .and. &
         (adiabatic_option_switch /= adiabatic_option_zero) ) then
! <TT
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
! TT: variable step AB3
!          dtinv = 1./dtime
          dtinv = 1./dtime(1)
! <TT
          do ik=1,naky
! TT: correction of W invariant for yavg electron response
             if (ik==1 .and. adiabatic_option_switch == adiabatic_option_yavg) &
                  cycle
             do it = 1,nakx
                do ig=-ntgrid, ntgrid-1

                   !WARNING: THIS IS NOT UPDATED FOR TIME DECENTERING!!!!! GGH
                   phi_m   = (phinew(ig,it,ik) - phi(ig,it,ik))*dtinv

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik,lenergy_dot) = hk(it,ik,lenergy_dot) + &
                        fphi * real(conjg(phinew(ig,it,ik))*phi_m) &
                        * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)

                end do
             end do
          end do
       endif
    endif !END Correction to E_dot for single species runs---------------------
 
!GGH New E_dot calc
! TT: variable step AB3
!    dtinv = 1./dtime
    dtinv = 1./dtime(1)
! <TT
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2

          do ig=-ntgrid, ntgrid
             
             !Calculate old fluctuating energy de
             havg = g(ig,isgn,iglo)

             j0phiavg = phi(ig,it,ik) * aj0(iglo) * fphi * spec(is)%zt

             phi_avg = phi(ig,it,ik) * fphi * spec(is)%zt

             de = 0.5 * spec(is)%temp * spec(is)%dens * &
                  ( conjg(havg) * havg &
                  + conjg(phi_avg) * phi_avg &
                  - conjg(j0phiavg) * havg &
                  - conjg(havg) * j0phiavg )

            !Calculate new fluctuating energy denew
             havg = gnew(ig,isgn,iglo)

             j0phiavg = phinew(ig,it,ik) * aj0(iglo) * fphi * spec(is)%zt

             phi_avg = phinew(ig,it,ik) * fphi * spec(is)%zt

             denew=0.5*spec(is)%temp*spec(is)%dens * &
                  ( conjg(havg) * havg &
                  + conjg(phi_avg) * phi_avg &
                  - conjg(j0phiavg) * havg &
                  - conjg(havg) * j0phiavg) 

             !Set g0 as the change of energy (denew-de)/dt
             g0(ig,isgn,iglo) = (denew - de) *dtinv

!GGH             g0(ig,isgn,iglo) = ((1.-fexp(is))*denew-fexp(is)*de) * 2.*dtinv

          end do
       end do
    end do
    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: E dot computation finished?" !DEBUG-KDN-110718
    !GGH -END new e_dot calc

    call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          do it = 1, nakx
             do is = 1, nspec
                hk(it,ik,lenergy_dot) = hk(it,ik,lenergy_dot) + &
                     & sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
             end do
          end do
       end do
      if(debug) write(*,*) "dist_fn:get_heat: E dot summation done." !DEBUG-KDN-110718
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------------
! Gradient Contributions to Heating-----------------------------------------
!-----------------------------------------------------------------------------

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)

       do isgn=1,2
          do ig=-ntgrid, ntgrid

             chi = phinew(ig,it,ik) * aj0(iglo) * fphi &
                  - aparnew(ig,it,ik) * aj0(iglo) * vpac(isgn,iglo) * spec(is)%stm * fapar & 
                  + bparnew(ig,it,ik) * aj1vp2(iglo) * spec(is)%tz * fbpar

             havg = gnew(ig,isgn,iglo)

             g0(ig,isgn,iglo) = zi * wstar(ik,ie,is) * dtinv * conjg(havg) * chi * spec(is)%dens
            
          end do
       end do
    end do

    call integrate_moment (g0, tot)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, nakx
                hk(it,ik,lgradients(is)) = &
                     & sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
             end do
          end do
       end do
      if(debug) write(*,*) "dist_fn:get_heat: gradient contribution done." !DEBUG-KDN-110718
    end if

! GGH 02APR08- Commented out hyperviscosity and hyperresistivity
!-----------------------------------------------------------------------------
! Hyperviscosity------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!    if (D_v > epsilon(0.)) then
!
!       do iglo=g_lo%llim_proc, g_lo%ulim_proc
!          is = is_idx(g_lo, iglo)
!          it = it_idx(g_lo, iglo)
!          ik = ik_idx(g_lo, iglo)
!          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
!           do isgn=1,2
!             do ig=-ntgrid, ntgrid
!                
!                havg = gnew (ig,isgn,iglo)
!
!                j0phiavg = phinew(ig,it,ik) * aj0(iglo) * fphi * spec(is)%zt
!
!                j1bparavg= bparnew(ig,it,ik) * aj1vp2(iglo) * fbpar
!
!!Set g0 for hyperviscous heating
!                g0(ig,isgn,iglo) = spec(is)%dens * spec(is)%temp * D_v * akperp4 * &
!                     ( conjg(havg) * havg - conjg(havg) * j0phiavg - conjg(havg) * j1bparavg)
!
!             end do
!          end do
!       end do
!
!       call integrate_moment (g0, tot)
!       if (proc0) then
!          do ik = 1, naky
!             do it = 1, nakx
!                do is = 1, nspec
!                   hk(it,ik,lhypervisc(is)) = sum(real(tot(-ntgrid:ntgrid,it,ik,is)) * wgt(-ntgrid:ntgrid))
!                end do
!             end do
!          end do
!       end if
!
!    end if !End Hyperviscous Heating Calculation
!
!
!-----------------------------------------------------------------------------
! Hyperresistivity------------------------------------------------------------
!-----------------------------------------------------------------------------
! 
!    if (D_eta > epsilon(0.)) then
!
!       do iglo=g_lo%llim_proc, g_lo%ulim_proc
!          is = is_idx(g_lo, iglo)
!          it = it_idx(g_lo, iglo)
!          ik = ik_idx(g_lo, iglo)
!          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
!           do isgn=1,2
!             do ig=-ntgrid, ntgrid
!                
!                havg = gnew(ig,it,ik)
!
!                j0aparavg = aparnew(ig,it,ik) * aj0(iglo) * fapar * spec(is)%zstm * vpac(isgn,iglo)
!
!!Set g0 for hyperresistive heating
!                g0(ig,isgn,iglo) = spec(is)%dens * spec(is)%temp * D_eta * akperp4 * &
!                     conjg(havg) * j0aparavg
!
!             end do
!          end do
!       end do
!
!       call integrate_moment (g0, tot)
!       if (proc0) then
!          do ik = 1, naky
!             do it = 1, nakx
!                do is = 1, nspec
!                   hk(it,ik,lhyperres(is)) = sum(real(tot(-ntgrid:ntgrid,it,ik,is))*wgt(-ntgrid:ntgrid))
!                end do
!             end do
!          end do
!       end if
!
!    end if !End Hyperresistivity Heating Calculation
!
!-----------------------------------------------------------------------------
!Finish Energy-------------------------------------------------------------
!-----------------------------------------------------------------------------

!GGH Calculate hs2-------------------------------------------------------------
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid, ntgrid
             g0(ig,isgn,iglo) = 0.5 * spec(is)%temp * spec(is)%dens * &
                  & cabs(gnew(ig,isgn,iglo))**2
          end do
       end do
    end do
    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: 3-fold do loop in hs2 done." !DEBUG-KDN-110718

    call integrate_moment (g0, tot)
    if(debug.and.proc0) write(*,*) "dist_fn:get_heat: integrate_moment for hs2 done." !DEBUG-KDN-110718

    if (proc0) then
       do ik = 1, naky
          do it = 1, nakx
             do is = 1, nspec             
                !hs2 = int_r int_v T/F0 hs^2/2
                hk(it,ik,lhs2(is)) = sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
             end do
          end do
       end do
      if(debug) write(*,*) "dist_fn:get_heat: hs2 summation done." !DEBUG-KDN-110718
    end if

    if (store_eq) then
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          do isgn=1,2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = 0.5 * spec(is)%temp * spec(is)%dens * &
                     & cabs(gnew(ig,isgn,iglo)-g_eq(ig,isgn,iglo))**2
             end do
          end do
       end do

       call integrate_moment (g0, tot)

       if (proc0) then
          do ik = 1, naky
             do it = 1, nakx
                do is = 1, nspec             
                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik,lhs2_pt(is)) = sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
                end do
             end do
          end do
       end if
    end if

!Calculate phis2-------------------------------------------------------------
    if (proc0) then
       do ik=1,naky
          do it = 1,nakx
             do ig=-ntgrid, ntgrid-1
                do is = 1, nspec
                   phi_avg = phinew(ig,it,ik) * fphi * spec(is)%zt

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik,lphis2(is)) = hk(it,ik,lphis2(is))  &
                        + 0.5 * spec(is)%temp * spec(is)%dens * cabs(phi_avg)**2 &
                        * wgt(ig)
                enddo
             end do
          end do
       end do
       if (store_eq) then
          do ik=1,naky
             do it = 1,nakx
                do ig=-ntgrid, ntgrid-1
                   do is = 1, nspec
                      phi_pt_avg = (phinew(ig,it,ik)-phi_eq(ig,it,ik)) &
                           & * fphi * spec(is)%zt
                      hk(it,ik,lphis2_pt(is)) = hk(it,ik,lphis2_pt(is))  &
                           & + 0.5 * spec(is)%temp * spec(is)%dens &
                           & * cabs(phi_pt_avg)**2 * wgt(ig)
                   end do
                end do
             end do
          end do
       endif
    endif

! Calculate u2 (kinetic energy)-----------------------------------------------
    allocate(upar(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=upar,v="upar")
    allocate(ux(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=ux,v="ux")
    allocate(uy(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=uy,v="uy")
    allocate(dens(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=dens,v="dens")
    allocate(pxx(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=pxx,v="pxx")
    allocate(pyy(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=pyy,v="pyy")
    allocate(pzz(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8 (c4=pzz,v="pzz")
    upar=0.; ux=0.; uy=0.
    dens=0.; pxx=0.; pyy=0.; pzz=0.
    call getmoms_r(dens=dens,uz=upar,ux=ux,uy=uy,pxx=pxx,pyy=pyy,pzz=pzz,adjust_in=.false.)
    if (proc0) then
       do ik=1,naky
          do it = 1,nakx
             do ig=-ntgrid, ntgrid-1
                do is = 1, nspec
                   akperp2=aky(ik)**2+akx(it)**2
                   hk(it,ik,lkphi2(is)) = hk(it,ik,lkphi2(is))  &
                        & + wgt(ig) * spec(is)%dens * spec(is)%mass * .25 &
                        & * akperp2*cabs(phinew(ig,it,ik))**2*fphi
                   hk(it,ik,luperp2(is)) = hk(it,ik,luperp2(is))  &
                        & + wgt(ig) * spec(is)%dens * spec(is)%mass &
                        & * (cabs(ux(ig,it,ik,is))**2+cabs(uy(ig,it,ik,is))**2)
                   hk(it,ik,lupar2(is)) = hk(it,ik,lupar2(is))  &
                        & + wgt(ig) * spec(is)%dens * spec(is)%mass &
                        & * cabs(upar(ig,it,ik,is))**2
                   hk(it,ik,ltxx(is)) = hk(it,ik,ltxx(is))  &
                        & + wgt(ig) * spec(is)%temp * &
                        & ( conjg(dens(ig,it,ik,is)) * (pxx(ig,it,ik,is) - &
                        & spec(is)%temp*dens(ig,it,ik,is) ) + &
                        & dens(ig,it,ik,is) * conjg(pxx(ig,it,ik,is) - &
                        & spec(is)%temp*dens(ig,it,ik,is) ) )
                   hk(it,ik,ltyy(is)) = hk(it,ik,ltyy(is))  &
                        & + wgt(ig) * spec(is)%temp * &
                        & ( conjg(dens(ig,it,ik,is)) * (pyy(ig,it,ik,is) - &
                        & spec(is)%temp*dens(ig,it,ik,is) ) + &
                        & dens(ig,it,ik,is) * conjg(pyy(ig,it,ik,is) - &
                        & spec(is)%temp*dens(ig,it,ik,is) ) )
                   hk(it,ik,ltzz(is)) = hk(it,ik,ltzz(is))  &
                        & + wgt(ig) * spec(is)%temp * &
                        & ( conjg(dens(ig,it,ik,is)) * (pzz(ig,it,ik,is) - &
                        & spec(is)%temp*dens(ig,it,ik,is) ) + &
                        & dens(ig,it,ik,is) * conjg(pzz(ig,it,ik,is) - &
                        & spec(is)%temp*dens(ig,it,ik,is) ) )
                end do
             end do
          end do
       end do
    end if

    if (store_eq) then
       gnew=gnew-g_eq
       call getmoms_r(uz=upar,ux=ux,uy=uy,adjust_in=.false.)
       gnew=gnew+g_eq
       if (proc0) then
          do ik=1,naky
             do it = 1,nakx
                do ig=-ntgrid, ntgrid-1
                   do is = 1, nspec
                      akperp2=aky(ik)**2+akx(it)**2
                      hk(it,ik,lkphi2_pt(is)) = hk(it,ik,lkphi2_pt(is))  &
                           & + wgt(ig) * spec(is)%dens * spec(is)%mass * .25 &
                           & * akperp2*cabs(phinew(ig,it,ik)-phi_eq(ig,it,ik))**2*fphi
                      hk(it,ik,luperp2_pt(is)) = hk(it,ik,luperp2_pt(is))  &
                           & + wgt(ig) * spec(is)%dens * spec(is)%mass &
                           & * (cabs(ux(ig,it,ik,is))**2+cabs(uy(ig,it,ik,is))**2)
                      hk(it,ik,lupar2_pt(is)) = hk(it,ik,lupar2_pt(is))  &
                           & + wgt(ig) * spec(is)%dens * spec(is)%mass &
                           & * cabs(upar(ig,it,ik,is))**2

                      hk(it,ik,ltxx_pt(is)) = hk(it,ik,ltxx_pt(is))  &
                           & + wgt(ig) * spec(is)%temp * &
                           & ( conjg(dens(ig,it,ik,is)) * (pxx(ig,it,ik,is) - &
                           & spec(is)%temp*dens(ig,it,ik,is) ) + &
                           & dens(ig,it,ik,is) * conjg(pxx(ig,it,ik,is) - &
                           & spec(is)%temp*dens(ig,it,ik,is) ) )
                      hk(it,ik,ltyy_pt(is)) = hk(it,ik,ltyy_pt(is))  &
                           & + wgt(ig) * spec(is)%temp * &
                           & ( conjg(dens(ig,it,ik,is)) * (pyy(ig,it,ik,is) - &
                           & spec(is)%temp*dens(ig,it,ik,is) ) + &
                           & dens(ig,it,ik,is) * conjg(pyy(ig,it,ik,is) - &
                           & spec(is)%temp*dens(ig,it,ik,is) ) )
                      hk(it,ik,ltzz_pt(is)) = hk(it,ik,ltzz_pt(is))  &
                           & + wgt(ig) * spec(is)%temp * &
                           & ( conjg(dens(ig,it,ik,is)) * (pzz(ig,it,ik,is) - &
                           & spec(is)%temp*dens(ig,it,ik,is) ) + &
                           & dens(ig,it,ik,is) * conjg(pzz(ig,it,ik,is) - &
                           & spec(is)%temp*dens(ig,it,ik,is) ) )
                   end do
                end do
             end do
          end do
       end if
    end if
    call dealloc8 (c4=upar, v="upar")
    call dealloc8 (c4=ux, v="ux")
    call dealloc8 (c4=uy, v="uy")
    call dealloc8 (c4=dens, v="dens")
    call dealloc8 (c4=pxx, v="pxx")
    call dealloc8 (c4=pyy, v="pyy")
    call dealloc8 (c4=pzz, v="pzz")

! Calculate delfs2 (rest of energy)-----------------------------------------------

!GGH  Include response of Boltzmann species for single species runs
    if ( (.not. has_electron_species(spec)) .and. &
         (adiabatic_option_switch /= adiabatic_option_zero) ) then
! TT: variable step AB3
!       dtinv = 1./dtime
       dtinv = 1./dtime(1)
! <TT
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             ! zonal component correction for yavg option
             if ( adiabatic_option_switch == adiabatic_option_yavg &
                  .and. aky(ik) < epsilon(0.0) ) cycle
             do it = 1,nakx
                do ig=-ntgrid, ntgrid-1

                   phi_avg = phinew(ig,it,ik) * fphi

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik,lenergy) = hk(it,ik,lenergy) + &
                        fphi * cabs(phi_avg)**2 &
                        * 0.5 * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)

                end do
             end do
          end do
       endif
    endif !END Correction to energy for single species runs---------------------

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2

          do ig=-ntgrid, ntgrid
             
             havg = gnew(ig,isgn,iglo)

             j0phiavg = phinew(ig,it,ik) * aj0(iglo) * fphi * spec(is)%zt

             phi_avg =  phinew(ig,it,ik) * fphi * spec(is)%zt

             g0(ig,isgn,iglo) = 0.5 * spec(is)%temp * spec(is)%dens &
                  * (conjg(havg) * havg &
                  + conjg(phi_avg) * phi_avg &
                  - conjg(j0phiavg) * havg &
                  - conjg(havg) * j0phiavg)
          end do
       end do
    end do

    call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          do it = 1, nakx
             do is = 1, nspec             
                hk(it,ik,lenergy)=hk(it,ik,lenergy)+sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)

                !Delfs2 = int_r int_v T/F0 dfs^2/2
                hk(it,ik,ldelfs2(is))=sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
             end do
          end do
       end do
    end if

    if (store_eq) then
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          do isgn=1,2

             do ig=-ntgrid, ntgrid
             
                havg = gnew(ig,isgn,iglo)-g_eq(ig,isgn,iglo)
                j0phiavg = (phinew(ig,it,ik)-phi_eq(ig,it,ik)) &
                     & * aj0(iglo) * fphi * spec(is)%zt
                phi_avg =  (phinew(ig,it,ik)-phi_eq(ig,it,ik)) &
                     & * fphi * spec(is)%zt
                g0(ig,isgn,iglo) = 0.5 * spec(is)%temp * spec(is)%dens &
                     * (conjg(havg) * havg &
                     + conjg(phi_avg) * phi_avg &
                     - conjg(j0phiavg) * havg &
                     - conjg(havg) * j0phiavg)
             end do
          end do
       end do

       call integrate_moment (g0, tot)

       if (proc0) then
          do ik = 1, naky
             do it = 1, nakx
                do is = 1, nspec             
                   hk(it,ik,ldelfs2_pt(is))= &
                        & sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
                end do
             end do
          end do
       end if
    end if

!!! n0 * T0 * (q phi/T0)^2 (1-gamma0) 
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       do isgn=1,2
          do ig=-ntgrid, ntgrid
             g0(ig,isgn,iglo) = 0.5 * spec(is)%temp * spec(is)%dens &
                  & * cabs(spec(is)%zt * phinew(ig,it,ik) * fphi)**2 &
                  & * (1.-aj0(iglo)**2)
          end do
       end do
    end do

    call integrate_moment (g0, tot)

    if (proc0) then
       do ik = 1, naky
          do it = 1, nakx
             do is = 1, nspec             
                hk(it,ik,lg0phi2(is))= &
                     & sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
             end do
          end do
       end do
    end if

    if (store_eq) then
       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          do isgn=1,2
             do ig=-ntgrid, ntgrid
                g0(ig,isgn,iglo) = 0.5 * spec(is)%temp * spec(is)%dens &
                     & * cabs(spec(is)%zt * (phinew(ig,it,ik)-phi_eq(ig,it,ik)) * fphi)**2 &
                     & * (1.-aj0(iglo)**2)
             end do
          end do
       end do

       call integrate_moment (g0, tot)

       if (proc0) then
          do ik = 1, naky
             do it = 1, nakx
                do is = 1, nspec             
                   hk(it,ik,lg0phi2_pt(is))= &
                        & sum(real(tot(-ntgrid:ntgrid-1,it,ik,is)) * wgt)
                end do
             end do
          end do
       end if
    end if

    call dealloc8 (c4=tot, v="tot")

!!
!! Put g, gnew back to their usual meanings
!!
    call g_adjust (g,    phi,    bpar,    -fphi, -fbpar)
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar, is_g_gnew)
    if (store_eq) call g_adjust (g_eq, phi_eq, bpar_eq, -fphi, -fbpar)

    ! RN> moved from subroutine heating in agk_diagnostics
    if (proc0) then
       !Calculate collisional heating
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, nakx
                !Sum heating by k over all z points (ig)
                hk(it, ik, lcollisions(is)) = spec(is)%temp*spec(is)%dens* &
                     sum(real(c_rate(-ntgrid:ntgrid-1,it,ik,is,1))*wgt)

                hk(it, ik, lhypercoll(is)) = spec(is)%temp*spec(is)%dens* &
                     sum(real(c_rate(-ntgrid:ntgrid-1,it,ik,is,2))*wgt)
                
                hk(it, ik, limp_colls(is)) = spec(is)%temp*spec(is)%dens* &
                     sum(real(c_rate(-ntgrid:ntgrid-1,it,ik,is,3))*wgt)
             end do
          end do
       end do
    end if

    ! all branching for it=ik=1, and ik=1 are done here
    if(nonlin) hk(1,1,:)=0.
    hk(1:nakx,2:naky,:)=.5*hk(1:nakx,2:naky,:)
    
!    if(proc0) deallocate (wgt)
    if(proc0) call dealloc8 (r1=wgt,v="wgt")

  end subroutine get_heat
!-----------------------------------------------------------------------------

  subroutine reset_init

    use dist_fn_arrays, only: gnew, g
    initializing  = .true.
    initialized = .false.
    
    wdrift = 0.
    a = 0. ; b = 0. ; r = 0. ; ainv = 0.
    g = 0. ; gnew = 0.  ;  g0 = 0.

  end subroutine reset_init

  subroutine reset_physics

    call init_wstar

  end subroutine reset_physics

!>MAB
  subroutine get_verr (errest, erridx, phi, bpar)

    use mp, only: proc0
    use le_grids, only: integrate_test, integrate_species
    use le_grids, only: eint_error, lint_error
    use le_grids, only: nlambda, negrid, ecut
    use le_grids, only: vgrid, nesub
    use egrid, only: x0
    use theta_grid, only: ntgrid
    use kgrids, only: nakx, naky, aky, akx
    use species, only: nspec, spec
    use dist_fn_arrays, only: gnew, aj0, vpa
    use run_parameters, only: fphi, fapar, fbpar, beta
    use agk_layouts, only: g_lo
    use nonlinear_terms, only: nonlin
    use agk_mem, only: alloc8, dealloc8
    implicit none

    integer :: ig, it, ik, is, ipt, iglo, isgn, np

    integer, dimension (:,:), intent (out) :: erridx
    real, dimension (:,:), intent (out) :: errest
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar

    real, dimension (:), allocatable :: wgt
    real, dimension (:), allocatable :: errtmp
    integer, dimension (:), allocatable :: idxtmp
    real, dimension (:,:), allocatable, save :: kmax
    complex, dimension (:,:,:), allocatable :: phi_app, apar_app
    complex, dimension (:,:,:,:), allocatable :: phi_e, phi_l
    complex, dimension (:,:,:,:), allocatable :: apar_e, apar_l

    real :: errcut_phi, errcut_apar

    allocate (wgt(nspec)); call alloc8(r1=wgt,v="wgt")
    allocate (errtmp(2)); call alloc8(r1=errtmp,v="errtmp")
    allocate (idxtmp(3)); call alloc8(i1=idxtmp,v="idetmp")

    if (vgrid) then
       np = nesub
    else
       if (ecut > 0.0) then
          np = negrid - 1
       else
          np = negrid
       end if
    end if

    if (fphi > epsilon(0.0)) then
       allocate(phi_app(-ntgrid:ntgrid,nakx,naky)) ; call alloc8 (c3=phi_app, v="phi_app")
       allocate(phi_e(-ntgrid:ntgrid,nakx,naky,np)) ; call alloc8 (c4=phi_e, v="phi_e")  ! TT
       allocate(phi_l(-ntgrid:ntgrid,nakx,naky,nlambda)) ; call alloc8 (c4=phi_l, v="phi_l")
    end if

    if (fapar > epsilon(0.0)) then
       allocate(apar_app(-ntgrid:ntgrid,nakx,naky)) ; call alloc8 (c3=apar_app, v="apar_app")
       allocate(apar_e(-ntgrid:ntgrid,nakx,naky,np)) ; call alloc8 (c4=apar_e, v="apar_e")  ! TT
       allocate(apar_l(-ntgrid:ntgrid,nakx,naky,nlambda)) ; call alloc8 (c4=apar_l, v="apar_l") 
    end if

! first call to g_adjust converts gyro-averaged dist. fn. (g)
! into nonadiabatic part of dist. fn. (h)

    call g_adjust (gnew, phi, bpar, fphi, fbpar, is_g_gnew)

! take gyro-average of h at fixed total position (not g.c. position)
    if (fphi > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(iglo)*gnew(:,isgn,iglo)
          end do
       end do
       wgt = spec%z * spec%dens
       call integrate_species (g0, wgt, phi_app)
       ! now phi_app = qn int J0*h dv

!    call integrate_test (g0, wgt, phi_app, istep)  ! only around for testing

! integrates dist fn of each species over v-space
! after dropping an energy grid point and returns
! phi_e, which contains the integral approximations
! to phi for each point dropped

       call eint_error (g0, wgt, phi_e)

! integrates dist fn of each species over v-space
! after dropping a lambda grid point and returns phi_l.
! phi_l contains nlambda approximations for the integral over lambda that
! come from dropping different pts from the gaussian quadrature grid

       call lint_error (g0, wgt, phi_l)

    end if

    if (fapar > epsilon(0.0)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g0(:,isgn,iglo) = aj0(iglo)*vpa(isgn,iglo)*gnew(:,isgn,iglo)
          end do
       end do
       
       wgt = 2.0 * beta * spec%z * spec%dens * sqrt(spec%temp/spec%mass)
       call integrate_species (g0, wgt, apar_app)
       ! apar_app = 2 beta qn T/m int vpar*J0*h dv

!    call integrate_test (g0, wgt, apar_app, istep)  ! only around for testing

       call eint_error (g0, wgt, apar_e)
       call lint_error (g0, wgt, apar_l)

    end if

! second call to g_adjust converts from h back to g

    call g_adjust (gnew, phi, bpar, -fphi, -fbpar, is_g_gnew)

    if (.not. allocated(kmax)) then
       allocate (kmax(nakx, naky)); call alloc8(r2=kmax,v="kmax")
       do ik = 1, naky
          do it = 1, nakx
             kmax(it,ik) = max(abs(akx(it)),abs(aky(ik)))
          end do
       end do
    end if
    
    errest = 0.0
    erridx = 0
    
    if (fphi > epsilon(0.0)) then
!!$       errcut_phi = 0.0
!!$       do ig = -ntgrid, ntgrid
!!$          errcut_phi = max(errcut_phi, maxval(cabs(phi_app(ig,:,:))*kmax))
!!$       end do
!!$       errcut_phi = errcut_phi/100
       ! TT: do the same thing
       errcut_phi = maxval( maxval(abs(phi_app),1) * kmax ) / 100.

       call estimate_error (phi_app, phi_e, kmax, errcut_phi, errtmp, idxtmp)
       errest(1,:) = errtmp
       erridx(1,:) = idxtmp
       
       call estimate_error (phi_app, phi_l, kmax, errcut_phi, errtmp, idxtmp)
       errest(2,:) = errtmp
       erridx(2,:) = idxtmp

       call dealloc8 (c3=phi_app, v="phi_app")
       call dealloc8 (c4=phi_e, v="phi_e")
       call dealloc8 (c4=phi_l, v="phi_l")
    end if
    
    if (fapar > epsilon(0.0)) then
!!$       errcut_apar = 0.0
!!$       do ig = -ntgrid, ntgrid
!!$          errcut_apar = max(errcut_apar, maxval(cabs(apar_app(ig,:,:))*kmax))
!!$       end do
!!$       errcut_apar = errcut_apar/100
       ! TT: do the same thing
       errcut_apar = maxval( maxval(cabs(apar_app),1) * kmax ) / 100.

       call estimate_error (apar_app, apar_e, kmax, errcut_apar, errtmp, idxtmp)
       errest(3,:) = errtmp
       erridx(3,:) = idxtmp
       
       call estimate_error (apar_app, apar_l, kmax, errcut_apar, errtmp, idxtmp)
       errest(4,:) = errtmp
       erridx(4,:) = idxtmp

       call dealloc8 (c3=apar_app, v="apar_app")
       call dealloc8 (c4=apar_e, v="apar_e")
       call dealloc8 (c4=apar_l, v="apar_l")
    end if

!    deallocate (wgt, errtmp, idxtmp)
    call dealloc8(r1=wgt,v="wgt")
    call dealloc8(r1=errtmp,v="errtmp")
    call dealloc8(i1=idxtmp,v="idxtmp")

  end subroutine get_verr

  subroutine estimate_error (app1, app2, kmax, errcut, errest, erridx)

    use kgrids, only: naky, nakx
    use theta_grid, only: ntgrid

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: app1
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: app2
    real, dimension (:,:), intent (in) :: kmax
    real, intent (in) :: errcut
    real, dimension (:), intent (out) :: errest
    integer, dimension (:), intent (out) :: erridx

    integer :: ik, it, ig, ipt
    integer :: igmax, ikmax, itmax, gpcnt
    real :: gdsum, gdmax, gpavg, gnsum, gsmax, gpsum, gptmp
! TT>
    real, parameter :: eps = 1.e-10
! <TT

    igmax = 0; ikmax = 0; itmax = 0
    gdsum = 0.0; gdmax = 0.0; gpavg = 0.0; gnsum = 0.0; gsmax = 0.0

    do ik = 1, naky
       do it = 1, nakx
          do ig=-ntgrid,ntgrid
             gpcnt = 0; gpsum = 0.0
             ! errcut = max (kmax*{phi,apar}_app) / 100
             if ( kmax(it,ik)*abs(app1(ig,it,ik)) &
                  < max(errcut, 10*epsilon(0.0)) ) cycle
             do ipt=1, size(app2(0,1,1,:))
                gptmp = kmax(it,ik)*cabs(app1(ig,it,ik) - app2(ig,it,ik,ipt))
                gpsum = gpsum + gptmp
                gpcnt = gpcnt + 1
             end do
             gpavg = gpsum/gpcnt
             ! gpavg is the average absolute error of kmax * int J0*h dv
             ! for a given (ig,it,ik)

             if (gpavg > gdmax) then
                ! gdmax is the maximum of gpavg for computed (ig,it,ik)
                igmax = ig
                ikmax = ik
                itmax = it
                gdmax = gpavg
                gsmax = kmax(it,ik)*cabs(app1(ig,it,ik))
             end if

             gnsum = gnsum + gpavg
             gdsum = gdsum + kmax(it,ik)*cabs(app1(ig,it,ik))

             if (test .and. ig==0) then
                print *, 'it, ik, gpavg, gnsum, gdsum = ', it, ik, gpavg, gnsum, gdsum
                print *, 'app1, app2 = ', abs(app1(ig,it,ik)), abs(app2(ig,it,ik,:))
             end if

          end do
       end do
    end do

    if(gsmax /= 0.) gdmax = gdmax/gsmax
    ! gsmax is kmax * int J0*h dv corresponding to max_(ig,it,ik) gpavg

    erridx(1) = igmax
    erridx(2) = ikmax
    erridx(3) = itmax
    errest(1) = gdmax
! TT: lowest limit for relative error
!    errest(2) = gnsum/gdsum
    errest(2) = gnsum
    if (gdsum > eps) errest(2) = errest(2)/gdsum
! <TT

  end subroutine estimate_error
!<MAB

  subroutine write_f (last)

    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use agk_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx, ie_idx
    use agk_layouts, only: idx_local, proc_id
    use le_grids, only: al, e, negrid, nlambda
    use theta_grid, only: ntgrid
    use agk_time, only: time
    use dist_fn_arrays, only: gnew  ! changed from g to gnew - MAB
    use run_parameters, only: fphi, fbpar   ! MAB
    use fields_arrays, only: phinew, bparnew
    use agk_mem, only: alloc8

    integer :: iglo, ik, it, is
    integer :: ie, il, ig
    integer, save :: unit
    real :: vpa, vpe
    complex, dimension(2) :: gtmp
    logical :: first = .true.
    logical, intent(in)  :: last 

! xpts and ypts are only temporary additions so that Matlab script works
! in current form (designed for gs2)
    real, dimension (:), allocatable, save :: xpts, ypts

    if (.not. allocated(xpts)) then
       allocate(xpts(negrid)); call alloc8(r1=xpts,v="xpts")
    end if
    if (.not. allocated(ypts)) then
       allocate(ypts(nlambda)); call alloc8(r1=ypts,v="ypts")
    end if

    xpts = 0.0
    ypts = 0.0

    if (proc0) then
       if (first) then 
          call get_unused_unit (unit)
          call open_output_file (unit, ".dist")
          write (unit,*) negrid*nlambda
          first = .false.
       end if
    endif

    call g_adjust (gnew, phinew, bparnew, fphi, fbpar, is_g_gnew)

    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo, iglo) ; if (ik /= 1) cycle
       it = it_idx(g_lo, iglo) ; if (it /= 1) cycle
       is = is_idx(g_lo, iglo) !; if (is /= 1) cycle
       ie = ie_idx(g_lo, iglo) 
       ig = 0
       il = il_idx(g_lo, iglo)
       if (idx_local (g_lo, ik, it, il, ie, is)) then
          if (proc0) then 
             gtmp = gnew(ig,:,iglo)  ! replaced g with gnew - MAB
          else
             call send (gnew(ig,:,iglo), 0)   ! replaced g with gnew - MAB
          endif
       else if (proc0) then
          call receive (gtmp, proc_id(g_lo, iglo))
       endif
       if (proc0) then
          vpa = sqrt(e(ie,is)*max(0.0, 1.0-al(il)))
          vpe = sqrt(e(ie,is)*al(il))
!             write (unit, "(8(1x,e13.6))") vpe,vpa,gtmp(1),-vpa,gtmp(2)             
!          gtmp = gtmp* exp(-e(ie,is)) ! can be a problem for large ecut
          write (unit, "(8(1x,e13.6),1x,i2)") vpe,vpa,gtmp(1),-vpa,gtmp(2),time,is
! changed to work with Matlab script for generating 2D plots of g(v)
!          write (unit, "(8(1x,e13.6))") vpa, vpe, e(ie,is), al(il), &
!               xpts(ie), ypts(il), real(gtmp(1)), real(gtmp(2))
       end if
    end do
    if (proc0) write (unit, *)
    if (last .and. proc0) call close_output_file (unit)
    
    call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar, is_g_gnew)

  end subroutine write_f
!------------------------------------------------------------------------------
! Velocity Plane Analysis---GGH 23JAN08
!------------------------------------------------------------------------------
! Routine to output slices of velocity space
! For each (it,ik) Fourier mode, it produces a separate file at ig=0 plane
!     for g       runname.vpx(it)y(ik)g
!     for h       runname.vpx(it)y(ik)h
  subroutine write_vp 
    use mp, only: proc0, send, receive, sum_reduce
    use le_grids, only: al, e, negrid, nlambda
    use species, only: nspec
    use kgrids, only: nakx, naky
    use nonlinear_terms, only: nonlin
    use file_utils, only: open_output_file, close_output_file, get_unused_unit
    use dist_fn_arrays, only: gnew  
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, bparnew
    use agk_layouts, only: idx, proc_id, idx_local, g_lo
    use agk_mem, only: alloc8, dealloc8 
    implicit none
    !Local
    integer :: inn             !1=g, 2=h
    integer :: it              !nakx index
    integer :: ik              !naky index
    integer :: is              !nspec index
    integer :: ie              !negrid index
    integer :: ig              !ntgrid index
    integer :: isgn            !sign of vpar index
    integer :: il              !nlambda index
    integer :: iglo            !iglo g layout index
    integer :: iproc           !Processor number for current meshpoint
    complex, pointer, dimension(:,:,:) :: gv  !v-plane of dist function
    complex, pointer, dimension(:) :: gv2  !v-plane of dist function
    real :: vpar,vperp         !Parallel and perpendicular velocity
    character(10) :: suffix    !Output filename suffix  .vpx##y##
    integer :: unit           !Output File unitnumber
    complex, dimension(2) :: gtmp   !Temporary g (isgn=1,2)

    !Allocate a variable  to contain the slice in velocity space for the given mode
!    if (proc0) write(*,'(a,i4,a,i3,a,i2,a)')'Attempting to allocate a complex variable of size (', 2*nlambda, ',', negrid, ',', nspec, ')'
    
    allocate(gv(2*nlambda,negrid,nspec)) !!! RN> unable to account memory for pointers
    gv=0. 
    allocate(gv2(1:2*nlambda*negrid*nspec)) !!! RN> unable to account memory for pointers
    gv2=0. 

!    if (proc0) write(*,'(a)')'Allocation successful.'

    !Do analysis both for g (inn=1) and h (inn=2)
    do inn=1,2
       if (inn .eq. 2) call g_adjust (gnew, phinew, bparnew, fphi, fbpar, is_g_gnew)

    !Loop over all Fourier modes to get v-space slice for each one
    do ik=1,naky
       do it=1,nakx
          !Skip the constant (0,0) mode if nonlinear run
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          !DEBUG
          if (proc0) write(*,'(a,i2,a,i2,a)')'Analyzing Fourier mode (it,ik)= (',it,',',ik,')'

          !Get slice in v-space-------------------------------------------------------
          if (.false.) then  !Individual send and receive (INEFFICIENT)
          !Loop over species, energy, and pitch angle to assemble data for slice
          do is=1,nspec
             do ie=1,negrid
                do il=1,nlambda
                   !Choose a single point along z
                   ig= 0    !Midplane for the moment (ig runs from -ntgrid to ntgrid)
                   !For this data point (ik,it,is,ie,il) find processor
                   iglo=idx(g_lo,ik,it,il,ie,is)    !Get iglo
                   iproc=proc_id(g_lo,iglo)         !Get proc number for this iglo
                   !If local, get data point and send to gtmp on proc0 
                   if (idx_local (g_lo, ik, it, il, ie, is)) then
                      if (proc0) then 
                         gtmp = gnew(ig,:,iglo)  
                      else
                         call send (gnew(ig,:,iglo), 0)  
                      endif
                   !Receive data point if proc0
                   else if (proc0) then
                      call receive (gtmp, iproc)
                   endif
                   !Put the values into the correct place in gv
                   if (proc0) then
                      gv(ivpar(1,il),ie,is)=cmplx(gtmp(1))
                      gv(ivpar(2,il),ie,is)=cmplx(gtmp(2))
                   endif
                enddo
             enddo
          enddo
          else !Sum_reduce method (EFFICIENT)
          gv=0.  !Initialize slice variable
          !Loop over species, energy, and pitch angle to assemble data for slice
          do is=1,nspec
             do ie=1,negrid
                do il=1,nlambda
                   !Choose a single point along z
                   ig= 0    !Midplane for the moment (ig runs from -ntgrid to ntgrid)
                   !For this data point (ik,it,is,ie,il) find processor
                   iglo=idx(g_lo,ik,it,il,ie,is)    !Get iglo
                   !If local, put data point and send to gtmp on proc0 
                   if (idx_local (g_lo, ik, it, il, ie, is)) then
                      gtmp = gnew(ig,:,iglo)  
                      gv(ivpar(1,il),ie,is)=cmplx(gtmp(1))
                      gv(ivpar(2,il),ie,is)=cmplx(gtmp(2))
                   endif
                enddo
             enddo
          enddo
          !Sum_reduce to proc0 to get all values of array
          !Reshape to 1-D array
          gv2=reshape(gv,(/ 2*nlambda*negrid*nspec /))
          !Call sum_reduce on 1-D array to get all variables to proc0
          call sum_reduce(gv2,0)
          !Reshape back to 3-D array
          gv=reshape(gv2,(/ 2*nlambda,negrid,nspec /))
          endif

          !Output data to file runname.vpx##y##g or  .vpx##y##h   ------------------
          if (proc0) then
            !Create output file suffix
             if (inn .eq. 1) then
                write(suffix,'(a4,i2.2,a1,i2.2,a1)')'.vpx',it,'y',ik,'g'
             elseif (inn .eq. 2) then
                write(suffix,'(a4,i2.2,a1,i2.2,a1)')'.vpx',it,'y',ik,'h'
             endif
             !DEBUG
             write(*,'(a,i2,a,i2,a,a)')'Output for Fourier mode (it,ik)= (',it,',',ik,') to ',suffix
 
             !Open file runname.vpx##y##(g or h)
             call get_unused_unit (unit)
             call open_output_file (unit,suffix)
             
             do is=1,nspec
                do ie=1,negrid
                   isgn=2
                   do il=1,nlambda
                      !Calculate vpar and vperp for this value
                      vpar  = -sqrt(e(ie,is)*max(0.0, 1.0-al(il)))
                      vperp = sqrt(e(ie,is)*al(il))
                      write(unit,'(4es14.6,4i4)')vpar,vperp,gv(ivpar(1,il),ie,is),is,ie,isgn,il
                   enddo
                   isgn=1
                   do il=nlambda,1,-1
                      !Calculate vpar and vperp for this value
                      vpar  = sqrt(e(ie,is)*max(0.0, 1.0-al(il)))
                      vperp = sqrt(e(ie,is)*al(il))
                      write(unit,'(4es14.6,4i4)')vpar,vperp,gv(ivpar(1,il),ie,is),is,ie,isgn,il
                   enddo
                   write(unit,'(a)')' '
                enddo
             enddo

             !Close file
             call close_output_file (unit)
          end if
          
       enddo
    enddo !END Loop over all Fourier modes
    if (inn .eq. 2) call g_adjust (gnew, phinew, bparnew, -fphi, -fbpar, is_g_gnew)
    enddo !Loop over inn (1=g, 2=h)

    !Deallocate
    deallocate (gv)  !!! RN> unable to account memory for pointers
    deallocate (gv2) !!! RN> unable to account memory for pointers

  end subroutine write_vp
!------------------------------------------------------------------------------
! ivpar  GGH 23JAN08
!------------------------------------------------------------------------------
! Function to get sorted single index for the pitch angle and sign for write_vp
  integer function ivpar(isgn,il)
    use le_grids, only: nlambda
    implicit none
    integer :: isgn            !sign of vpar index
    integer :: il              !nlambda index

    if (isgn .eq. 1) then
       ivpar=nlambda+il
    else
       ivpar=nlambda-il+1
    endif
  end function ivpar
!------------------------------------------------------------------------------
! This subroutine only returns epar correctly for linear runs.
  subroutine get_epar (phi, apar, phinew, aparnew, epar)
    use theta_grid, only: ntgrid, delthet, gradpar
    use run_parameters, only: fphi, fapar
    use agk_time, only: dtime
    use kgrids, only: naky, nakx
    complex, dimension(-ntgrid:,:,:) :: phi, apar, phinew, aparnew, epar
    complex :: phi_m, apar_m

    integer :: ig, ik, it

    do ik = 1, naky
       do it = 1, nakx
          do ig = -ntgrid, ntgrid-1
             ! ignoring decentering in time and space for now
             phi_m = 0.5*(phi(ig+1,it,ik)-phi(ig,it,ik) + &
                  phinew(ig+1,it,ik)-phinew(ig,it,ik))*fphi
             apar_m = 0.5*(aparnew(ig+1,it,ik)+aparnew(ig,it,ik) & 
                  -apar(ig+1,it,ik)-apar(ig,it,ik))*fapar
             
! TT: variable step AB3
!             epar(ig,it,ik) = -phi_m/delthet*gradpar - apar_m/dtime
             epar(ig,it,ik) = -phi_m/delthet*gradpar - apar_m/dtime(1)
! <TT
          end do
       end do
    end do    

  end subroutine get_epar

! TT: Gabe's velocity space diagnostics
  subroutine write_vspec (g0, igomega, ekp)

    use mp, only: proc0
    use file_utils, only: open_output_file, flush_output_file
    use le_grids, only: npgrid!, e
    use theta_grid, only: ntgrid
    use kgrids, only: nakx, naky
    use species, only: nspec
    use agk_layouts, only: g_lo!, ie_idx, is_idx
    use agk_time, only: time
    use dist_fn_arrays, only: gp
! TT> for test
    use dist_fn_arrays, only: vperp2
! <TT
    use agk_mem, only: alloc8, dealloc8

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g0
    integer, intent (in) :: igomega
    real, dimension (:,:,:), intent (out) :: ekp ! (nakx,naky,npgrid)
    logical, save :: first=.true.
    integer :: ip!, iglo, ie, is
    integer, save :: unit
    complex, dimension (nakx, naky, npgrid, nspec) :: gkp!, gkpe
! TT> for test
    complex, dimension (:,:,:), allocatable :: gtmp
! <TT

    if (first) then
       ! open file
       if (proc0) then
          call open_output_file (unit, '.vspec')
          write (unit, *) '# time      p          sum_k E(p)      E(0,p)      sum_k E(p) w/ exp'
       end if
       first=.false.
    end if

    if (test) then

       allocate (gtmp(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)); call alloc8(c3=gtmp,v="gtmp")
       ! test functions; exp(vperp2)/2.0 is the normalization factor
       ! H(1/v) = 1/p
       gtmp(igomega,:,:) = spread(exp(vperp2)/sqrt(vperp2),1,2) / 2.0
       ! H(e^(-2v^2)) = e^(-p^2/8)/4
!       gtmp(igomega,:,:) = spread(exp(-vperp2),1,2) / 2.0
       call Hankel_transform (gtmp, igomega, gkp)
       do ip=1, npgrid
          write (unit,*) gp(ip), abs(gkp(1,1,ip,1))
       end do
!       deallocate (gtmp)
       call dealloc8(c3=gtmp,v="gtmp")

    else

       ! gkp is a Hankel transform of g
       call Hankel_transform (g0, igomega, gkp)

       ! gkpe is a Hankel transform of g/F_0
!!$       allocate (gtmp(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
!!$       do iglo = g_lo%llim_proc, g_lo%ulim_proc
!!$          ie = ie_idx(g_lo,iglo)
!!$          is = is_idx(g_lo,iglo)
!!$          gtmp(:,:,iglo) = g0(:,:,iglo) * exp(e(ie,is))
!!$       end do
!!$       call Hankel_transform (gtmp, igomega, gkpe)
!!$       deallocate (gtmp)

       if (proc0) then ! write .vspec file
          ! get ekp --- last index picking up 1st species
          ekp(:,1,:) = real( gkp(:,1,:,1)*conjg(gkp(:,1,:,1)) ) * spread(gp,1,nakx)
          ekp(:,2:,:) = real( gkp(:,2:,:,1)*conjg(gkp(:,2:,:,1)) ) * 0.5 &
               * spread( spread(gp,1,naky-1), 1, nakx )
          ! sum w.r.t. kx & ky
          do ip=1, npgrid
             write (unit,'(2f10.5,3es15.5e3)') time, gp(ip), &
                  & sum(ekp(:,:,ip)), ekp(1,1,ip)!!$, &
!!$                  & gp(ip) * ( sum( real( gkpe(:,1,ip,1)*conjg(gkpe(:,1,ip,1)) ) ) &
!!$                  & + sum( real( gkpe(:,2:,ip,1)*conjg(gkpe(:,2:,ip,1)) ) ) * 0.5 )
          end do
          write (unit,*)
          call flush_output_file (unit)
       end if

    end if

  end subroutine write_vspec

  subroutine Hankel_transform (g0, igomega, gkp, all)
    use le_grids, only: negrid, e, integrate_species, vgrid, nesub, npgrid
    use kgrids, only: aky
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: vperp2, gp
    use agk_layouts, only: g_lo
    use spfunc, only: j0
    use agk_mem, only: alloc8, dealloc8
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g0
    integer, intent(in) :: igomega
    complex, dimension (:,:,:,:), intent (out) :: gkp
    integer, optional, intent (in) :: all

    logical, save :: pzero=.false.
    integer :: ip, iv, iglo
    real :: jold, jnew, v, delv=0.05, eps=1.e-5, dky
    real, dimension (:,:), allocatable, save :: j0p
!    complex, dimension (2, g_lo%llim_proc:g_lo%ulim_alloc, negrid) :: g1
    complex, dimension (:,:,:), allocatable :: g1

    if (.not.allocated(gp)) then

       if (pzero) then
          ! determine p in terms of Bessel zeroes
          allocate (gp(npgrid)); call alloc8(r1=gp,v="gp")
          iv = 0
          jnew = 1.0
          do ip=1, npgrid
             do  ! infinite loop
                iv = iv + 1
                jold = jnew
                v = delv*iv
                jnew = j0(v)
                if (jnew*jold < 0.) exit
             end do
             call root (v, eps)
             gp(ip) = v
          end do
          if (vgrid) then
             gp = gp / sqrt(e(nesub,1))  ! got rescaled p(1:npgrid)
          else
             gp = gp / sqrt(e(negrid,1))  ! got rescaled p(1:npgrid)
          end if
       else
          ! define p simply by ky's
          dky = aky(2) - aky(1)
          allocate (gp(npgrid)); call alloc8(r1=gp,v="gp")
          if (dky > 0.1) then
             ! linear grid
             gp = (/ (ip*dky, ip=1, npgrid) /)
          else
             ! log grid
             gp = (/ (exp( (ip-1)*(log(10.0)-log(dky))/(npgrid-1) + log(dky) ), &
                  & ip=1, npgrid) /)
          end if
       end if

       ! create J0(p*vperp) array
       allocate (j0p(g_lo%llim_proc:g_lo%ulim_alloc,npgrid)); call alloc8(r2=j0p,v="j0p")
       do ip=1, npgrid
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             j0p(iglo,ip) = j0(sqrt(vperp2(iglo)) * gp(ip))
          end do
       end do ! got j0p(g_lo,1:npgrid) array
    end if

    ! Hankel transform
    allocate (g1(2,g_lo%llim_proc:g_lo%ulim_alloc,npgrid)); call alloc8(c3=g1,v="g1")
    g1 = spread(g0(igomega,:,:),3,npgrid) * spread(j0p,1,2)
    ! TT: This should actually be called integrate_moments, but...
    if (present(all)) then
       call integrate_species (g1, gkp, all)
    else
       call integrate_species (g1, gkp)
    end if
    call dealloc8 (c3=g1,v="g1")

  contains
    subroutine root (v, eps)
      use file_utils, only: error_unit
      use spfunc, only: j0, j1
      integer :: stp, maxstp=30
      real, intent (inout) :: v
      real, intent (in) :: eps
      do stp=1, maxstp
         v = v + j0(v) / j1(v) / v ! remember that j1(v) = J1(v)/v
         if (abs(j0(v))<eps) exit
      end do
      if (stp>maxstp) write(error_unit(),*) &
           & 'WARNING: max step exceeded in root finding of Hankel_transform'
    end subroutine root
  end subroutine Hankel_transform

  subroutine dot (a, anew, adot, fac)

! Get the time derivative of a field.
! 
    use agk_time, only: dtime
    use kgrids, only: naky, nakx
    use theta_grid, only: ntgrid

    implicit none
    complex, intent (in), dimension (-ntgrid:,:,:) :: a, anew
    complex, intent (out), dimension (-ntgrid:,:,:) :: adot
    real, intent (in) :: fac
    real :: dtinv
    integer :: ig, it, ik

! TT: variable step AB3
!    dtinv = 1./dtime
    dtinv = 1./dtime(1)
! <TT
    do ik=1,naky
       do it=1,nakx
          do ig=-ntgrid,ntgrid
             adot(ig,it,ik) = fac*(anew(ig,it,ik) - a(ig,it,ik))*dtinv
          end do
       end do
    end do
    
  end subroutine dot

  subroutine init_eq
    use run_parameters, only: store_eq
    use theta_grid, only: ntgrid
    use agk_layouts, only: g_lo
    use species, only: nspec
    use kgrids, only: nakx, naky
    use dist_fn_arrays, only: g_eq
    use dist_fn_arrays, only: dens_eq, ux_eq, uy_eq, uz_eq
    use dist_fn_arrays, only: pxx_eq, pyy_eq, pzz_eq, pxy_eq, pyz_eq, pzx_eq
    use agk_mem, only: alloc8

    logical, save :: initialized = .false.

    if (.not.store_eq .or. initialized) return
    
    allocate (g_eq(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc)) ; call alloc8 (c3=g_eq, v="g_eq")
    g_eq=0.
    
    allocate(dens_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=dens_eq,v='dens_eq')
    allocate(ux_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=ux_eq,v='ux_eq')
    allocate(uy_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=uy_eq,v='uy_eq')
    allocate(uz_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=uz_eq,v='uz_eq')
    allocate(pxx_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=pxx_eq,v='pxx_eq')
    allocate(pyy_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=pyy_eq,v='pyy_eq')
    allocate(pzz_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=pzz_eq,v='pzz_eq')
    allocate(pxy_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=pxy_eq,v='pxy_eq')
    allocate(pyz_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=pyz_eq,v='pyz_eq')
    allocate(pzx_eq(-ntgrid:ntgrid,nakx,naky,nspec)) ; call alloc8(c4=pzx_eq,v='pzx_eq')
    dens_eq=0.
    ux_eq=0.; uy_eq=0.; uz_eq=0.
    pxx_eq=0.; pyy_eq=0.; pzz_eq=0.; pxy_eq=0.; pyz_eq=0.; pzx_eq=0.

    initialized = .true.
  end subroutine init_eq

  subroutine init_mom_coeff
    use agk_layouts, only: g_lo
    use species, only: nspec
    use kgrids, only: nakx, naky, kperp2
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: aj0,aj1vp2
    use dist_fn_arrays, only: vpa, vperp2
    use le_grids, only: integrate_moment
    use species, only: spec
    use agk_mem, only: alloc8, dealloc8
    integer :: i
    integer :: isgn, iglo
    integer :: it,ik,is
    complex, allocatable :: coeff0(:,:,:,:)
    complex, allocatable :: gtmp(:,:,:)
    real :: wgt
!    logical, parameter :: analytical = .false.
    logical, parameter :: analytical = .true.
    logical, save :: initialized=.false.
    real :: bsq

    if(initialized) return

    if(.not.allocated(mom_coeff)) then
       allocate(mom_coeff(nakx,naky,nspec,ncnt_mom_coeff)); call alloc8(r4=mom_coeff,v="mom_coeff")
    end if
    if(.not.allocated(mom_coeff_npara)) then
       allocate(mom_coeff_npara(nakx,naky,nspec)); call alloc8(r3=mom_coeff_npara,v="mom_coeff_npara")
    end if
    if(.not.allocated(mom_coeff_nperp)) then
       allocate(mom_coeff_nperp(nakx,naky,nspec)); call alloc8(r3=mom_coeff_nperp,v="mom_coeff_nperp")
    end if
    if(.not.allocated(mom_coeff_tpara)) then
       allocate(mom_coeff_tpara(nakx,naky,nspec)); call alloc8(r3=mom_coeff_tpara,v="mom_coeff_tpara")
    end if
    if(.not.allocated(mom_coeff_tperp)) then
       allocate(mom_coeff_tperp(nakx,naky,nspec)); call alloc8(r3=mom_coeff_tperp,v="mom_coeff_tperp")
    end if
    if(.not.allocated(mom_shift_para)) then
       allocate(mom_shift_para(nakx,naky,nspec)); call alloc8(r3=mom_shift_para,v="mom_shift_para")
    end if
    if(.not.allocated(mom_shift_perp)) then
       allocate(mom_shift_perp(nakx,naky,nspec)); call alloc8(r3=mom_shift_perp,v="mom_shift_perp")
    end if

    mom_coeff(:,:,:,:)=0.
    mom_coeff_npara(:,:,:)=0. ; mom_coeff_nperp(:,:,:)=0.
    mom_coeff_tpara(:,:,:)=0. ; mom_coeff_tperp(:,:,:)=0.
    mom_shift_para(:,:,:)=0.  ; mom_shift_perp(:,:,:)=0.

    allocate(coeff0(-ntgrid:ntgrid,nakx,naky,nspec)); call alloc8(c4=coeff0,v="coeff0")
    allocate(gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)); call alloc8(c3=gtmp,v="gtmp")
    coeff0(:,:,:,:)=cmplx(0.,0.)
    gtmp(:,:,:)=cmplx(0.,0.)

    if (analytical) then
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                bsq=.25*spec(is)%smz**2*kperp2(it,ik)
                mom_coeff(it,ik,is,1) = exp(-bsq)
                mom_coeff(it,ik,is,2) = exp(-bsq) *.5
                mom_coeff(it,ik,is,3) = exp(-bsq) *(1.-bsq)
                mom_coeff(it,ik,is,4) = exp(-bsq) *.75
                mom_coeff(it,ik,is,5) = exp(-bsq) *(1.-bsq)*.5
                mom_coeff(it,ik,is,6) = exp(-bsq) *.5
                mom_coeff(it,ik,is,7) = exp(-bsq) *.25
                mom_coeff(it,ik,is,8) = exp(-bsq) *(1.-.5*bsq)
                mom_shift_perp(it,ik,is)=1.-bsq !A20/A00
             end do
          end do
       end do
       mom_shift_para(:,:,:)=.5 !A02/A00

       mom_coeff_npara(:,:,:)=1. !2*A02/A00
       mom_coeff_nperp(:,:,:)=1. !2*B20/A00

       mom_coeff_tperp(:,:,:)=0. ! (A22-S_perp A02)/(B40-S_perp B20)
       mom_coeff_tpara(:,:,:)=0. ! (B22-S_para B20)/(A04-S_para A02)
    else
       do i = 1, ncnt_mom_coeff
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do isgn = 1,2
                if(i==1) wgt=aj0(iglo) ! A00
                if(i==2) wgt=aj0(iglo)*vpa(isgn,iglo)**2 ! A02
                if(i==3) wgt=aj0(iglo)*vperp2(iglo) ! A20
                if(i==4) wgt=aj0(iglo)*vpa(isgn,iglo)**4 ! A04
                if(i==5) wgt=aj0(iglo)*vpa(isgn,iglo)**2*vperp2(iglo) ! A22
                if(i==6) wgt=.5*aj1vp2(iglo) ! B20
                if(i==7) wgt=.5*aj1vp2(iglo)*vpa(isgn,iglo)**2 !B22
                if(i==8) wgt=.5*aj1vp2(iglo)*vperp2(iglo) ! B40
                gtmp(-ntgrid:ntgrid,isgn,iglo) = wgt*cmplx(1.,0.)
             end do
          end do
          call integrate_moment(gtmp,coeff0,1)
          where(real(coeff0(0,1:nakx,1:naky,1:nspec)) == 0.)
             mom_coeff(1:nakx,1:naky,1:nspec,i)=1.
          elsewhere
             mom_coeff(1:nakx,1:naky,1:nspec,i)= &
                  & coeff0(0,1:nakx,1:naky,1:nspec)
          end where
       end do

       mom_shift_para(:,:,:)=mom_coeff(:,:,:,2)/mom_coeff(:,:,:,1) !A02/A00
       mom_shift_perp(:,:,:)=mom_coeff(:,:,:,3)/mom_coeff(:,:,:,1) !A20/A00

       mom_coeff_npara(:,:,:)=2.*mom_coeff(:,:,:,2)/mom_coeff(:,:,:,1) !2*A02/A00
       mom_coeff_nperp(:,:,:)=2.*mom_coeff(:,:,:,6)/mom_coeff(:,:,:,1) !2*B20/A00

       mom_coeff_tperp(:,:,:)= & ! (A22-S_perp A02)/(B40-S_perp B20)
            & (mom_coeff(:,:,:,5)-mom_shift_perp(:,:,:)*mom_coeff(:,:,:,2)) / &
            & (mom_coeff(:,:,:,8)-mom_shift_perp(:,:,:)*mom_coeff(:,:,:,6))
       mom_coeff_tpara(:,:,:)= & ! (B22-S_para B20)/(A04-S_para A02)
            & (mom_coeff(:,:,:,7)-mom_shift_para(:,:,:)*mom_coeff(:,:,:,6)) / &
            & (mom_coeff(:,:,:,4)-mom_shift_para(:,:,:)*mom_coeff(:,:,:,2))
    endif

!    deallocate(gtmp,coeff0)
    call dealloc8(c3=gtmp,v="gtmp")
    call dealloc8(c4=coeff0,v="coeff0")
    
    initialized=.true.
  end subroutine init_mom_coeff

  subroutine pass_nwrite (nwrite_in, navg_in)

    implicit none

    integer, intent (in) :: nwrite_in, navg_in

    nwrite = nwrite_in ; navg = navg_in

  end subroutine pass_nwrite

end module dist_fn

