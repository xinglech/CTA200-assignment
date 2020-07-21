#include "macros.h"

program Scalar_1D
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use utils, only : newunit
  use gaussianRandomField  ! remove this to fluctuations module
  use Fluctuations
  use eom
  use integrator
  use bubble_extraction, only : count_bubbles, mean_cos

  implicit none

  real(dl), dimension(:,:), pointer :: fld
  real(dl), pointer :: time
  real(dl) :: dtout_, dt_  ! Figure out how to get rid of this horrible nonlocality (used in output, and elsewhere)

  integer :: i

  real(dl) :: alph, t_cross
  integer :: n_cross
  integer :: nSamp
  
  type SimParams
     real(dl) :: dx, dt, dtout
     integer :: nLat
  end type SimParams
  type(SimParams) :: sim

  type ScalarLattice
     real(dl), dimension(:,:), allocatable :: flds
  end type ScalarLattice
  type(ScalarLattice) :: simulation

! What needs fixing, set phi0 out here, allow m^2 to vary from vacuum value, etc.
  call set_lattice_params(1024,50._dl,2)
! (1024,50._dl,2) = (n, l, nf) = (nLat, len, nFld) 
  call set_model_params(0.5_dl,100._dl)  ! A default for the double well
  ! call set_model_params(1.2_dl,1._dl) = (lam, m2)
  
  fld(1:nLat,1:4) => yvec(1:2*nLat*nFld)
  time => yvec(2*nLat*nFld+1)
  alph = 4._dl; n_cross = 4

  call initialize_rand(87,18)  ! Seed for random field generation.  Adjust to make a new field realisation
  call setup(nVar)

  nSamp = 1
  do i=1,nSamp
     call initialise_fields(fld,nLat/4+1,0.5*twopi)
     call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross)
  enddo
 
!  call forward_backward_evolution(0.4_dl/omega,10000,100)
  
contains

  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup
  !> before starting the program
  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
    call initialize_transform_1d(tPair,nLat)  ! nonlocality
  end subroutine setup

  subroutine time_evolve(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no
    integer :: i,j, outsize, nums
    integer :: b_file

    open(unit=newunit(b_file),file='bubble-count.dat')
    
    print*,"dx is ", dx, "dt is ",dt, "dx/dt is ",dx/dt
    if (dt > dx) print*,"Warning, violating Courant condition"
    
    outsize = ns/no; nums = ns/outsize
    print*,"dt out is ",dt*outsize
    dt_ = dt; dtout_ = dt*outsize  ! Used here again
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_fields(fld)
       write(b_file,*) count_bubbles(fld(:,1)), mean_cos(fld(:,1))
    enddo
    write(b_file,*)
  end subroutine time_evolve

  subroutine initialise_fields(fld,kmax,phi,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax
    real(dl), intent(in), optional :: phi
    integer, intent(in), optional :: klat

    integer :: i; real(dl) :: dt, theta
    integer :: kc, nn
    real(dl) :: phiL
    
    nn = size(fld(:,1))/2+1
    kc = nn; if (present(klat)) kc = klat
    phiL = 0.5_dl*twopi; if (present(phi)) phiL = phi

    call initialise_mean_fields(fld)
    yvec(size(fld)+1) = 0._dl ! Add a tcur pointer here
    call initialize_vacuum_fluctuations(fld(:,1:2),len,get_m2eff(1),kmax,phiL,kc)
    call initialize_vacuum_fluctuations(fld(:,3:4),len,get_m2eff(2),kmax,phiL,kc)
  end subroutine initialise_fields

  function light_cross_time(len) result(tmax)
    real(dl), intent(in) :: len
    real(dl) :: tmax
    tmax = 0.5_dl*len
  end function light_cross_time

  function convert_t_to_nstep(dt,dtout,tend) result(ns)
    real(dl), intent(in) :: dt,dtout, tend
    integer, dimension(1:2) :: ns

    ns(1) = int(tend/dt)
    ns(2) = int(dtout/dt)
  end function convert_t_to_nstep

  subroutine convert_tstep_to_int(dt,dtout,tend,ns,nout)
    real(dl), intent(in) :: dt, dtout, tend
    integer, intent(out) :: ns, nout

    ns = int(tend/dt)
    nout = int(dtout/dt)
  end subroutine convert_tstep_to_int

  subroutine extend_grid(fld_old,fld_new,us)
    real(dl), dimension(:,:), intent(in) :: fld_old
    real(dl), dimension(:,:), intent(out) :: fld_new
    integer, intent(in) :: us  ! upsample ratio
  end subroutine extend_grid

  subroutine resample(fld_old,fld_new,us)
    real(dl), dimension(:,:), intent(in) :: fld_old
    real(dl), dimension(:,:), intent(out) :: fld_new
    integer, intent(in) :: us
  end subroutine resample
  
  !>@brief
  !> Evolve a collection of ns field trajectories holding the long-wavelength part of the field fixed while varying the short wavelengths
  subroutine vary_high_k_modes(phi_l,ns)
    real(dl), dimension(:,:), intent(in) :: phi_l
    integer, intent(in) :: ns

    real(dl), dimension(1:nlat) :: df
    integer :: i
!    call initialise_fields(phi_l,nlat/8)
    do i=1,ns
       ! call generate_1dGRF(df)
       fld(:,1) = phi_l(:,1) + df
       ! call generate_1dGRF(df)
       fld(:,2) = phi_l(:,2) + df
       !call time_evolve()
    enddo
  end subroutine vary_high_k_modes

  !>@brief
  !> Evolve a collection of ns field trajectories holding the short-wavelength part of the field fixed while varying the long wavelengths
  subroutine vary_low_k_modes(phi_s,ns)
    real(dl), dimension(:), intent(in) :: phi_s
    integer, intent(in) :: ns
  end subroutine vary_low_k_modes
  
  ! Fix nLat nonlocality here
  subroutine forward_evolution(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no

    call initialise_fields(fld,nLat/8)
    call setup(nVar)
    call output_fields(fld)
    call time_evolve(dt,ns,no)
    call write_checkpoint(fld,time,dx,nLat)
  end subroutine forward_evolution
  
  subroutine forward_backward_evolution(dt,ns,no,amp)
    real(dl), intent(in) :: dt
    integer,intent(in) :: ns,no
    real(dl), intent(in), optional :: amp
    
!    call initialize_rand(72,18)
    call initialise_fields(fld,nLat/8)
    call setup(nVar)
    call output_fields(fld)
    call time_evolve(dt,ns,no)
    call write_checkpoint(fld,time,dx,nLat)

    ! now time reverse by flipping sign of time derivative
    fld(:,2) = -fld(:,2)
    if (present(amp)) call initialise_new_fluctuations(fld,amp)  

    call time_evolve(dt,ns,no)
  end subroutine forward_backward_evolution

  !>@brief
  !> Reverse the time flow of a simulation by flipping the sign of phidot
  subroutine reverse_time(fld)
    real(dl), intent(inout), dimension(:,:) :: fld
    fld(:,2) = -fld(:,2)
  end subroutine reverse_time
  
  !>@brief
  !> Initialise the field to have mean value given by the false vacuum and no mean velocity 
  subroutine initialise_mean_fields(fld)
    real(dl), dimension(:,:), intent(out) :: fld
    real(dl), dimension(1:nFld) :: fv
    integer :: i, ii
    fv = phi_fv()
    do i=1,nFld
       ii = 2*(i-1)+1
       fld(:,ii) = fv(i)
       fld(:,ii+1) = 0._dl
    enddo
  end subroutine initialise_mean_fields
  
  !!!! Fix this thing up
  !>@brief
  !> Initialise the field fluctuations
  subroutine initialise_new_fluctuations(fld,amp)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: amp
    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j
    
    spec = 0._dl
    do i=2,nLat/2
       spec(i) = (0.5_dl)**0.5 / (sqrt(len))
    enddo
    call generate_1dGRF(df,spec(1:128),.false.)
    fld(:,1) = fld(:,1) + amp*df(:)
    call generate_1dGRF(df,spec(1:128),.false.)
    fld(:,2) = fld(:,2) + amp*df(:)
  end subroutine initialise_new_fluctuations

! Add smoothing from my other repository
  subroutine output_fields(fld)
    real(dl), dimension(:,:), intent(in) :: fld
    logical :: o; integer :: i, ii
    integer, save :: oFile
    real(dl), dimension(1:nLat,1:nFld) :: gsq, gsq_fd

!    if (.true.) return
    
    inquire(file='fields.dat',opened=o)
    if (.not.o) then
       open(unit=oFile,file='fields.dat')
       write(oFile,*) "# Lattice Parameters"
       write(oFile,*) "# n = ",nLat," dx = ",dx
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_
       write(oFile,*) "#"
       write(oFile,*) "# Phi  PhiDot  Chi  ChiDot  GradPhi^2 (FD)  V(phi)  GradPhi^2 (Spec) V_quad"
    endif

    gsq_fd = 0.
    do i=1,nFld
       ii = 2*(i-1)+1
       gsq_fd(1,i) = 0.5_dl*( (fld(nLat,ii)-fld(1,ii))**2+(fld(2,ii)-fld(1,ii))**2 )
       gsq_fd(nLat,i) = 0.5_dl*( (fld(nLat-1,ii)-fld(nLat,ii))**2+(fld(nLat,ii)-fld(1,ii))**2  )
       gsq_fd(2:nLat-1,i) = 0.5_dl*( (fld(1:nLat-2,ii)-fld(2:nLat-1,ii))**2+(fld(3:nLat,ii)-fld(2:nlat-1,ii))**2 )
    enddo
    gsq_fd = gsq_fd / dx**2
      
#ifdef FOURIER
    do i=1,nFld
       tPair%realSpace(:) = fld(1:nLat,2*(i-1)+1)
       call gradsquared_1d_wtype(tPair,dk)
       gsq(:,i) = tPair%realSpace(:)
    enddo
#else
    gsq = 0._dl  ! tPair isn't created unless doing Fourier transforms
#endif
    ! Fix this if I change array orderings
    do i=1,size(fld(:,1))
       write(oFile,*) fld(i,:), gsq_fd(i,:), v(fld(i,1), fld(i,3)), gsq(i,:)
    enddo
    write(oFile,*)
    
    print*,"conservation :", &
         0.5_dl*sum( 0.5_dl*gsq) + sum(v(fld(:,1), fld(:,3))) + 0.5_dl*sum(fld(:,2::2)**2), &
         0.5_dl*sum(gsq_fd) + sum(v(fld(:,1),fld(:,3))) + 0.5_dl*sum(fld(:,2::2)**2)
  end subroutine output_fields

  !>@brief
  !> Write a checkpoint file with all information for restarting the simulation.
  subroutine write_checkpoint(fld,tcur,dx,n)
    real(dl), intent(in) :: fld(:,:), tcur, dx
    integer, intent(in) :: n

    integer :: fn, i
    open(unit=newunit(fn),file='flds.chk')
    write(fn,*) n, dx, tcur
    do i=1,n
       write(fn,*) fld(i,:)
    enddo
    close(fn)
  end subroutine write_checkpoint

  !>@brief
  !> Read in a previously produced checkpoint file to initialise the simulation.
  subroutine read_checkpoint(fld,tcur,nLat,fName)
    real(dl), intent(out) :: fld(:,:), tcur
    integer, intent(out) :: nLat
    character(*), intent(in) :: fName

    integer :: fn
    integer :: i, n
    real(dl) :: dx, tc
    
    open(unit=newunit(fn),file=fName)
    read(fn,*) n, dx, tc
    print*,"Reading checkpoint file: N = ",n," dx = ",dx," t = ",tc
    ! Add a check that n matches the parameter used for the sim
    do i=1,n
       read(fn,*) fld(i,:)
    enddo
    close(fn)
  end subroutine read_checkpoint

  subroutine initialise_from_file(fName,n)
    character(*), intent(in) :: fName
    integer, intent(in) :: n
    integer :: i; real(dl) :: f,df
    integer :: fNum

    open(unit=newunit(fNum),file=fName)
    do i=1,n
       read(fNum,*) f,df; fld(i,1) = f; fld(i,2) = df
    enddo
    close(fNum)
  end subroutine initialise_from_file
  
end program Scalar_1D
