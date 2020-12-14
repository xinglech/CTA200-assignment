!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> This module provides storage and equations of motion for a relativistic scalar evolving in one spatial dimension
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  implicit none

  integer :: nLat, nFld, nVar
  real(dl), dimension(:), allocatable, target :: yvec

  real(dl) :: len, dx, dk
  real(dl) :: lambda !, m2eff, m2eff2
  real(dl), dimension(:), allocatable :: m2eff
  
#ifdef FOURIER
  type(transformPair1D) :: tPair
#endif
  
contains

  ! Fix the paramter nature above (in particular the declaration of yvec
  subroutine set_lattice_params(n,l,nf)
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: l
    nLat = n; len = l !; nFld = nf
    nFld = 2  ! Hardcoding number of fields, since otherwise things break
    
    nVar = 2*nFld*nLat + 1; allocate(yvec(1:nVar))
    allocate(m2eff(1:nFld))
    dx = len / dble(nLat); dk = twopi/len
  end subroutine set_lattice_params

  ! Add appropriate subroutine calls here to potential derivs, etc. here
  subroutine set_model_params(lam,m2)
    real(dl), intent(in) :: lam, m2
    
    lambda = lam
    m2eff(1) = 2._dl*(1._dl-lam); m2eff(2) = 0._dl
  end subroutine set_model_params

  function get_m2eff(fInd) result(m2)
    real(dl) :: m2
    integer, intent(in) :: fInd
    select case (fInd)
    case (1)
       m2 = 2._dl*(1._dl-lambda)
       ! m2 = -1._dl I cannot have negative mass?
    case(2)
       m2 = 0._dl
    case default
       m2 = 0._dl
    end select
  end function get_m2eff
  
  real(dl) elemental function v(phi, chi)
    real(dl), intent(in) :: phi, chi
    v = 0.25_dl*(phi**2-1._dl)**2+lambda*((1._dl/3._dl)*phi**3-phi+2._dl/3._dl) 
!      + 0.25_dl*(chi**2-1._dl)**2 + lambda*((1._dl/3._dl)*chi**3-chi+2._dl/3._dl)
  end function v

  real(dl) elemental function vp(phi)
    real(dl), intent(in) :: phi
    vp =  (phi+lambda)*(phi**2-1._dl)
  end function vp

!  real(dl) elemental function vpp(phi)
!    real(dl), intent(in) :: phi
!    vpp = 3._dl*phi**2+2*phi*lambda-1._dl
!  end function vpp
  
  function phi_fv()
    real(dl), dimension(1:nFld) :: phi_fv
    phi_fv = (/ -0.88784039_dl, 0._dl /)
  end function phi_fv

! Lechun: introduce chi and v2(chi), vp2(chi), vpp2(chi)
! where did I call m2eff in #set_model_params(lam,m2)? Ln 50
  real(dl) elemental function v2(phi, chi)
    real(dl), intent(in) :: phi, chi
    v2 = lambda*chi
  end function v2

  real(dl) elemental function vp2(phi,chi)
    real(dl), intent(in) :: phi, chi
    vp2 =  lambda
  end function vp2

  real(dl) elemental function vpp2(phi,chi)
    real(dl), intent(in) :: phi, chi
    vpp2 = 0._dl
  end function vpp2

  !>@brief
  !> Compute the derivatives of the scalar field in the effective time-independent potential
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
#ifdef DISCRETE
    real(dl) :: lNorm 
    lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl 
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -((yc(FLD)+lambda)*(yc(FLD)**2-1._dl))
    yp(FLD2) = yc(DFLD2)
    yp(DFLD2) = 0._dl !-((yc(FLD2)+lambda)*(yc(FLD2)**2-1._dl))
    
#ifdef DIFF
#ifdef DISCRETE
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )

    yp(3*nlat+1) = yp(3*nlat+1) + lNorm * ( yc(3*nlat) - 2._dl*yc(2*nLat+1) + yc(2*nLat+2) )
    yp(4*nlat) = yp(4*nlat) + lNorm * ( yc(3*nlat-1) - 2._dl*yc(3*nlat) + yc(2*nLat+1) )
    
    yp(3*nlat+2:4*nlat-1) = yp(3*nlat+2:4*nlat-1) + lNorm*( yc(2*nLat+1:3*nlat-2) - 2._dl*yc(2*nLat+2:3*nlat-1) + yc(2*nLat+3:3*nlat) )
#else
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)

    tPair%realSpace(:) = yc(FLD2)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD2) = yp(DFLD2) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs
  
end module eom
