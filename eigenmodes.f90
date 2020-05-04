!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> Canadian Institute for Theoretical Astrophysics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module eigenmodes
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi

  implicit none

contains

  !>@brief
  !> Compute the bubble eigenmodes from an "initial" late time profile
  subroutine eigenmodes(fld,eval,evec,spec)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), dimension(:), intent(out) :: eval
    real(dl), dimension(:,:), intent(out) :: evec
    logical, intent(in) :: spec

    ! 1. Start by constructing the Laplacian matrix, provide options for way this is done
    ! 2. Compute the (diagonal) effective mass matrix using the field profile and potential
    ! 3. Use LAPACK to compute the eigenvalues and eigenvectors
  end subroutine eigenmodes

  subroutine initialise_eigenmodes(fld)
    real(dl), dimension(:), intent(in) :: fld

    ! 1. Call subroutine to obtain the eigenmodes
    ! 2. Generate random deviates to multiply onto eigenmodes
    ! 3. Add realisation of eigenmodes onto original field
  end subroutine initialise_eigenmodes

end module eigenmodes
