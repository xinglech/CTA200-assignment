module Procedures
  use constants, only : dl
  implicit none

  abstract interface
    real(dl) function null_spec(k)
      real(dl), intent(in) :: k
    end function null_spec
  end interface

contains

  subroutine

end module Procedures
