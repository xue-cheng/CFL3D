module module_kinds
  use,intrinsic :: iso_fortran_env
  implicit none
  private
#ifdef DBLE_PRECSN
  integer,parameter,public :: wp = real64  !! Real precision
#else
  integer,parameter,public :: wp = real32  !! Real precision
#endif
  public :: get_free_unit
  contains 
  function get_free_unit() result(u)
    implicit none
    integer :: u
    logical :: e
    do u=100,9999
      inquire(unit=u, opened=e)
      if (.not.e) exit
    enddo
  end function get_free_unit
end module module_kinds