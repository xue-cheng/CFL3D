module module_kinds
  use,intrinsic :: iso_fortran_env
  implicit none
  private
#ifdef DBLE_PRECSN
  integer,parameter,public :: wp = real64  !! Real precision
#else
  integer,parameter,public :: wp = real32  !! Real precision
#endif
end module module_kinds