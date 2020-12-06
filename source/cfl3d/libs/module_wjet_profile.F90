module module_wjet_profile
  use bspline_module
  use module_kinds, only: wp, get_free_unit
  implicit none
  type(bspline_1d),pointer,save :: bsp_wjet
contains
  subroutine read_wjet_profile(iflag)
    implicit none
    integer :: iounit, np, order, i, iflag
    logical :: exists
    common /bcwjet/ jet_profile
    integer :: jet_profile
    real(wp), allocatable :: t(:), v(:)
    if (jet_profile == 0) return
    iounit = get_free_unit()
    open(unit=iounit, file='walljet.profile', status='old')
    read(iounit, *) ! title
    read(iounit, *) np, order
    allocate(t(np), v(np))
    do i=1,np
      read(iounit, *) t(i), v(i)
    enddo
    bsp_wjet%initialize(t, v, order, iflag, .true.)
  
  end subroutine read_wjet_profile

  subroutine wjet_vel(t, v, iflag)
    implicit none
    real(wp) :: t, v
    integer :: iflag
    call bsp_wjet%evaluate(t, 0, v, iflag)
    v = max(0.0, min(1.0, v))
  end subroutine
  
end module