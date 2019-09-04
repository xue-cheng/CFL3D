module module_gust
  use bspline_module
  use bspline_kinds_module, only: wp
  implicit none

  type gust_profile
    integer:: nx, ny, nz, kx, ky, kz
    real(wp):: xmin,xmax,ymin,ymax,zmin,zmax
    real(wp), allocatable:: x(:), y(:), z(:)
    real(wp), allocatable:: u(:,:,:), v(:,:,:), w(:,:,:)
  end type gust_profile

  type(gust_profile),save :: gprf
  type(bspline_3d),save :: gitp_u, gitp_v, gitp_w

contains

  subroutine read_gust_profile(iflag)
    implicit none
    integer, intent(out)::iflag
    integer :: ug = 955
    logical :: exists
    integer :: i,j,k
    integer :: ip3dgrd,ialph
    common /igrdtyp/ ip3dgrd,ialph
    if (ialph .ne. 0) then 
      iflag = 3
      return 
    end if
    iflag = 0 ! succeed
    inquire(file="gust_profile.inp", exist=exists)
    if (.not. exists) then
      iflag = 1 ! file not exists
      return
    end if
    inquire(unit=ug, opened=exists)
    if (exists) then
      iflag = 2 ! unit has been opened
      return
    end if
    open(unit=ug, file="gust_profile.inp", status='old')
    read(ug, *) !title
    read(ug, *) gprf%nx,gprf%ny,gprf%nz,gprf%kx,gprf%ky,gprf%kz
    call check_dir_order(gprf%nx,gprf%ny,gprf%nz,gprf%kx,gprf%ky,gprf%kz,iflag)
    if (iflag .ne. 0) return
    if (ialph .eq. 0) then 
      ! alpha in x-z plane
      allocate(gprf%x(gprf%nx), gprf%y(gprf%ny), gprf%z(gprf%nz))
      allocate(gprf%u(gprf%nx,gprf%ny,gprf%nz))
      allocate(gprf%v(gprf%nx,gprf%ny,gprf%nz))
      allocate(gprf%w(gprf%nx,gprf%ny,gprf%nz))
      ! read profile data
      read(ug, *) (gprf%x(i), i=1,gprf%nx)
      read(ug, *) (gprf%y(j), j=1,gprf%ny)
      read(ug, *) (gprf%z(k), k=1,gprf%nz)
      read(ug, *) (((gprf%u(i,j,k), i=1,gprf%nx), j=1,gprf%ny), k=1,gprf%nz)
      read(ug, *) (((gprf%v(i,j,k), i=1,gprf%nx), j=1,gprf%ny), k=1,gprf%nz)
      read(ug, *) (((gprf%w(i,j,k), i=1,gprf%nx), j=1,gprf%ny), k=1,gprf%nz)
    end if
    gprf%xmin = minval(gprf%x)
    gprf%xmax = maxval(gprf%x)
    gprf%ymin = minval(gprf%y)
    gprf%ymax = maxval(gprf%y)
    gprf%zmin = minval(gprf%z)
    gprf%zmax = maxval(gprf%z)
    close(ug)
    contains
    subroutine check_dir_order(nx,ny,nz,kx,ky,kz,iflag)
      implicit none
      integer::nx,ny,nz,kx,ky,kz,iflag
      if(nx.lt.3) then
        iflag = 11
      else if (ny.lt.3) then 
        iflag = 12
      else if (nz.lt.3) then 
        iflag = 13
      else if (kx.lt.2 .or. kx.gt.6) then 
        iflag = 14
      else if (ky.lt.2 .or. ky.gt.6) then 
        iflag = 15
      else if (kz.lt.2 .or. kz.gt.6) then 
        iflag = 16
      end if
    end subroutine
  end subroutine read_gust_profile



  subroutine build_gust_interpolator(iflag)
    implicit none
    integer, intent(out)::iflag
    iflag = 0
    call gitp_u%initialize(gprf%x,gprf%y,gprf%z,gprf%u,gprf%kx,gprf%ky,gprf%kz,iflag)
    if (iflag.ne.0) then
      iflag = pack_dir_err(1, iflag)
      return
    end if
    call gitp_v%initialize(gprf%x,gprf%y,gprf%z,gprf%v,gprf%kx,gprf%ky,gprf%kz,iflag)
    if (iflag.ne.0) then
      iflag = pack_dir_err(2, iflag)
      return
    end if
    call gitp_w%initialize(gprf%x,gprf%y,gprf%z,gprf%w,gprf%kx,gprf%ky,gprf%kz,iflag)
    if (iflag.ne.0) then
      iflag = pack_dir_err(3, iflag)
      return
    end if
  end subroutine build_gust_interpolator


  

  subroutine gust_vel(xc, yc, zc, time, qiv, ug, vg, wg, iflag)
    real::xc, yc, zc, time, qiv(5), ug, vg, wg
    real(wp)::xg,yg,zg
    integer::iflag
    xg = xc - qiv(2)*time
    yg = yc - qiv(3)*time
    zg = zc - qiv(4)*time
    if ( &
    xg .gt. gprf%xmax .or. xg .lt. gprf%xmin .or. &
    yg .gt. gprf%ymax .or. yg .lt. gprf%ymin .or. &
    zg .gt. gprf%zmax .or. zg .lt. gprf%zmin ) then
      ug = 0.
      vg = 0.
      wg = 0.
    else
      call gitp_u%evaluate(xg,yg,zg,0,0,0,ug,iflag)
      if (iflag .ne. 0) then
        iflag = pack_dir_err(1, iflag)
        return
      endif 
      call gitp_v%evaluate(xg,yg,zg,0,0,0,vg,iflag)
      if (iflag .ne. 0) then
        iflag = pack_dir_err(2, iflag)
        return
      endif 
      call gitp_w%evaluate(xg,yg,zg,0,0,0,wg,iflag)
      if (iflag .ne. 0) then
        iflag = pack_dir_err(3, iflag)
        return
      endif
    end if
  end subroutine gust_vel


  function pack_dir_err(idir, iflag) result(flg)
    implicit none
    integer:: idir, iflag, flg
    flg = sign(idir,iflag)*10000+iflag
  end function pack_dir_err


  function get_err_msg(iflag) result(msg)
    integer:: iflag, idir, ifspline
    character(len=:),allocatable :: msg    !! status message associated with the flag
    character(len=10) :: istr   !! for integer to string conversion
    integer           :: istat  !! for write statement
    integer           :: ip3dgrd,ialph
    common /igrdtyp/ ip3dgrd,ialph
    idir = abs(iflag / 10000)
    ifspline = mod(iflag, 10000)
    select case (idir)
    case(0)
      select case (ifspline)
      case(0)
        msg = 'Succeeded'
      case(1)
        msg = 'Error in read_gust_profile: "gust_profile.inp" does not exist'
      case(2)
        msg = 'Error in read_gust_profile: unit numer has been used occupied'
      case(3)
        msg = 'Error in read_gust_profile: gust module does not compatible with ialph=1'
      case(11)
        msg = 'Error in read_gust_profile: at least 3 points in x-direction'
      case(12)
        msg = 'Error in read_gust_profile: at least 3 points in y-direction'
      case(13)
        msg = 'Error in read_gust_profile: at least 3 points in z-direction'
      case(14)
        msg = 'Error in read_gust_profile: invlaid kx, 2<=kx<=6'
      case(15)
        msg = 'Error in read_gust_profile: invlaid ky, 2<=ky<=6'
      case(16)
        msg = 'Error in read_gust_profile: invlaid kz, 2<=kz<=6'
      case default
        write(istr,fmt='(I10)',iostat=istat) iflag
        msg = 'Error in "gust_profile.dat": unknown status flag: '//trim(adjustl(istr))
      end select
    case(1)
      msg = 'Error in ug: ' // get_status_message(ifspline)
    case(2)
      if(ialph .ne. 0) then
        msg = 'Error in vg: ' // get_status_message(ifspline)
      else
        msg = 'Error in wg: ' // get_status_message(ifspline)
      end if
    case(3)
      if(ialph .ne. 0) then
        msg = 'Error in wg: ' // get_status_message(ifspline)
      else
        msg = 'Error in vg: ' // get_status_message(ifspline)
      end if
    case default
      write(istr,fmt='(I10)',iostat=istat) iflag
      msg = 'Unknown status flag: '//trim(adjustl(istr))
    end select
  end function get_err_msg

end module module_gust