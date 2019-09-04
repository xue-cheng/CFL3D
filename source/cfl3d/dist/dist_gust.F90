module dist_gust
  use module_gust
  use mpi
  use bspline_kinds_module, only: wp

  implicit none
  contains
  subroutine distribute_gust_profile(iflag)
    implicit none
    integer :: iflag
    integer :: nnodes,myhost,myid,mycomm
    integer :: nnsize
    integer :: MY_REAL_PRECISION
    common /mydist2/ nnodes,myhost,myid,mycomm
#   ifdef DBLE_PRECSN
  MY_REAL_PRECISION = MPI_DOUBLE_PRECISION
#   else
  MY_REAL_PRECISION = MPI_REAL
#   endif
    call MPI_Bcast(gprf%nx, 1, MPI_INTEGER, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%ny, 1, MPI_INTEGER, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%nz, 1, MPI_INTEGER, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%kx, 1, MPI_INTEGER, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%ky, 1, MPI_INTEGER, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%kz, 1, MPI_INTEGER, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%xmin, 1, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%xmax, 1, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%ymin, 1, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%ymax, 1, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%zmin, 1, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%zmax, 1, MY_REAL_PRECISION, myhost, mycomm, iflag)
    if ( myid .ne. myhost ) then
      allocate(gprf%x(gprf%nx), gprf%y(gprf%ny), gprf%z(gprf%nz))
      allocate(gprf%u(gprf%nx,gprf%ny,gprf%nz))
      allocate(gprf%v(gprf%nx,gprf%ny,gprf%nz))
      allocate(gprf%w(gprf%nx,gprf%ny,gprf%nz))
    end if
    nnsize = gprf%nx*gprf%ny*gprf%nz
    call MPI_Bcast(gprf%x, gprf%nx, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%y, gprf%ny, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%z, gprf%nz, MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%u, nnsize , MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%v, nnsize , MY_REAL_PRECISION, myhost, mycomm, iflag)
    call MPI_Bcast(gprf%w, nnsize , MY_REAL_PRECISION, myhost, mycomm, iflag)
    end subroutine distribute_gust_profile
end module dist_gust