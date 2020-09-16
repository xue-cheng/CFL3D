module module_des_save


use module_kinds, only: wp
use,intrinsic :: iso_fortran_env, only: error_unit
implicit none

#ifdef CMPLX
#  define DESV complex(wp)
#else
#  define DESV real(wp)
#endif

  type, public :: des_geom_t
    DESV, allocatable :: hmax(:,:,:)   !jd1,kd1,id1
    DESV, allocatable :: hmin(:,:,:)   !jd1,kd1,id1
    DESV, allocatable :: rn(:,:,:,:,:) !3,8,jd1,kd1,id1
  contains
    procedure, public, pass :: update_h => des_geom_upd_h
    procedure, public, pass :: update_r => des_geom_upd_r
  end type des_geom_t
  DESV,public,save :: des_small
  type(des_geom_t),allocatable,public,save::desgeo(:)
contains

  subroutine init_des_info(nblock)
    implicit none
    integer,intent(in):: nblock
    !common
    common /zero/ iexp
    integer :: iexp
    des_small = 10.**(-iexp+1)
    ! just allocate here, use lazy init
    if (.not.allocated(desgeo)) allocate(desgeo(nblock))
  end subroutine init_des_info

  subroutine des_geom_upd_h(this,jdim,kdim,idim,sj,sk,si,vol,is2d,isunst)
  implicit none
  class(des_geom_t), intent(inout):: this
  integer, intent(in) :: jdim,kdim,idim
  logical, intent(in) :: isunst, is2d
  DESV :: sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5)
  DESV :: vol(jdim,kdim,idim-1)
  !local
  integer :: j,k,i,jd1,kd1,id1
  DESV :: hj, hk, hi,hmax,hmin
  logical :: need_update
  ! function 
  DESV ::ccmax,ccmin
  need_update = isunst
  if (.not.allocated(this%hmax)) then 
    need_update=.true.
    allocate(this%hmax(jdim-1,kdim-1,idim-1))
  endif
  if (.not.allocated(this%hmin)) then 
    need_update=.true.
    allocate(this%hmin(jdim-1,kdim-1,idim-1))
  endif
  if (need_update) then
    jd1 = jdim-1
    kd1 = kdim-1
    id1 = idim-1
    do i=1,idim-1
      do j=1,jdim-1
        do k=1,kdim-1
          hj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
          hk = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
          hmax = ccmax(hj, hk)
          hmin = ccmin(hj, hk)
          if(.not.is2d) then 
            hi = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
            hmax = ccmax(hmax, hi)
            hmin = ccmin(hmin, hi)
          endif
          this%hmax(j,k,i) = hmax
          this%hmin(j,k,i) = hmin
        enddo
      enddo
    enddo
  endif
  end subroutine des_geom_upd_h

  
  subroutine des_geom_upd_r(this,jdim,kdim,idim,x,y,z,isunst)
    implicit none
    class(des_geom_t), intent(inout):: this
    integer, intent(in) :: jdim,kdim,idim
    logical, intent(in) :: isunst
    DESV :: x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
    !local
    integer :: j,k,i,jd1,kd1,id1
    logical :: need_update
    DESV :: cx, cy, cz
    ! function 
    need_update = isunst
    if (.not.allocated(this%rn)) then 
      need_update=.true.
      allocate(this%rn(3,8,jdim-1,kdim-1,idim-1))
    endif
    if (need_update) then
      jd1 = jdim-1
      kd1 = kdim-1
      id1 = idim-1
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            cx=.125*(x(j,k  ,i  )+x(j+1,k  ,i  )&
                    +x(j,k+1,i  )+x(j+1,k+1,i  )&
                    +x(j,k  ,i+1)+x(j+1,k  ,i+1)&
                    +x(j,k+1,i+1)+x(j+1,k+1,i+1))
            cy=.125*(y(j,k  ,i  )+y(j+1,k  ,i  )&
                    +y(j,k+1,i  )+y(j+1,k+1,i  )&
                    +y(j,k  ,i+1)+y(j+1,k  ,i+1)&
                    +y(j,k+1,i+1)+y(j+1,k+1,i+1))
            cz=.125*(z(j,k  ,i  )+z(j+1,k  ,i  )&
                    +z(j,k+1,i  )+z(j+1,k+1,i  )&
                    +z(j,k  ,i+1)+z(j+1,k  ,i+1)&
                    +z(j,k+1,i+1)+z(j+1,k+1,i+1))
            this%rn(1,1,j,k,i) = x(j  ,k  ,i  )-cx
            this%rn(1,2,j,k,i) = x(j+1,k  ,i  )-cx
            this%rn(1,3,j,k,i) = x(j  ,k+1,i  )-cx
            this%rn(1,4,j,k,i) = x(j+1,k+1,i  )-cx
            this%rn(1,5,j,k,i) = x(j  ,k  ,i+1)-cx
            this%rn(1,6,j,k,i) = x(j+1,k  ,i+1)-cx
            this%rn(1,7,j,k,i) = x(j  ,k+1,i+1)-cx
            this%rn(1,8,j,k,i) = x(j+1,k+1,i+1)-cx
            this%rn(2,1,j,k,i) = y(j  ,k  ,i  )-cy
            this%rn(2,2,j,k,i) = y(j+1,k  ,i  )-cy
            this%rn(2,3,j,k,i) = y(j  ,k+1,i  )-cy
            this%rn(2,4,j,k,i) = y(j+1,k+1,i  )-cy
            this%rn(2,5,j,k,i) = y(j  ,k  ,i+1)-cy
            this%rn(2,6,j,k,i) = y(j+1,k  ,i+1)-cy
            this%rn(2,7,j,k,i) = y(j  ,k+1,i+1)-cy
            this%rn(2,8,j,k,i) = y(j+1,k+1,i+1)-cy
            this%rn(3,1,j,k,i) = z(j  ,k  ,i  )-cz
            this%rn(3,2,j,k,i) = z(j+1,k  ,i  )-cz
            this%rn(3,3,j,k,i) = z(j  ,k+1,i  )-cz
            this%rn(3,4,j,k,i) = z(j+1,k+1,i  )-cz
            this%rn(3,5,j,k,i) = z(j  ,k  ,i+1)-cz
            this%rn(3,6,j,k,i) = z(j+1,k  ,i+1)-cz
            this%rn(3,7,j,k,i) = z(j  ,k+1,i+1)-cz
            this%rn(3,8,j,k,i) = z(j+1,k+1,i+1)-cz
          enddo
        enddo
      enddo
    endif
  end subroutine des_geom_upd_r
  
end module module_des_save