!TODO:
! 1. there is no need to copy data in `sub_j` & `sub_k`


module module_fds_blend
  use module_kinds, only: wp
  use,intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  type, public :: fds_sigma_t
    integer :: iblk, level, jdim, kdim, idim
    type(fds_sigma_t), pointer :: next
#ifdef CMPLX
    complex(wp), allocatable:: sigma(:,:,:), sig_i(:,:,:), sig_j(:,:,:), sig_k(:,:,:)
#else 
    real(wp)   , allocatable:: sigma(:,:,:), sig_i(:,:,:), sig_j(:,:,:), sig_k(:,:,:)
#endif
  contains
    procedure, pass, public  :: init => fds_sigma_init
    procedure, pass, public  :: update => fds_sigma_update
    procedure, pass, public  :: sub_i => fds_sigma_sub_i
    procedure, pass, public  :: sub_j => fds_sigma_sub_j
    procedure, pass, public  :: sub_k => fds_sigma_sub_k
    procedure, pass, public  :: dump =>fds_sigma_dump
    procedure, pass, private :: interpolate => fds_sigma_interpolate
    procedure, pass, private :: update_faces => fds_sigma_update_faces
    procedure, pass, private :: update_i => fds_sigma_update_i
    procedure, pass, private :: update_j => fds_sigma_update_j
    procedure, pass, private :: update_k => fds_sigma_update_k
    final :: fds_sigma_dtor
  end type fds_sigma_t
  type(fds_sigma_t), public, allocatable, target, save :: fds_sigma(:)
#ifdef CMPLX
    complex(wp), allocatable, public, save:: wk_sigma(:)
#else 
    real(wp)   , allocatable, public, save:: wk_sigma(:)
#endif
  integer, public, save :: nwk_sigma
contains
  
  subroutine init_fds_blend(nblock, levelg, mblk2nd, jdimg, kdimg, idimg)
  implicit none
  integer,intent(in):: nblock
  integer,intent(in):: levelg(nblock),mblk2nd(nblock)
  integer,intent(in):: jdimg(nblock),kdimg(nblock),idimg(nblock)
  common /fdsblend/ ifdsblend, cfdsblend, cfds_smin, cfds_smax
  integer :: ifdsblend
#ifdef CMPLX
  complex(wp) :: cfdsblend, cfds_smin, cfds_smax
#else 
  real(wp)    :: cfdsblend, cfds_smin, cfds_smax
#endif  
  ! local
  integer :: iblk


  if (.not.allocated(fds_sigma)) then 
    allocate(fds_sigma(nblock))
    nwk_sigma = 0
    do iblk = 1,nblock
      call fds_sigma(iblk)%init(iblk,&
                                jdimg(iblk),kdimg(iblk),idimg(iblk),&
                                levelg(iblk),mblk2nd(iblk))
    enddo
    allocate(wk_sigma(nwk_sigma))
  endif
  end subroutine init_fds_blend

  subroutine final_fds_sigma()
  implicit none
  if (allocated(fds_sigma)) deallocate(fds_sigma)
  if (allocated(wk_sigma)) deallocate(wk_sigma)
  end subroutine final_fds_sigma


!===============================================================================
! Methods for fds_sigma_t
!===============================================================================
  subroutine fds_sigma_init(this,iblk,jdim,kdim,idim,lvblk,host)
    implicit none
    ! args & return
    class(fds_sigma_t), intent(inout):: this
    integer,intent(in)  :: iblk,jdim,kdim,idim,lvblk,host
    ! common
    common /mydist2/ nnodes,myhost,myid,mycomm
    integer :: nnodes,myhost,myid,mycomm
    common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
    integer :: levt,kode,mode,ncyc,mtt,icyc,level,lglobal
    common /fdsblend/ ifdsblend, cfdsblend, cfds_smin, cfds_smax
    integer :: ifdsblend
#ifdef CMPLX
    complex(wp) :: cfdsblend, cfds_smin, cfds_smax
#else 
    real(wp)    :: cfdsblend, cfds_smin, cfds_smax
#endif  
    !update max # cells for working array
    nwk_sigma = max(nwk_sigma, jdim*kdim*idim)
    this%iblk   = iblk
    this%jdim   = jdim
    this%kdim   = kdim
    this%idim   = idim
    this%level  = lvblk
    if (host.eq.myid.and.ifdsblend.gt.1) then 
      allocate(this%sigma(jdim,kdim,idim))
      allocate(this%sig_i(jdim,kdim,idim))
      allocate(this%sig_j(jdim,kdim,idim))
      allocate(this%sig_k(jdim,kdim,idim))
      this%sigma  = cfds_smax
      this%sig_i  = cfds_smax
      this%sig_j  = cfds_smax
      this%sig_k  = cfds_smax
    endif
    if (lvblk.lt.lglobal) fds_sigma(iblk-1)%next => fds_sigma(iblk)
  end subroutine fds_sigma_init

  subroutine fds_sigma_dtor(this)
  implicit none
    type(fds_sigma_t):: this
    if (allocated(this%sigma)) deallocate(this%sigma)
    if (allocated(this%sig_i)) deallocate(this%sig_i)
    if (allocated(this%sig_j)) deallocate(this%sig_j)
    if (allocated(this%sig_k)) deallocate(this%sig_k)
  end subroutine fds_sigma_dtor

  recursive subroutine fds_sigma_interpolate(this)
    implicit none
    class(fds_sigma_t), intent(inout) :: this
    ! local
    integer :: i2d,id1,jd1,kd1,nc,i,j,k
    if (associated(this%next)) then
      id1     = this%idim-1
      jd1     = this%jdim-1
      kd1     = this%kdim-1
      if ((id1.gt.1.and.mod(id1,2).ne.0).or.mod(jd1,2).ne.0.or.mod(kd1,2).ne.0) then 
        write(error_unit, *) "[",__FILE__, __LINE__, "] ", &
        "Dimension Error, BLOCK #", this%iblk, " I=", this%idim, " J=", this%jdim, " K=", this%kdim
        stop 1
      endif
      i2d = 0
      if (id1.eq.1) i2d = 1
      nc = 2**(3-i2d)
      do i=1,id1,2
        do k=1,kd1,2
          do j=1,jd1,2
            this%next%sigma(j/2+1,k/2+1,i/2+1) = sum(this%sigma(j:j+1,k:k+1,i:i+1-i2d))/nc
          enddo
        enddo
      enddo
      call this%next%update_faces()
      call this%next%interpolate()
    endif
  end subroutine fds_sigma_interpolate

  subroutine fds_sigma_update(this, sigma)
    implicit none
    class(fds_sigma_t), intent(inout) :: this
#ifdef CMPLX
    complex(wp), intent(in) :: sigma(this%jdim,this%kdim,this%idim)
#else 
    real(wp)   , intent(in) :: sigma(this%jdim,this%kdim,this%idim)
#endif  
    this%sigma = sigma
    call this%update_faces()
    call this%interpolate()
  end subroutine fds_sigma_update

  subroutine fds_sigma_dump(this, nt, x, y, z)
    implicit none
    class(fds_sigma_t), intent(in) :: this
    integer, intent(in) :: nt
#ifdef CMPLX
    complex(wp), dimension(this%jdim,this%kdim,this%idim) :: x,y,z
#else
    real(wp)   , dimension(this%jdim,this%kdim,this%idim) :: x,y,z
#endif
    !local
    integer :: uout, i, j, k
    real    :: tname
    character(len=50) :: name, stmp
    write(name, '(I0.5,''.'',I0.5,''.sigma'')') nt, this%iblk
    open(newunit=uout, file=trim(name), form='unformatted')
    write(uout) this%idim, this%jdim, this%kdim
    write(uout) this%level, nt
    write(uout) (((x(j,k,i), i=1,this%idim),j=1,this%jdim),k=1,this%kdim)
    write(uout) (((y(j,k,i), i=1,this%idim),j=1,this%jdim),k=1,this%kdim)
    write(uout) (((z(j,k,i), i=1,this%idim),j=1,this%jdim),k=1,this%kdim)
    write(uout) (((this%sigma(j,k,i), i=1,this%idim-1),j=1,this%jdim-1),k=1,this%kdim-1)
    close(uout)
  end subroutine fds_sigma_dump

  subroutine fds_sigma_update_faces(this)
    implicit none
    class(fds_sigma_t), intent(inout) :: this
    call this%update_i()
    call this%update_j()
    call this%update_k()
  end subroutine fds_sigma_update_faces

  subroutine fds_sigma_update_i(this)
    implicit none
    class(fds_sigma_t), intent(inout) :: this
    integer :: j,k,i,id0,id1,id2
#ifdef CMPLX
    complex(wp):: ccmax
#else
    real(wp)   :: ccmax
#endif
    id0 = this%idim
    id1 = id0-1
    id2 = id1-1
    if (id2.le.0) then
      this%sig_i = this%sigma
    else
      do i = 2,id1
        do k = 1,this%kdim
          do j = 1,this%jdim
            this%sig_i(j,k,i) = ccmax(this%sigma(j,k,i), this%sigma(j,k,i-1))
          enddo
        enddo
      enddo
      this%sig_i(:,:,  1) = this%sigma(:,:,  1)
      this%sig_i(:,:,id0) = this%sigma(:,:,id1)
    endif
  end subroutine fds_sigma_update_i

  subroutine fds_sigma_update_j(this)
    class(fds_sigma_t), intent(inout) :: this
    integer :: j,k,i,jd0,jd1,jd2
#ifdef CMPLX
    complex(wp):: ccmax
#else
    real(wp)   :: ccmax
#endif
    jd0 = this%jdim
    jd1 = jd0-1
    jd2 = jd1-1
    if(jd2.le.0) then
      this%sig_j = this%sigma
    else
      do i = 1,this%idim
        do k = 1,this%kdim
          do j = 2,jd0
            this%sig_j(j,k,i) = ccmax(this%sigma(j,k,i), this%sigma(j-1,k,i))
          enddo
        enddo
      enddo
      this%sig_j(  1,:,:) = this%sigma(  1,:,:)
      this%sig_j(jd0,:,:) = this%sigma(jd1,:,:)
    endif
  end subroutine fds_sigma_update_j

  subroutine fds_sigma_update_k(this)
    class(fds_sigma_t), intent(inout) :: this
    integer :: j,k,i,kd0,kd1,kd2
#ifdef CMPLX
    complex(wp):: ccmax
#else
    real(wp)   :: ccmax
#endif
    kd0 = this%kdim
    kd1 = kd0-1
    kd2 = kd1-1
    if (kd2.le.0) then
      this%sig_k = this%sigma
    else
      do i = 1,this%idim
        do k = 2,kd0
          do j = 1,this%jdim
            this%sig_k(j,k,i) = ccmax(this%sigma(j,k,i), this%sigma(j,k-1,i))
          enddo
        enddo
      enddo
      this%sig_k(:,  1,:)= this%sigma(:,  1,:)
      this%sig_k(:,kd0,:)= this%sigma(:,kd1,:)
    endif
  end subroutine fds_sigma_update_k

  subroutine fds_sigma_sub_i(this, k, nk, sig)
    class(fds_sigma_t), intent(in) :: this
    integer, intent(in) :: k, nk
#ifdef CMPLX
    complex(wp), intent(out) :: sig(this%jdim, nk, this%idim)
#else 
    real(wp)   , intent(out) :: sig(this%jdim, nk, this%idim)
#endif
    sig = this%sig_i(:,k:k+nk-1,:)
  end subroutine

  subroutine fds_sigma_sub_j(this, i, ni, sig)
    class(fds_sigma_t), intent(in) :: this
    integer, intent(in) :: i, ni
#ifdef CMPLX
    complex(wp), intent(out) :: sig(this%jdim, this%kdim, ni)
#else 
    real(wp)   , intent(out) :: sig(this%jdim, this%kdim, ni)
#endif
    sig = this%sig_j(:,:,i:i+ni-1)
  end subroutine fds_sigma_sub_j
  
  subroutine fds_sigma_sub_k(this, i, ni, sig)
    class(fds_sigma_t), intent(in) :: this
    integer, intent(in) :: i, ni
#ifdef CMPLX
    complex(wp), intent(out) :: sig(this%jdim, this%kdim, ni)
#else 
    real(wp)   , intent(out) :: sig(this%jdim, this%kdim, ni)
#endif
    sig = this%sig_k(:,:,i:i+ni-1)
  end subroutine fds_sigma_sub_k
  
end module module_fds_blend