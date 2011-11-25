module psb_d_base_vect_mod
  
  use psb_const_mod
  use psb_error_mod

  type psb_d_base_vect_type
    real(psb_dpk_), allocatable :: v(:)
  contains
    procedure, pass(x) :: get_nrows => d_base_get_nrows
    procedure, pass(x) :: sizeof   => d_base_sizeof
    procedure, pass(x) :: dot_v    => d_base_dot_v
    procedure, pass(x) :: dot_a    => d_base_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => d_base_axpby_v
    procedure, pass(y) :: axpby_a  => d_base_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    procedure, pass(y) :: mlt_v    => d_base_mlt_v
    procedure, pass(y) :: mlt_a    => d_base_mlt_a
    procedure, pass(z) :: mlt_a_2  => d_base_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => d_base_mlt_v_2
    procedure, pass(z) :: mlt_va   => d_base_mlt_va
    procedure, pass(z) :: mlt_av   => d_base_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2, mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: scal     => d_base_scal
    procedure, pass(x) :: nrm2     => d_base_nrm2
    procedure, pass(x) :: amax     => d_base_amax
    procedure, pass(x) :: asum     => d_base_asum
    procedure, pass(x) :: all      => d_base_all
    procedure, pass(x) :: zero     => d_base_zero
    procedure, pass(x) :: asb      => d_base_asb
    procedure, pass(x) :: sync     => d_base_sync
    procedure, pass(x) :: gthab    => d_base_gthab
    procedure, pass(x) :: gthzv    => d_base_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => d_base_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => d_base_free
    procedure, pass(x) :: ins      => d_base_ins
    procedure, pass(x) :: bld_x    => d_base_bld_x
    procedure, pass(x) :: bld_n    => d_base_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: getCopy  => d_base_getCopy
    procedure, pass(x) :: cpy_vect => d_base_cpy_vect
    generic, public    :: assignment(=) => cpy_vect, set_scal
    procedure, pass(x) :: set_scal => d_base_set_scal
    procedure, pass(x) :: set_vect => d_base_set_vect
    generic, public    :: set      => set_vect, set_scal
  end type psb_d_base_vect_type

  public  :: psb_d_base_vect
  private :: constructor, size_const
  interface psb_d_base_vect
    module procedure constructor, size_const
  end interface psb_d_base_vect

contains
  
  subroutine d_base_bld_x(x,this)
    use psb_realloc_mod
    real(psb_dpk_), intent(in) :: this(:)
    class(psb_d_base_vect_type), intent(inout) :: x
    integer :: info

    call psb_realloc(size(this),x%v,info)
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,'base_vect_bld')
      return
    end if
    x%v(:)  = this(:)

  end subroutine d_base_bld_x
    
  
  subroutine d_base_bld_n(x,n)
    use psb_realloc_mod
    integer, intent(in) :: n
    class(psb_d_base_vect_type), intent(inout) :: x
    integer :: info
    
    call psb_realloc(n,x%v,info)
    call x%asb(n,info)

  end subroutine d_base_bld_n
    
  function  d_base_getCopy(x) result(res)
    class(psb_d_base_vect_type), intent(in)  :: x
    real(psb_dpk_), allocatable              :: res(:)
    integer :: info
    
    allocate(res(x%get_nrows()),stat=info) 
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,'base_getCopy')
      return
    end if
    res(:) = x%v(:)
  end function d_base_getCopy
    
  subroutine d_base_cpy_vect(res,x)
    real(psb_dpk_), allocatable, intent(out) :: res(:)
    class(psb_d_base_vect_type), intent(in)  :: x
    integer :: info

    res = x%v 

  end subroutine d_base_cpy_vect

  subroutine d_base_set_scal(x,val)
    class(psb_d_base_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val
        
    integer :: info
    x%v = val
    
  end subroutine d_base_set_scal

  subroutine d_base_set_vect(x,val)
    class(psb_d_base_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val(:)
    integer :: nr
    integer :: info
    if (allocated(x%v)) then 
      nr = min(size(x%v),size(val))
      x%v(1:nr) = val(1:nr)
    else
      x%v = val
    end if
  end subroutine d_base_set_vect
    
  
  function constructor(x) result(this)
    real(psb_dpk_)   :: x(:)
    type(psb_d_base_vect_type) :: this
    integer :: info

    this%v = x
    call this%asb(size(x),info)
  end function constructor
    
  
  function size_const(n) result(this)
    integer, intent(in) :: n
    type(psb_d_base_vect_type) :: this
    integer :: info

    call this%asb(n,info)

  end function size_const
    
  function d_base_get_nrows(x) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(in) :: x
    integer :: res
    res = 0
    if (allocated(x%v)) res = size(x%v)
  end function d_base_get_nrows

  function d_base_sizeof(x) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(in) :: x
    integer(psb_long_int_k_) :: res
    res = psb_sizeof_dp*x%get_nrows()
  end function d_base_sizeof

  function d_base_dot_v(n,x,y) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x, y
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: ddot
    
    res = dzero
    !
    ! Note: this is the base implementation.
    !  When we get here, we are sure that X is of
    !  TYPE psb_d_base_vect
    !
    select type(yy => y)
    type is (psb_d_base_vect_type)
      res = ddot(n,x%v,1,y%v,1)
    class default
      res = y%dot(n,x%v)
    end select

  end function d_base_dot_v

  function d_base_dot_a(n,x,y) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x
    real(psb_dpk_), intent(in)    :: y(:)
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: ddot
    
    res = ddot(n,y,1,x%v,1)

  end function d_base_dot_a
    
  subroutine d_base_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer, intent(in)               :: m
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer, intent(out)              :: info
    
    select type(xx => x)
    type is (psb_d_base_vect_type)
      call psb_geaxpby(m,alpha,x%v,beta,y%v,info)
    class default
      call y%axpby(m,alpha,x%v,beta,info)
    end select

  end subroutine d_base_axpby_v

  subroutine d_base_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer, intent(in)               :: m
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer, intent(out)              :: info
    
    call psb_geaxpby(m,alpha,x,beta,y%v,info)
    
  end subroutine d_base_axpby_a

    
  subroutine d_base_mlt_v(x, y, info, xconj)
    use psi_serial_mod
    use psb_string_mod
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x
    class(psb_d_base_vect_type), intent(inout) :: y
    integer, intent(out)                       :: info    
    character, intent(in), optional            :: xconj
    integer   :: i, n
    character :: xconj_

    info = 0
    if (present(xconj)) then 
      xconj_ = (psb_toupper(xconj))
    else
      xconj_ = 'N'
    end if

    select type(xx => x)
    type is (psb_d_base_vect_type)
      n = min(size(y%v), size(xx%v))
      do i=1, n 
        y%v(i) = y%v(i)*xx%v(i)
      end do
    class default
      call y%mlt(x%v,info)
    end select

  end subroutine d_base_mlt_v

  subroutine d_base_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: y
    integer, intent(out)              :: info
    integer :: i, n

    info = 0
    n = min(size(y%v), size(x))
    do i=1, n 
      y%v(i) = y%v(i)*x(i)
    end do
    
  end subroutine d_base_mlt_a


  subroutine d_base_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:)
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info
    integer :: i, n

    info = 0    
    n = min(size(z%v), size(x), size(y))
!!$    write(0,*) 'Mlt_a_2: ',n
    if (alpha == dzero) then 
      if (beta == done) then 
        return 
      else
        do i=1, n
          z%v(i) = beta*z%v(i)
        end do
      end if
    else
      if (alpha == done) then 
        if (beta == dzero) then 
          do i=1, n 
            z%v(i) = y(i)*x(i)
          end do
        else if (beta == done) then 
          do i=1, n 
            z%v(i) = z%v(i) + y(i)*x(i)
          end do
        else 
          do i=1, n 
            z%v(i) = beta*z%v(i) + y(i)*x(i)
          end do
        end if
      else if (alpha == -done) then 
        if (beta == dzero) then 
          do i=1, n 
            z%v(i) = -y(i)*x(i)
          end do
        else if (beta == done) then 
          do i=1, n 
            z%v(i) = z%v(i) - y(i)*x(i)
          end do
        else 
          do i=1, n 
            z%v(i) = beta*z%v(i) - y(i)*x(i)
          end do
        end if
      else
        if (beta == dzero) then 
          do i=1, n 
            z%v(i) = alpha*y(i)*x(i)
          end do
        else if (beta == done) then 
          do i=1, n 
            z%v(i) = z%v(i) + alpha*y(i)*x(i)
          end do
        else 
          do i=1, n 
            z%v(i) = beta*z%v(i) + alpha*y(i)*x(i)
          end do
        end if
      end if
    end if
  end subroutine d_base_mlt_a_2

  subroutine d_base_mlt_v_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: y
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    
    call z%mlt(alpha,x%v,y%v,beta,info)

  end subroutine d_base_mlt_v_2

  subroutine d_base_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: y
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    
    call z%mlt(alpha,x,y%v,beta,info)

  end subroutine d_base_mlt_av

  subroutine d_base_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:)
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    
    call z%mlt(alpha,y,x,beta,info)

  end subroutine d_base_mlt_va

  subroutine d_base_scal(alpha, x)
    use psi_serial_mod
    implicit none 
    class(psb_d_base_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent (in)       :: alpha
    
    if (allocated(x%v)) x%v = alpha*x%v

  end subroutine d_base_scal


  function d_base_nrm2(n,x) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: dnrm2
    
    res =  dnrm2(n,x%v,1)

  end function d_base_nrm2
  
  function d_base_amax(n,x) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    
    res =  maxval(abs(x%v(1:n)))

  end function d_base_amax

  function d_base_asum(n,x) result(res)
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    
    res =  sum(abs(x%v(1:n)))

  end function d_base_asum
  
  subroutine d_base_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in)               :: n
    class(psb_d_base_vect_type), intent(out)    :: x
    integer, intent(out)              :: info
    
    call psb_realloc(n,x%v,info)
    
  end subroutine d_base_all

  subroutine d_base_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_d_base_vect_type), intent(inout)    :: x
    
    if (allocated(x%v)) x%v=dzero

  end subroutine d_base_zero

  subroutine d_base_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in)              :: n
    class(psb_d_base_vect_type), intent(inout) :: x
    integer, intent(out)             :: info
    
    if (x%get_nrows() < n) &
         & call psb_realloc(n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')

  end subroutine d_base_asb

  subroutine d_base_sync(x)
    implicit none 
    class(psb_d_base_vect_type), intent(inout) :: x
    
    !
    ! The base version does nothing, it's just
    ! a placeholder.
    ! 
    
  end subroutine d_base_sync

  subroutine d_base_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer :: n, idx(:)
    real(psb_dpk_) :: alpha, beta, y(:)
    class(psb_d_base_vect_type) :: x
    
    call x%sync()
    call psi_gth(n,idx,alpha,x%v,beta,y)

  end subroutine d_base_gthab

  subroutine d_base_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer :: n, idx(:)
    real(psb_dpk_) ::  y(:)
    class(psb_d_base_vect_type) :: x
    
    call x%sync()
    call psi_gth(n,idx,x%v,y)

  end subroutine d_base_gthzv

  subroutine d_base_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer :: n, idx(:)
    real(psb_dpk_) :: beta, x(:)
    class(psb_d_base_vect_type) :: y
    
    call y%sync()
    call psi_sct(n,idx,x,beta,y%v)

  end subroutine d_base_sctb

  subroutine d_base_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_d_base_vect_type), intent(inout)  :: x
    integer, intent(out)              :: info
    
    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (info /= 0) call & 
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')
        
  end subroutine d_base_free

  subroutine d_base_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_d_base_vect_type), intent(inout)  :: x
    integer, intent(in)               :: n, dupl
    integer, intent(in)               :: irl(:)
    real(psb_dpk_), intent(in)        :: val(:)
    integer, intent(out)              :: info

    integer :: i

    info = 0
    if (psb_errstatus_fatal()) return 

    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
    else if (n > min(size(irl),size(val))) then 
      info = psb_err_invalid_input_

    else 
      select case(dupl) 
      case(psb_dupl_ovwrt_) 
        do i = 1, n
          !loop over all val's rows

          ! row actual block row 
          if (irl(i) > 0) then
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(irl(i)) = val(i)
          end if
        enddo

      case(psb_dupl_add_) 

        do i = 1, n
          !loop over all val's rows

          if (irl(i) > 0) then
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(irl(i)) = x%v(irl(i)) +  val(i)
          end if
        enddo

      case default
        info = 321
!!$      call psb_errpush(info,name)
!!$      goto 9999
      end select
    end if
    if (info /= 0) then 
      call psb_errpush(info,'base_vect_ins')
      return
    end if

  end subroutine d_base_ins

end module psb_d_base_vect_mod