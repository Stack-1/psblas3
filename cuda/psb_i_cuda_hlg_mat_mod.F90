!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
  

module psb_i_cuda_hlg_mat_mod

  use iso_c_binding
  use psb_i_mat_mod 
  use psb_i_hll_mat_mod


  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_i_hll_sparse_mat) :: psb_i_cuda_hlg_sparse_mat
    !
    ! ITPACK/HLL format, extended.
    ! We are adding here the routines to create a copy of the data
    ! into the GPU. 
    ! If HAVE_SPGPU is undefined this is just
    ! a copy of HLL, indistinguishable.
    ! 
    type(c_ptr) :: deviceMat = c_null_ptr
    integer     :: devstate  = is_host

  contains
    procedure, nopass  :: get_fmt       => i_cuda_hlg_get_fmt
    procedure, pass(a) :: sizeof        => i_cuda_hlg_sizeof
    procedure, pass(a) :: vect_mv       => psb_i_cuda_hlg_vect_mv
    procedure, pass(a) :: csmm          => psb_i_cuda_hlg_csmm
    procedure, pass(a) :: csmv          => psb_i_cuda_hlg_csmv
    procedure, pass(a) :: in_vect_sv    => psb_i_cuda_hlg_inner_vect_sv
    procedure, pass(a) :: scals         => psb_i_cuda_hlg_scals
    procedure, pass(a) :: scalv         => psb_i_cuda_hlg_scal
    procedure, pass(a) :: reallocate_nz => psb_i_cuda_hlg_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_i_cuda_hlg_allocate_mnnz
    ! Note: we do *not* need the TO methods, because the parent type
    ! methods will work. 
    procedure, pass(a) :: cp_from_coo   => psb_i_cuda_cp_hlg_from_coo
    procedure, pass(a) :: cp_from_fmt   => psb_i_cuda_cp_hlg_from_fmt
    procedure, pass(a) :: mv_from_coo   => psb_i_cuda_mv_hlg_from_coo
    procedure, pass(a) :: mv_from_fmt   => psb_i_cuda_mv_hlg_from_fmt
    procedure, pass(a) :: free          => i_cuda_hlg_free
    procedure, pass(a) :: mold          => psb_i_cuda_hlg_mold
    procedure, pass(a) :: is_host       => i_cuda_hlg_is_host
    procedure, pass(a) :: is_dev        => i_cuda_hlg_is_dev
    procedure, pass(a) :: is_sync       => i_cuda_hlg_is_sync
    procedure, pass(a) :: set_host      => i_cuda_hlg_set_host
    procedure, pass(a) :: set_dev       => i_cuda_hlg_set_dev
    procedure, pass(a) :: set_sync      => i_cuda_hlg_set_sync
    procedure, pass(a) :: sync          => i_cuda_hlg_sync
    procedure, pass(a) :: from_gpu      => psb_i_cuda_hlg_from_gpu
    procedure, pass(a) :: to_gpu        => psb_i_cuda_hlg_to_gpu
    final              :: i_cuda_hlg_finalize
  end type psb_i_cuda_hlg_sparse_mat

  private :: i_cuda_hlg_get_nzeros, i_cuda_hlg_free,  i_cuda_hlg_get_fmt, &
       & i_cuda_hlg_get_size, i_cuda_hlg_sizeof, i_cuda_hlg_get_nz_row


  interface 
    subroutine psb_i_cuda_hlg_vect_mv(alpha,a,x,beta,y,info,trans) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_, psb_i_base_vect_type, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(in)    :: a
      integer(psb_ipk_), intent(in)                 :: alpha, beta
      class(psb_i_base_vect_type), intent(inout) :: x
      class(psb_i_base_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)            :: trans
    end subroutine psb_i_cuda_hlg_vect_mv
  end interface

  interface 
    subroutine psb_i_cuda_hlg_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_i_cuda_hlg_sparse_mat, psb_ipk_,  psb_i_base_vect_type
      class(psb_i_cuda_hlg_sparse_mat), intent(in)    :: a
      integer(psb_ipk_), intent(in)                 :: alpha, beta
      class(psb_i_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)            :: trans
    end subroutine psb_i_cuda_hlg_inner_vect_sv
  end interface

  interface
    subroutine  psb_i_cuda_hlg_reallocate_nz(nz,a) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in)              :: nz
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
    end subroutine psb_i_cuda_hlg_reallocate_nz
  end interface

  interface
    subroutine  psb_i_cuda_hlg_allocate_mnnz(m,n,a,nz) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional    :: nz
    end subroutine psb_i_cuda_hlg_allocate_mnnz
  end interface

  interface 
    subroutine psb_i_cuda_hlg_mold(a,b,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_i_base_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(in)                  :: a
      class(psb_i_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_i_cuda_hlg_mold
  end interface

  interface 
    subroutine psb_i_cuda_hlg_from_gpu(a,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_hlg_from_gpu
  end interface

  interface 
    subroutine psb_i_cuda_hlg_to_gpu(a,info, nzrm) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: nzrm
    end subroutine psb_i_cuda_hlg_to_gpu
  end interface

  interface 
    subroutine psb_i_cuda_cp_hlg_from_coo(a,b,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_i_coo_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      class(psb_i_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_cp_hlg_from_coo
  end interface
  
  interface 
    subroutine psb_i_cuda_cp_hlg_from_fmt(a,b,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_i_base_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      class(psb_i_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_cp_hlg_from_fmt
  end interface
  
  interface 
    subroutine psb_i_cuda_mv_hlg_from_coo(a,b,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_i_coo_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      class(psb_i_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_mv_hlg_from_coo
  end interface
  

  interface 
    subroutine psb_i_cuda_mv_hlg_from_fmt(a,b,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_i_base_sparse_mat, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout)  :: a
      class(psb_i_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_i_cuda_mv_hlg_from_fmt
  end interface
  
  interface 
    subroutine psb_i_cuda_hlg_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)              :: alpha, beta, x(:)
      integer(psb_ipk_), intent(inout)           :: y(:)
      integer(psb_ipk_), intent(out)          :: info
      character, optional, intent(in)         :: trans
    end subroutine psb_i_cuda_hlg_csmv
  end interface
  interface 
    subroutine psb_i_cuda_hlg_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)              :: alpha, beta, x(:,:)
      integer(psb_ipk_), intent(inout)           :: y(:,:)
      integer(psb_ipk_), intent(out)          :: info
      character, optional, intent(in)         :: trans
    end subroutine psb_i_cuda_hlg_csmm
  end interface
  
  interface 
    subroutine psb_i_cuda_hlg_scal(d,a,info, side) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in)                 :: d(:)
      integer(psb_ipk_), intent(out)             :: info
      character, intent(in), optional            :: side
    end subroutine psb_i_cuda_hlg_scal
  end interface
  
  interface
    subroutine psb_i_cuda_hlg_scals(d,a,info) 
      import :: psb_i_cuda_hlg_sparse_mat, psb_ipk_, psb_ipk_
      class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in)                 :: d
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_hlg_scals
  end interface
  

contains 

  ! == ===================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ===================================

  
  function i_cuda_hlg_sizeof(a) result(res)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res

    
    if (a%is_dev()) call a%sync()
    res = 8 
    res = res + psb_sizeof_ip  * size(a%val)
    res = res + psb_sizeof_ip * size(a%irn)
    res = res + psb_sizeof_ip * size(a%idiag)
    res = res + psb_sizeof_ip * size(a%hkoffs)
    res = res + psb_sizeof_ip * size(a%ja)
    ! Should we account for the shadow data structure
    ! on the GPU device side? 
    ! res = 2*res
      
  end function i_cuda_hlg_sizeof

  function i_cuda_hlg_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'HLG'
  end function i_cuda_hlg_get_fmt
  


  ! == ===================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ===================================  

  subroutine  i_cuda_hlg_free(a) 
    use hlldev_mod
    implicit none 
    integer(psb_ipk_) :: info
    class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a

    if (c_associated(a%deviceMat)) &
         & call freeHllDevice(a%deviceMat)
    a%deviceMat = c_null_ptr
    call a%psb_i_hll_sparse_mat%free()
    
    return

  end subroutine i_cuda_hlg_free

  
  subroutine  i_cuda_hlg_sync(a) 
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), target, intent(in) :: a
    class(psb_i_cuda_hlg_sparse_mat), pointer :: tmpa
    integer(psb_ipk_) :: info

    tmpa => a
    if (tmpa%is_host()) then 
      call tmpa%to_gpu(info)
    else if (tmpa%is_dev()) then 
      call tmpa%from_gpu(info)
    end if
    call tmpa%set_sync()
    return

  end subroutine i_cuda_hlg_sync

  subroutine i_cuda_hlg_set_host(a)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
    
    a%devstate = is_host
  end subroutine i_cuda_hlg_set_host

  subroutine i_cuda_hlg_set_dev(a)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
    
    a%devstate = is_dev
  end subroutine i_cuda_hlg_set_dev

  subroutine i_cuda_hlg_set_sync(a)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a
    
    a%devstate = is_sync
  end subroutine i_cuda_hlg_set_sync

  function i_cuda_hlg_is_dev(a) result(res)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(in) :: a
    logical  :: res
  
    res = (a%devstate == is_dev)
  end function i_cuda_hlg_is_dev
  
  function i_cuda_hlg_is_host(a) result(res)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(in) :: a
    logical  :: res

    res = (a%devstate == is_host)
  end function i_cuda_hlg_is_host

  function i_cuda_hlg_is_sync(a) result(res)
    implicit none 
    class(psb_i_cuda_hlg_sparse_mat), intent(in) :: a
    logical  :: res

    res = (a%devstate == is_sync)
  end function i_cuda_hlg_is_sync


  subroutine  i_cuda_hlg_finalize(a) 
    use hlldev_mod
    implicit none 
    type(psb_i_cuda_hlg_sparse_mat), intent(inout) :: a

    if (c_associated(a%deviceMat)) &
         &  call freeHllDevice(a%deviceMat)
    a%deviceMat = c_null_ptr
    
    return
  end subroutine i_cuda_hlg_finalize

end module psb_i_cuda_hlg_mat_mod
