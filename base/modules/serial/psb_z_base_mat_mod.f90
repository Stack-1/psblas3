!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
!    
!
!
module psb_z_base_mat_mod
  
  use psb_base_mat_mod
  use psb_z_base_vect_mod


  !> \namespace  psb_base_mod  \class  psb_z_base_sparse_mat
  !! \extends psb_base_mat_mod::psb_base_sparse_mat
  !! The psb_z_base_sparse_mat type, extending psb_base_sparse_mat,
  !! defines a middle level  complex(psb_dpk_) sparse matrix object.
  !! This class object itself does not have any additional members
  !! with respect to those of the base class. Most methods cannot be fully
  !! implemented at this level, but we can define the interface for the
  !! computational methods requiring the knowledge of the underlying
  !! field, such as the matrix-vector product; this interface is defined,
  !! but is supposed to be overridden at the leaf level.
  !!
  !! About the method MOLD: this has been defined for those compilers
  !! not yet supporting ALLOCATE( ...,MOLD=...); it's otherwise silly to
  !! duplicate "by hand" what is specified in the language (in this case F2008)
  !!
  type, extends(psb_base_sparse_mat) :: psb_z_base_sparse_mat
  contains
    !
    ! Data management methods: defined here, but (mostly) not implemented.
    !    
    procedure, pass(a) :: csput_a       => psb_z_base_csput_a
    procedure, pass(a) :: csput_v       => psb_z_base_csput_v 
    generic, public    :: csput         => csput_a,  csput_v
    procedure, pass(a) :: csgetrow      => psb_z_base_csgetrow
    procedure, pass(a) :: csgetblk      => psb_z_base_csgetblk
    procedure, pass(a) :: get_diag      => psb_z_base_get_diag
    generic, public    :: csget         => csgetrow, csgetblk 
    procedure, pass(a) :: tril          => psb_z_base_tril
    procedure, pass(a) :: triu          => psb_z_base_triu
    procedure, pass(a) :: csclip        => psb_z_base_csclip 
    procedure, pass(a) :: cp_to_coo     => psb_z_base_cp_to_coo   
    procedure, pass(a) :: cp_from_coo   => psb_z_base_cp_from_coo 
    procedure, pass(a) :: cp_to_fmt     => psb_z_base_cp_to_fmt   
    procedure, pass(a) :: cp_from_fmt   => psb_z_base_cp_from_fmt 
    procedure, pass(a) :: mv_to_coo     => psb_z_base_mv_to_coo   
    procedure, pass(a) :: mv_from_coo   => psb_z_base_mv_from_coo 
    procedure, pass(a) :: mv_to_fmt     => psb_z_base_mv_to_fmt   
    procedure, pass(a) :: mv_from_fmt   => psb_z_base_mv_from_fmt 
    procedure, pass(a) :: mold          => psb_z_base_mold 
    procedure, pass(a) :: clone         => psb_z_base_clone
    procedure, pass(a) :: make_nonunit  => psb_z_base_make_nonunit
    procedure, pass(a) :: clean_zeros   => psb_z_base_clean_zeros
    !
    ! Convert internal indices
    !
    procedure, pass(a) :: cp_to_lcoo     => psb_z_base_cp_to_lcoo   
    procedure, pass(a) :: cp_from_lcoo   => psb_z_base_cp_from_lcoo 
    procedure, pass(a) :: cp_to_lfmt     => psb_z_base_cp_to_lfmt   
    procedure, pass(a) :: cp_from_lfmt   => psb_z_base_cp_from_lfmt 
    procedure, pass(a) :: mv_to_lcoo     => psb_z_base_mv_to_lcoo   
    procedure, pass(a) :: mv_from_lcoo   => psb_z_base_mv_from_lcoo 
    procedure, pass(a) :: mv_to_lfmt     => psb_z_base_mv_to_lfmt   
    procedure, pass(a) :: mv_from_lfmt   => psb_z_base_mv_from_lfmt 

    
    !
    ! Transpose methods: defined here but not implemented. 
    !    
    procedure, pass(a) :: transp_1mat => psb_z_base_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_z_base_transp_2mat
    procedure, pass(a) :: transc_1mat => psb_z_base_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_z_base_transc_2mat
    
    !
    ! Computational methods: defined here but not implemented. 
    !    
    procedure, pass(a) :: vect_mv     => psb_z_base_vect_mv
    procedure, pass(a) :: csmv        => psb_z_base_csmv
    procedure, pass(a) :: csmm        => psb_z_base_csmm
    generic, public    :: spmm        => csmm, csmv, vect_mv
    procedure, pass(a) :: in_vect_sv  => psb_z_base_inner_vect_sv
    procedure, pass(a) :: inner_cssv  => psb_z_base_inner_cssv    
    procedure, pass(a) :: inner_cssm  => psb_z_base_inner_cssm
    generic, public    :: inner_spsm  => inner_cssm, inner_cssv, in_vect_sv
    procedure, pass(a) :: vect_cssv   => psb_z_base_vect_cssv
    procedure, pass(a) :: cssv        => psb_z_base_cssv
    procedure, pass(a) :: cssm        => psb_z_base_cssm
    generic, public    :: spsm        => cssm, cssv, vect_cssv
    procedure, pass(a) :: scals       => psb_z_base_scals
    procedure, pass(a) :: scalv       => psb_z_base_scal
    generic, public    :: scal        => scals, scalv
    procedure, pass(a) :: maxval      => psb_z_base_maxval
    procedure, pass(a) :: spnmi       => psb_z_base_csnmi
    procedure, pass(a) :: spnm1       => psb_z_base_csnm1
    procedure, pass(a) :: rowsum      => psb_z_base_rowsum
    procedure, pass(a) :: arwsum      => psb_z_base_arwsum
    procedure, pass(a) :: colsum      => psb_z_base_colsum
    procedure, pass(a) :: aclsum      => psb_z_base_aclsum
  end type psb_z_base_sparse_mat
  
  private :: z_base_mat_sync, z_base_mat_is_host, z_base_mat_is_dev, &
       & z_base_mat_is_sync, z_base_mat_set_host, z_base_mat_set_dev,&
       & z_base_mat_set_sync
  
  !> \namespace  psb_base_mod  \class  psb_z_coo_sparse_mat
  !! \extends psb_z_base_mat_mod::psb_z_base_sparse_mat
  !! 
  !! psb_z_coo_sparse_mat type and the related methods. This is the
  !! reference type for all the format transitions, copies and mv unless
  !! methods are implemented that allow the direct transition from one
  !! format to another. It is defined here since all other classes must
  !! refer to it per the MEDIATOR design pattern.
  !!
  type, extends(psb_z_base_sparse_mat) :: psb_z_coo_sparse_mat
    !> Number of nonzeros.
    integer(psb_ipk_) :: nnz
    !> Row indices.
    integer(psb_ipk_), allocatable :: ia(:)
    !> Column indices.
    integer(psb_ipk_), allocatable :: ja(:)
    !> Coefficient values. 
    complex(psb_dpk_), allocatable :: val(:)

    integer, private   :: sort_status=psb_unsorted_
    
  contains
    !
    ! Data management methods. 
    !    
    procedure, pass(a) :: get_size     => z_coo_get_size
    procedure, pass(a) :: get_nzeros   => z_coo_get_nzeros
    procedure, nopass  :: get_fmt      => z_coo_get_fmt
    procedure, pass(a) :: sizeof       => z_coo_sizeof
    procedure, pass(a) :: reallocate_nz => psb_z_coo_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_z_coo_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_z_cp_coo_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_z_cp_coo_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_z_cp_coo_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_z_cp_coo_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_z_mv_coo_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_z_mv_coo_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_z_mv_coo_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_z_mv_coo_from_fmt

    !
    ! Convert internal indices
    !
    procedure, pass(a) :: cp_to_lcoo     => psb_z_cp_coo_to_lcoo   
    procedure, pass(a) :: cp_from_lcoo   => psb_z_cp_coo_from_lcoo 
    
    procedure, pass(a) :: csput_a      => psb_z_coo_csput_a
    procedure, pass(a) :: get_diag     => psb_z_coo_get_diag
    procedure, pass(a) :: csgetrow     => psb_z_coo_csgetrow
    procedure, pass(a) :: csgetptn     => psb_z_coo_csgetptn
    procedure, pass(a) :: reinit       => psb_z_coo_reinit
    procedure, pass(a) :: get_nz_row   => psb_z_coo_get_nz_row
    procedure, pass(a) :: fix          => psb_z_fix_coo
    procedure, pass(a) :: trim         => psb_z_coo_trim
    procedure, pass(a) :: clean_zeros  => psb_z_coo_clean_zeros
    procedure, pass(a) :: print        => psb_z_coo_print
    procedure, pass(a) :: free         => z_coo_free
    procedure, pass(a) :: mold         => psb_z_coo_mold
    procedure, pass(a) :: is_sorted    => z_coo_is_sorted
    procedure, pass(a) :: is_by_rows   => z_coo_is_by_rows
    procedure, pass(a) :: is_by_cols   => z_coo_is_by_cols
    procedure, pass(a) :: set_by_rows  => z_coo_set_by_rows
    procedure, pass(a) :: set_by_cols  => z_coo_set_by_cols
    procedure, pass(a) :: set_sort_status => z_coo_set_sort_status
    procedure, pass(a) :: get_sort_status => z_coo_get_sort_status

    !
    ! This is COO specific
    !
    procedure, pass(a) :: set_nzeros   => z_coo_set_nzeros
    
    !
    ! Transpose methods. These are the base of all
    ! indirection in transpose, together with conversions
    ! they are sufficient for all cases. 
    !
    procedure, pass(a) :: transp_1mat => z_coo_transp_1mat
    procedure, pass(a) :: transc_1mat => z_coo_transc_1mat

    !
    ! Computational methods. 
    !    
    procedure, pass(a) :: csmm       => psb_z_coo_csmm
    procedure, pass(a) :: csmv       => psb_z_coo_csmv
    procedure, pass(a) :: inner_cssm => psb_z_coo_cssm
    procedure, pass(a) :: inner_cssv => psb_z_coo_cssv
    procedure, pass(a) :: scals      => psb_z_coo_scals
    procedure, pass(a) :: scalv      => psb_z_coo_scal
    procedure, pass(a) :: maxval     => psb_z_coo_maxval
    procedure, pass(a) :: spnmi      => psb_z_coo_csnmi
    procedure, pass(a) :: spnm1      => psb_z_coo_csnm1
    procedure, pass(a) :: rowsum     => psb_z_coo_rowsum
    procedure, pass(a) :: arwsum     => psb_z_coo_arwsum
    procedure, pass(a) :: colsum     => psb_z_coo_colsum
    procedure, pass(a) :: aclsum     => psb_z_coo_aclsum
    
  end type psb_z_coo_sparse_mat
  
  private :: z_coo_get_nzeros, z_coo_set_nzeros, &
       & z_coo_get_fmt,  z_coo_free, z_coo_sizeof, &
       & z_coo_transp_1mat, z_coo_transc_1mat
  
  
  !> \namespace  psb_base_mod  \class  psb_lz_base_sparse_mat
  !! \extends psb_lbase_mat_mod::psb_lbase_sparse_mat
  !! The psb_lz_base_sparse_mat type, extending psb_base_sparse_mat,
  !! defines a middle level  complex(psb_dpk_) sparse matrix object.
  !! This class object itself does not have any additional members
  !! with respect to those of the base class. Most methods cannot be fully
  !! implemented at this level, but we can define the interface for the
  !! computational methods requiring the knowledge of the underlying
  !! field, such as the matrix-vector product; this interface is defined,
  !! but is supposed to be overridden at the leaf level.
  !!
  !! About the method MOLD: this has been defined for those compilers
  !! not yet supporting ALLOCATE( ...,MOLD=...); it's otherwise silly to
  !! duplicate "by hand" what is specified in the language (in this case F2008)
  !!
  type, extends(psb_lbase_sparse_mat) :: psb_lz_base_sparse_mat
  contains
    !
    ! Data management methods: defined here, but (mostly) not implemented.
    !    
    procedure, pass(a) :: csput_a       => psb_lz_base_csput_a
    procedure, pass(a) :: csput_v       => psb_lz_base_csput_v 
    generic, public    :: csput         => csput_a,  csput_v
    procedure, pass(a) :: csgetrow      => psb_lz_base_csgetrow
    procedure, pass(a) :: csgetblk      => psb_lz_base_csgetblk
    procedure, pass(a) :: get_diag      => psb_lz_base_get_diag
    generic, public    :: csget         => csgetrow, csgetblk 
    procedure, pass(a) :: tril          => psb_lz_base_tril
    procedure, pass(a) :: triu          => psb_lz_base_triu
    procedure, pass(a) :: csclip        => psb_lz_base_csclip 
    procedure, pass(a) :: cp_to_coo     => psb_lz_base_cp_to_coo   
    procedure, pass(a) :: cp_from_coo   => psb_lz_base_cp_from_coo 
    procedure, pass(a) :: cp_to_fmt     => psb_lz_base_cp_to_fmt   
    procedure, pass(a) :: cp_from_fmt   => psb_lz_base_cp_from_fmt 
    procedure, pass(a) :: mv_to_coo     => psb_lz_base_mv_to_coo   
    procedure, pass(a) :: mv_from_coo   => psb_lz_base_mv_from_coo 
    procedure, pass(a) :: mv_to_fmt     => psb_lz_base_mv_to_fmt   
    procedure, pass(a) :: mv_from_fmt   => psb_lz_base_mv_from_fmt 
    procedure, pass(a) :: mold          => psb_lz_base_mold 
    procedure, pass(a) :: clone         => psb_lz_base_clone
    procedure, pass(a) :: make_nonunit  => psb_lz_base_make_nonunit
    procedure, pass(a) :: clean_zeros   => psb_lz_base_clean_zeros

    !
    ! Computational methods: defined here but not implemented. 
    !    
    procedure, pass(a) :: scals       => psb_lz_base_scals
    procedure, pass(a) :: scalv       => psb_lz_base_scal
    generic, public    :: scal        => scals, scalv
    procedure, pass(a) :: maxval      => psb_lz_base_maxval
    procedure, pass(a) :: spnmi       => psb_lz_base_csnmi
    procedure, pass(a) :: spnm1       => psb_lz_base_csnm1
    procedure, pass(a) :: rowsum      => psb_lz_base_rowsum
    procedure, pass(a) :: arwsum      => psb_lz_base_arwsum
    procedure, pass(a) :: colsum      => psb_lz_base_colsum
    procedure, pass(a) :: aclsum      => psb_lz_base_aclsum
    !
    ! Convert internal indices
    !
    procedure, pass(a) :: cp_to_icoo     => psb_lz_base_cp_to_icoo   
    procedure, pass(a) :: cp_from_icoo   => psb_lz_base_cp_from_icoo 
    procedure, pass(a) :: cp_to_ifmt     => psb_lz_base_cp_to_ifmt   
    procedure, pass(a) :: cp_from_ifmt   => psb_lz_base_cp_from_ifmt 
    procedure, pass(a) :: mv_to_icoo     => psb_lz_base_mv_to_icoo   
    procedure, pass(a) :: mv_from_icoo   => psb_lz_base_mv_from_icoo 
    procedure, pass(a) :: mv_to_ifmt     => psb_lz_base_mv_to_ifmt   
    procedure, pass(a) :: mv_from_ifmt   => psb_lz_base_mv_from_ifmt 
    
    !
    ! Transpose methods: defined here but not implemented. 
    !    
    procedure, pass(a) :: transp_1mat => psb_lz_base_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_lz_base_transp_2mat
    procedure, pass(a) :: transc_1mat => psb_lz_base_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_lz_base_transc_2mat
    
  end type psb_lz_base_sparse_mat
  
  private :: lz_base_mat_sync, lz_base_mat_is_host, lz_base_mat_is_dev, &
       & lz_base_mat_is_sync, lz_base_mat_set_host, lz_base_mat_set_dev,&
       & lz_base_mat_set_sync
  
  !> \namespace  psb_base_mod  \class  psb_lz_coo_sparse_mat
  !! \extends psb_lz_base_mat_mod::psb_lz_base_sparse_mat
  !! 
  !! psb_lz_coo_sparse_mat type and the related methods. This is the
  !! reference type for all the format transitions, copies and mv unless
  !! methods are implemented that allow the direct transition from one
  !! format to another. It is defined here since all other classes must
  !! refer to it per the MEDIATOR design pattern.
  !!
  type, extends(psb_lz_base_sparse_mat) :: psb_lz_coo_sparse_mat
    !> Number of nonzeros.
    integer(psb_lpk_) :: nnz
    !> Row indices.
    integer(psb_lpk_), allocatable :: ia(:)
    !> Column indices.
    integer(psb_lpk_), allocatable :: ja(:)
    !> Coefficient values. 
    complex(psb_dpk_), allocatable :: val(:)

    integer, private   :: sort_status=psb_unsorted_
    
  contains
    !
    ! Data management methods. 
    !    
    procedure, pass(a) :: get_size     => lz_coo_get_size
    procedure, pass(a) :: get_nzeros   => lz_coo_get_nzeros
    procedure, nopass  :: get_fmt      => lz_coo_get_fmt
    procedure, pass(a) :: sizeof       => lz_coo_sizeof
    procedure, pass(a) :: reallocate_nz => psb_lz_coo_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_lz_coo_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_lz_cp_coo_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_lz_cp_coo_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_lz_cp_coo_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_lz_cp_coo_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_lz_mv_coo_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_lz_mv_coo_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_lz_mv_coo_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_lz_mv_coo_from_fmt
    procedure, pass(a) :: cp_to_icoo   => psb_lz_cp_coo_to_icoo
    procedure, pass(a) :: cp_from_icoo => psb_lz_cp_coo_from_icoo
    
    procedure, pass(a) :: csput_a      => psb_lz_coo_csput_a
    procedure, pass(a) :: get_diag     => psb_lz_coo_get_diag
    procedure, pass(a) :: csgetrow     => psb_lz_coo_csgetrow
    procedure, pass(a) :: csgetptn     => psb_lz_coo_csgetptn
    procedure, pass(a) :: reinit       => psb_lz_coo_reinit
    procedure, pass(a) :: get_nz_row   => psb_lz_coo_get_nz_row
    procedure, pass(a) :: fix          => psb_lz_fix_coo
    procedure, pass(a) :: trim         => psb_lz_coo_trim
    procedure, pass(a) :: clean_zeros  => psb_lz_coo_clean_zeros
    procedure, pass(a) :: print        => psb_lz_coo_print
    procedure, pass(a) :: free         => lz_coo_free
    procedure, pass(a) :: mold         => psb_lz_coo_mold
    procedure, pass(a) :: is_sorted    => lz_coo_is_sorted
    procedure, pass(a) :: is_by_rows   => lz_coo_is_by_rows
    procedure, pass(a) :: is_by_cols   => lz_coo_is_by_cols
    procedure, pass(a) :: set_by_rows  => lz_coo_set_by_rows
    procedure, pass(a) :: set_by_cols  => lz_coo_set_by_cols
    procedure, pass(a) :: set_sort_status => lz_coo_set_sort_status
    procedure, pass(a) :: get_sort_status => lz_coo_get_sort_status

    
    ! Computational methods: defined here but not implemented. 
    !    
    procedure, pass(a) :: scals      => psb_lz_coo_scals
    procedure, pass(a) :: scalv      => psb_lz_coo_scal
    procedure, pass(a) :: maxval     => psb_lz_coo_maxval
    procedure, pass(a) :: spnmi      => psb_lz_coo_csnmi
    procedure, pass(a) :: spnm1      => psb_lz_coo_csnm1
    procedure, pass(a) :: rowsum     => psb_lz_coo_rowsum
    procedure, pass(a) :: arwsum     => psb_lz_coo_arwsum
    procedure, pass(a) :: colsum     => psb_lz_coo_colsum
    procedure, pass(a) :: aclsum     => psb_lz_coo_aclsum

    !
    ! This is COO specific
    !
    procedure, pass(a) :: set_nzeros   => lz_coo_set_nzeros
    
    !
    ! Transpose methods. These are the base of all
    ! indirection in transpose, together with conversions
    ! they are sufficient for all cases. 
    !
    procedure, pass(a) :: transp_1mat => lz_coo_transp_1mat
    procedure, pass(a) :: transc_1mat => lz_coo_transc_1mat

    
  end type psb_lz_coo_sparse_mat
  
  private :: lz_coo_get_nzeros, lz_coo_set_nzeros, &
       & lz_coo_get_fmt,  lz_coo_free, lz_coo_sizeof, &
       & lz_coo_transp_1mat, lz_coo_transc_1mat
  
  
  ! == =================
  !
  ! BASE interfaces
  !
  ! == =================

  !> Function  csput:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Insert coefficients. 
  !!
  !!
  !!         Given  a list of NZ triples
  !!           (IA(i),JA(i),VAL(i))
  !!         record a new coefficient in A such that
  !!            A(IA(1:nz),JA(1:nz)) = VAL(1:NZ).
  !!            
  !!         The internal components IA,JA,VAL are reallocated as necessary.
  !!         Constraints:
  !!         - If the matrix A is in the BUILD state, then the method will
  !!           only work for COO matrices, all other format will throw an error.
  !!           In this case coefficients are queued inside A for further processing.
  !!         - If the matrix A is in the UPDATE state, then it can be in any format;
  !!           the update operation will perform either
  !!               A(IA(1:nz),JA(1:nz)) = VAL(1:NZ)
  !!           or
  !!               A(IA(1:nz),JA(1:nz)) =  A(IA(1:nz),JA(1:nz))+VAL(1:NZ)
  !!           according to the value of DUPLICATE.
  !!         - Coefficients with (IA(I),JA(I)) outside the ranges specified by
  !!           IMIN:IMAX,JMIN:JMAX will be ignored. 
  !!           
  !!  \param nz    number of triples in input
  !!  \param ia(:)  the input row indices
  !!  \param ja(:)  the input col indices
  !!  \param val(:)  the input coefficients
  !!  \param imin  minimum row index 
  !!  \param imax  maximum row index 
  !!  \param jmin  minimum col index 
  !!  \param jmax  maximum col index 
  !!  \param info  return code
  !!  \param gtl(:) [none] an array to renumber indices   (iren(ia(:)),iren(ja(:))
  !!
  !
  interface 
    subroutine psb_z_base_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_csput_a
  end interface
  
  interface 
    subroutine psb_z_base_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_vect_type), intent(inout)  :: val
      class(psb_i_base_vect_type), intent(inout)  :: ia, ja
      integer(psb_ipk_), intent(in)             :: nz, imin, imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_csput_v
  end interface
  
  !
  !
  !> Function  csgetrow:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Get a (subset of) row(s)
  !!        
  !!        getrow is the basic method by which the other (getblk, clip) can
  !!        be implemented.
  !!        
  !!        Returns the set
  !!           NZ, IA(1:nz), JA(1:nz), VAL(1:NZ)
  !!         each identifying the position of a nonzero in A
  !!         between row indices IMIN:IMAX; 
  !!         IA,JA are reallocated as necessary.
  !!         
  !!  \param imin  the minimum row index we are interested in 
  !!  \param imax  the minimum row index we are interested in 
  !!  \param nz the number of output coefficients
  !!  \param ia(:)  the output row indices
  !!  \param ja(:)  the output col indices
  !!  \param val(:)  the output coefficients
  !!  \param info  return code
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_z_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_z_base_csgetrow
  end interface
  
  !
  !> Function  csgetblk:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Get a (subset of) row(s)
  !!        
  !!        getblk is very similar to getrow, except that the output
  !!        is packaged in a psb_z_coo_sparse_mat object
  !!         
  !!  \param imin  the minimum row index we are interested in 
  !!  \param imax  the minimum row index we are interested in 
  !!  \param b     the output (sub)matrix
  !!  \param info  return code
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_z_base_csgetblk(imin,imax,a,b,info,&
         & jmin,jmax,iren,append,rscale,cscale,chksz)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_z_base_csgetblk
  end interface
  
  !
  !
  !> Function  csclip:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Get a submatrix.
  !!        
  !!        csclip is practically identical to getblk.
  !!        One of them has to go away.....
  !!         
  !!  \param b     the output submatrix
  !!  \param info  return code
  !!  \param imin [1] the minimum row index we are interested in 
  !!  \param imax [a%get_nrows()] the minimum row index we are interested in 
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_z_base_csclip(a,b,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(out) :: b
      integer(psb_ipk_),intent(out)            :: info
      integer(psb_ipk_), intent(in), optional  :: imin,imax,jmin,jmax
      logical, intent(in), optional            :: rscale,cscale
    end subroutine psb_z_base_csclip
  end interface
  !
  !> Function  tril:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief  Copy the lower triangle, i.e. all entries
  !!         A(I,J) such that J-I <= DIAG
  !!         default value is DIAG=0, i.e. lower triangle up to
  !!         the main diagonal.
  !!         DIAG=-1 means copy the strictly lower triangle
  !!         DIAG= 1 means copy the lower triangle plus the first diagonal
  !!                 of the upper triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!         
  !!  \param l     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in 
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in 
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param u  [none]  copy of the complementary triangle
  !!           
  !
  interface 
    subroutine psb_z_base_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(out) :: l
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_ipk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_z_coo_sparse_mat), optional, intent(out) :: u
    end subroutine psb_z_base_tril
  end interface
  
  !
  !> Function  triu:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief  Copy the upper triangle, i.e. all entries
  !!         A(I,J) such that DIAG <= J-I
  !!         default value is DIAG=0, i.e. upper triangle from 
  !!         the main diagonal up.
  !!         DIAG= 1 means copy the strictly upper triangle
  !!         DIAG=-1 means copy the upper triangle plus the first diagonal
  !!                 of the lower triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!         Optionally copies the lower triangle at the same time
  !!         
  !!  \param u     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in 
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in 
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param l  [none]  copy of the complementary triangle
  !!           
  !
  interface 
    subroutine psb_z_base_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(out) :: u
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_ipk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_z_coo_sparse_mat), optional, intent(out) :: l
    end subroutine psb_z_base_triu
  end interface
  
  
  !
  !> Function  get_diag:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Extract the diagonal of A. 
  !!        
  !!   D(i) = A(i:i), i=1:min(nrows,ncols)
  !!
  !! \param d(:)  The output diagonal
  !! \param info  return code. 
  ! 
  interface 
    subroutine psb_z_base_get_diag(a,d,info) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_get_diag
  end interface
  
  !
  !> Function  mold:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Allocate a class(psb_z_base_sparse_mat) with the
  !!     same dynamic type as the input.
  !!     This is equivalent to allocate(  mold=  ) and is provided
  !!     for those compilers not yet supporting mold.
  !!   \param b The output variable
  !!   \param info return code
  ! 
  interface 
    subroutine psb_z_base_mold(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(in)                 :: a
      class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_z_base_mold
  end interface

  !
  !
  !> Function  clone:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Allocate and clone  a class(psb_z_base_sparse_mat) with the
  !!     same dynamic type as the input. 
  !!     This is equivalent to allocate( source=  ) except that
  !!     it should guarantee a deep copy wherever needed.
  !!     Should also be equivalent to calling mold and then copy,
  !!     but it can also be implemented by default using cp_to_fmt.
  !!   \param b The output variable
  !!   \param info return code
  ! 
  interface 
    subroutine psb_z_base_clone(a,b, info)
      import 
      implicit none 
      class(psb_z_base_sparse_mat), intent(inout)              :: a
      class(psb_z_base_sparse_mat), allocatable, intent(inout) :: b
      integer(psb_ipk_), intent(out)                           :: info      
    end subroutine psb_z_base_clone
  end interface


  !
  !
  !> Function  make_nonunit:
  !! \memberof  psb_z_base_make_nonunit
  !! \brief Given a matrix for which is_unit() is true, explicitly
  !!     store the unit diagonal and set is_unit() to false. 
  !!     This is needed e.g. when scaling
  ! 
  interface 
    subroutine psb_z_base_make_nonunit(a)
      import 
      implicit none 
      class(psb_z_base_sparse_mat), intent(inout) :: a
    end subroutine psb_z_base_make_nonunit
  end interface

  
  !
  !> Function  cp_to_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert to psb_z_coo_sparse_mat
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_to_coo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_to_coo
  end interface
  
  !
  !> Function  cp_from_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert from psb_z_coo_sparse_mat
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_from_coo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)     :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_from_coo
  end interface
  
  !
  !> Function  cp_to_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert to a class(psb_z_base_sparse_mat)
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%cp_to_coo(tmp) and then b%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_to_fmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_to_fmt
  end interface
  
  !
  !> Function  cp_from_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert from a class(psb_z_base_sparse_mat)
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%cp_to_coo(tmp) and then a%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_from_fmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_from_fmt
  end interface
  
  !
  !> Function  mv_to_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert to psb_z_coo_sparse_mat, freeing the source.
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_to_coo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_to_coo
  end interface
  
  !
  !> Function  mv_from_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert from psb_z_coo_sparse_mat, freeing the source.
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_from_coo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_from_coo
  end interface
  
  !
  !> Function  mv_to_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert to a class(psb_z_base_sparse_mat), freeing the source.
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%mv_to_coo(tmp) and then b%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_to_fmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_to_fmt
  end interface
  
  !
  !> Function  mv_from_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert from a class(psb_z_base_sparse_mat), freeing the source.
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%mv_to_coo(tmp) and then a%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_from_fmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_from_fmt
  end interface
  !
  !> Function  cp_to_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert to psb_z_coo_sparse_mat
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_to_lcoo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_to_lcoo
  end interface
  
  !
  !> Function  cp_from_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert from psb_z_coo_sparse_mat
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_from_lcoo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(in)     :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_from_lcoo
  end interface
  
  !
  !> Function  cp_to_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert to a class(psb_z_base_sparse_mat)
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%cp_to_coo(tmp) and then b%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_to_lfmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_to_lfmt
  end interface
  
  !
  !> Function  cp_from_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Copy and convert from a class(psb_z_base_sparse_mat)
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%cp_to_coo(tmp) and then a%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_cp_from_lfmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_cp_from_lfmt
  end interface
  
  !
  !> Function  mv_to_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert to psb_z_coo_sparse_mat, freeing the source.
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_to_lcoo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_to_lcoo
  end interface
  
  !
  !> Function  mv_from_coo:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert from psb_z_coo_sparse_mat, freeing the source.
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_from_lcoo(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_from_lcoo
  end interface
  
  !
  !> Function  mv_to_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert to a class(psb_z_base_sparse_mat), freeing the source.
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%mv_to_coo(tmp) and then b%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_to_lfmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_to_lfmt
  end interface
  
  !
  !> Function  mv_from_fmt:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Convert from a class(psb_z_base_sparse_mat), freeing the source.
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%mv_to_coo(tmp) and then a%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_z_base_mv_from_lfmt(a,b,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_mv_from_lfmt
  end interface


  !
  !> 
  !! \memberof  psb_z_base_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_clean_zeros
  !
  interface
    subroutine  psb_z_base_clean_zeros(a, info)
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_z_base_clean_zeros
  end interface
  
  !
  !> Function  transp:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        Copyout version
  !!   \param b The output variable
  !  
   interface 
    subroutine psb_z_base_transp_2mat(a,b)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_base_sparse_mat), intent(out)    :: b
    end subroutine psb_z_base_transp_2mat
  end interface
  
  !
  !> Function  transc:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Conjugate Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        Copyout version.
  !!   \param b The output variable
  !  
  interface  
    subroutine psb_z_base_transc_2mat(a,b)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      class(psb_base_sparse_mat), intent(out)    :: b
    end subroutine psb_z_base_transc_2mat
  end interface
  
  !
  !> Function  transp:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        In-place version.
  !  
  interface 
    subroutine psb_z_base_transp_1mat(a)
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
    end subroutine psb_z_base_transp_1mat
  end interface
  
  !
  !> Function  transc:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Conjugate Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        In-place version.
  !  
  interface 
    subroutine psb_z_base_transc_1mat(a)
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
    end subroutine psb_z_base_transc_1mat
  end interface
  
  !
  !> Function  csmm:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Product by a dense rank 2 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:,:) the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:,:) the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !!
  !
  interface 
    subroutine psb_z_base_csmm(alpha,a,x,beta,y,info,trans)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_base_csmm
  end interface
  
  !> Function  csmv:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Product by a dense rank 1 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:)   the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:)   the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !!
  !
  interface 
    subroutine psb_z_base_csmv(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_base_csmv
  end interface
  
  !> Function  vect_mv:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Product by an encapsulated array type(psb_z_vect_type)
  !!
  !!        Compute
  !!           Y = alpha*op(A)*X + beta*Y
  !!        Usually the unwrapping of the encapsulated vector is done
  !!        here, so that all the derived classes need only the
  !!        versions with the standard arrays.
  !!        Must be overridden explicitly in case of non standard memory
  !!        management; an example would be external memory allocation
  !!        in attached processors such as GPUs. 
  !!
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x      the input X
  !! \param beta   Scaling factor for y
  !! \param y      the input/output  Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !!
  !
  interface 
    subroutine psb_z_base_vect_mv(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)       :: alpha, beta
      class(psb_z_base_vect_type), intent(inout) :: x
      class(psb_z_base_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)  :: trans
    end subroutine psb_z_base_vect_mv
  end interface
  
  !
  !> Function  cssm:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Triangular system solve by a dense rank 2 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A^-1)*X + beta*Y
  !!
  !!        Internal workhorse called by cssm. 
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:,:) the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:,:) the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !!
  !
  interface 
    subroutine psb_z_base_inner_cssm(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_base_inner_cssm
  end interface
  
  
  !
  !> Function  cssv:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Triangular system solve by a dense rank 1 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A^-1)*X + beta*Y
  !!
  !!        Internal workhorse called by cssv. 
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:)   the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:)   the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !! \param scale  [N] Apply a scaling on Right (R) i.e. ADX
  !!               or on the Left (L)  i.e.  DAx
  !! \param D(:)   [none] Diagonal for scaling. 
  !!
  !
  interface 
    subroutine psb_z_base_inner_cssv(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_z_base_inner_cssv
  end interface
  
  !
  !> Function  inner_vect_cssv:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Triangular system solve by
  !!        an encapsulated array type(psb_z_vect_type)
  !!
  !!        Compute
  !!           Y = alpha*op(A^-1)*X + beta*Y
  !!
  !!        Internal workhorse called by vect_cssv. 
  !!        Must be overridden explicitly in case of non standard memory
  !!        management; an example would be external memory allocation
  !!        in attached processors such as GPUs. 
  !!
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x      the input dense X
  !! \param beta   Scaling factor for y
  !! \param y     the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !
  interface 
    subroutine psb_z_base_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)       :: alpha, beta
      class(psb_z_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)  :: trans
    end subroutine psb_z_base_inner_vect_sv
  end interface
  
  !
  !> Function  cssm:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Triangular system solve by a dense rank 2 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A^-1)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:,:) the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:,:) the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !! \param scale  [N] Apply a scaling on Right (R) i.e. ADX
  !!               or on the Left (L)  i.e.  DAx
  !! \param D(:)   [none] Diagonal for scaling. 
  !!
  !
  interface 
    subroutine psb_z_base_cssm(alpha,a,x,beta,y,info,trans,scale,d)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_z_base_cssm
  end interface
  
  !
  !> Function  cssv:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Triangular system solve by a dense rank 1 array.
  !!
  !!        Compute
  !!           Y = alpha*op(A^-1)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x(:)   the input dense X
  !! \param beta   Scaling factor for y
  !! \param y(:)   the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !! \param scale  [N] Apply a scaling on Right (R) i.e. ADX
  !!               or on the Left (L)  i.e.  DAx
  !! \param D(:)   [none] Diagonal for scaling. 
  !!
  !
  interface 
    subroutine psb_z_base_cssv(alpha,a,x,beta,y,info,trans,scale,d)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_z_base_cssv
  end interface
    
  !
  !> Function  vect_cssv:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Triangular system solve by
  !!        an encapsulated array type(psb_z_vect_type)
  !!
  !!        Compute
  !!           Y = alpha*op(A^-1)*X + beta*Y
  !!
  !! \param alpha  Scaling factor for Ax
  !! \param A      the input sparse matrix
  !! \param x      the input dense X
  !! \param beta   Scaling factor for y
  !! \param y     the input/output dense Y
  !! \param info   return code
  !! \param trans  [N] Whether to use A (N), its transpose (T)
  !!               or its conjugate transpose (C)
  !! \param scale  [N] Apply a scaling on Right (R) i.e. ADX
  !!               or on the Left (L)  i.e.  DAx
  !! \param D      [none] Diagonal for scaling. 
  !!
  !
  interface 
    subroutine psb_z_base_vect_cssv(alpha,a,x,beta,y,info,trans,scale,d)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)       :: alpha, beta
      class(psb_z_base_vect_type), intent(inout) :: x,y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)  :: trans, scale
      class(psb_z_base_vect_type), optional, intent(inout)   :: d
    end subroutine psb_z_base_vect_cssv
  end interface
  
  !
  !> Function  base_scals:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Scale a matrix by a single scalar value
  !!
  !! \param d      Scaling factor 
  !! \param info   return code
  !
  interface 
    subroutine psb_z_base_scals(d,a,info) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_base_scals
  end interface
  
  !
  !> Function  base_scal:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Scale a matrix by a vector
  !!
  !! \param d(:)   Scaling vector
  !! \param info   return code
  !! \param side   [L] Scale on the Left (rows) or on the Right (columns)
  !
  interface 
    subroutine psb_z_base_scal(d,a,info,side) 
      import 
      class(psb_z_base_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_z_base_scal
  end interface
  
  !
  !> Function  base_maxval:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Maximum absolute value of all coefficients;
  !! 
  !
  interface 
    function psb_z_base_maxval(a) result(res)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_base_maxval
  end interface
  
  !
  !
  !> Function  base_csnmi:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Operator infinity norm
  !! 
  !
  interface 
    function psb_z_base_csnmi(a) result(res)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_base_csnmi
  end interface

  !
  !
  !> Function  base_csnmi:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Operator 1-norm
  !! 
  !
  interface 
    function psb_z_base_csnm1(a) result(res)
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_base_csnm1
  end interface

  !
  !
  !> Function  base_rowsum:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Sum along the rows
  !! \param d(:) The output row sums
  !! 
  !
  interface 
    subroutine psb_z_base_rowsum(d,a) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_base_rowsum
  end interface

  !
  !> Function  base_arwsum:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Absolute value sum along the rows
  !! \param d(:) The output row sums
  !! 
  interface 
    subroutine psb_z_base_arwsum(d,a) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_base_arwsum
  end interface
  
  !
  !
  !> Function  base_colsum:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Sum along the columns
  !! \param d(:) The output col sums
  !! 
  !
  interface 
    subroutine psb_z_base_colsum(d,a) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_base_colsum
  end interface

  !
  !> Function  base_aclsum:
  !! \memberof  psb_z_base_sparse_mat
  !! \brief Absolute value sum along the columns
  !! \param d(:) The output col sums
  !! 
  interface 
    subroutine psb_z_base_aclsum(d,a) 
      import 
      class(psb_z_base_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_base_aclsum
  end interface

  
  ! == ===============
  !
  ! COO interfaces
  !
  ! == ===============

  !
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_reallocate_nz
  !
  interface
    subroutine  psb_z_coo_reallocate_nz(nz,a) 
      import 
      integer(psb_ipk_), intent(in) :: nz
      class(psb_z_coo_sparse_mat), intent(inout) :: a
    end subroutine psb_z_coo_reallocate_nz
  end interface
  
  !
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_reinit
  !
  interface 
    subroutine psb_z_coo_reinit(a,clear)
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_z_coo_reinit
  end interface
  !
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_trim
  !
  interface
    subroutine  psb_z_coo_trim(a)
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
    end subroutine psb_z_coo_trim
  end interface
  !
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_clean_zeros
  !
  interface
    subroutine  psb_z_coo_clean_zeros(a,info)
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_z_coo_clean_zeros
  end interface
  
  !
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_allocate_mnnz
  !
  interface
    subroutine  psb_z_coo_allocate_mnnz(m,n,a,nz) 
      import 
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_z_coo_allocate_mnnz
  end interface

  
  !> \memberof psb_z_coo_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface 
    subroutine psb_z_coo_mold(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(in)                  :: a
      class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_z_coo_mold
  end interface
  
  
  !
  !> Function print.
  !! \memberof  psb_z_coo_sparse_mat
  !! \brief Print the matrix to file in MatrixMarket format
  !!
  !! \param iout  The unit to write to
  !! \param iv    [none] Renumbering for both rows and columns
  !! \param head  [none] Descriptive header for the file
  !! \param ivr   [none] Row renumbering
  !! \param ivc   [none] Col renumbering
  !!
  !
  interface
    subroutine psb_z_coo_print(iout,a,iv,head,ivr,ivc)
      import 
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_z_coo_sparse_mat), intent(in) :: a   
      integer(psb_ipk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_z_coo_print
  end interface
  
  
  
  !
  !> Function get_nz_row.
  !! \memberof  psb_z_coo_sparse_mat
  !! \brief How many nonzeros in a row?
  !!
  !! \param idx  The row to search.
  !!
  !
  interface 
    function  psb_z_coo_get_nz_row(idx,a) result(res)
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: idx
      integer(psb_ipk_) :: res
    end function psb_z_coo_get_nz_row
  end interface
  
  
  !
  !> Funtion: fix_coo_inner
  !! \brief Make sure the entries are sorted and duplicates are handled.
  !!   Used internally by fix_coo
  !! \param nzin  Number of entries on input to be  handled
  !! \param dupl  What to do with duplicated entries.
  !! \param ia(:) Row indices
  !! \param ja(:) Col indices
  !! \param val(:) Coefficients
  !! \param nzout  Number of entries after sorting/duplicate handling
  !! \param info   return code
  !! \param idir [psb_row_major_] Sort in row major order or col major order 
  !! 
  !
  interface 
    subroutine psb_z_fix_coo_inner(nr,nc,nzin,dupl,ia,ja,val,nzout,info,idir) 
      import 
      integer(psb_ipk_), intent(in)           :: nr,nc,nzin,dupl
      integer(psb_ipk_), intent(inout)        :: ia(:), ja(:)
      complex(psb_dpk_), intent(inout) :: val(:)
      integer(psb_ipk_), intent(out)          :: nzout
      integer(psb_ipk_), intent(out)          :: info
      integer(psb_ipk_), intent(in), optional :: idir
    end subroutine psb_z_fix_coo_inner
  end interface
  
  !
  !> Function fix_coo
  !! \memberof  psb_z_coo_sparse_mat
  !! \brief Make sure the entries are sorted and duplicates are handled.
  !! \param info   return code
  !! \param idir [psb_row_major_] Sort in row major order or col major order 
  !!
  !
  interface 
    subroutine psb_z_fix_coo(a,info,idir) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), intent(in), optional :: idir
    end subroutine psb_z_fix_coo
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_to_coo
  interface 
    subroutine psb_z_cp_coo_to_coo(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_cp_coo_to_coo
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_from_coo
  interface 
    subroutine psb_z_cp_coo_from_coo(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_z_cp_coo_from_coo
  end interface
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_to_coo
  interface 
    subroutine psb_z_cp_coo_to_lcoo(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_cp_coo_to_lcoo
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_from_coo
  interface 
    subroutine psb_z_cp_coo_from_lcoo(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_z_cp_coo_from_lcoo
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_from_coo
  !! 
  interface 
    subroutine psb_z_cp_coo_to_fmt(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(in)   :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_z_cp_coo_to_fmt
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cp_from_fmt
  !! 
   interface 
    subroutine psb_z_cp_coo_from_fmt(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_z_cp_coo_from_fmt
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_to_coo
  interface 
    subroutine psb_z_mv_coo_to_coo(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_mv_coo_to_coo
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_from_coo
  interface 
    subroutine psb_z_mv_coo_from_coo(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_z_mv_coo_from_coo
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_to_fmt
  interface 
    subroutine psb_z_mv_coo_to_fmt(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_z_mv_coo_to_fmt
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_mv_from_fmt
  interface 
    subroutine psb_z_mv_coo_from_fmt(a,b,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout)  :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_z_mv_coo_from_fmt
  end interface
  
  interface 
    subroutine psb_z_coo_cp_from(a,b)
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      type(psb_z_coo_sparse_mat), intent(in)   :: b
    end subroutine psb_z_coo_cp_from
  end interface
  
  interface 
    subroutine psb_z_coo_mv_from(a,b)
      import 
      class(psb_z_coo_sparse_mat), intent(inout)  :: a
      type(psb_z_coo_sparse_mat), intent(inout) :: b
    end subroutine psb_z_coo_mv_from
  end interface
  
  
  !> Function csput
  !! \memberof  psb_z_coo_sparse_mat
  !! \brief  Add coefficients into the matrix.
  !!
  !! \param nz  Number of entries to be added
  !! \param ia(:) Row indices
  !! \param ja(:) Col indices
  !! \param val(:) Values
  !! \param imin  Minimum row index to accept
  !! \param imax  Maximum row index to accept
  !! \param jmin  Minimum col index to accept
  !! \param jmax  Maximum col index to accept
  !! \param info return code
  !! \param gtl [none] Renumbering for rows/columns
  !!
  !
  interface 
    subroutine psb_z_coo_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_ipk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_coo_csput_a
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_csgetptn
  interface 
    subroutine psb_z_coo_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_z_coo_csgetptn
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csgetrow
  interface 
    subroutine psb_z_coo_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      integer(psb_ipk_), intent(in)                  :: imin,imax
      integer(psb_ipk_), intent(out)                 :: nz
      integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_ipk_), intent(in), optional        :: iren(:)
      integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale,chksz
    end subroutine psb_z_coo_csgetrow
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cssv
  interface 
    subroutine psb_z_coo_cssv(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_coo_cssv
  end interface
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_cssm
  interface 
    subroutine psb_z_coo_cssm(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_coo_cssm
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csmv
  interface 
    subroutine psb_z_coo_csmv(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_coo_csmv
  end interface

  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csmm
  interface 
    subroutine psb_z_coo_csmm(alpha,a,x,beta,y,info,trans) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      integer(psb_ipk_), intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_z_coo_csmm
  end interface
  
    
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_maxval
  interface 
    function psb_z_coo_maxval(a) result(res)
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_coo_maxval
  end interface

  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csnmi
  interface 
    function psb_z_coo_csnmi(a) result(res)
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_coo_csnmi
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_csnm1
  interface 
    function psb_z_coo_csnm1(a) result(res)
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_z_coo_csnm1
  end interface

  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_rowsum
  interface 
    subroutine psb_z_coo_rowsum(d,a) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_coo_rowsum
  end interface
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_arwsum
  interface 
    subroutine psb_z_coo_arwsum(d,a) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_coo_arwsum
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_colsum
  interface 
    subroutine psb_z_coo_colsum(d,a) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_coo_colsum
  end interface

  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_aclsum
  interface 
    subroutine psb_z_coo_aclsum(d,a) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_z_coo_aclsum
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_get_diag
  interface 
    subroutine psb_z_coo_get_diag(a,d,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_coo_get_diag
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_scal
  interface 
    subroutine psb_z_coo_scal(d,a,info,side) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_z_coo_scal
  end interface
  
  !> 
  !! \memberof  psb_z_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_z_base_scals
  interface
    subroutine psb_z_coo_scals(d,a,info) 
      import 
      class(psb_z_coo_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_z_coo_scals
  end interface
  
  ! == =================
  !
  ! BASE interfaces
  !
  ! == =================

  !> Function  csput:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Insert coefficients. 
  !!
  !!
  !!         Given  a list of NZ triples
  !!           (IA(i),JA(i),VAL(i))
  !!         record a new coefficient in A such that
  !!            A(IA(1:nz),JA(1:nz)) = VAL(1:NZ).
  !!            
  !!         The internal components IA,JA,VAL are reallocated as necessary.
  !!         Constraints:
  !!         - If the matrix A is in the BUILD state, then the method will
  !!           only work for COO matrices, all other format will throw an error.
  !!           In this case coefficients are queued inside A for further processing.
  !!         - If the matrix A is in the UPDATE state, then it can be in any format;
  !!           the update operation will perform either
  !!               A(IA(1:nz),JA(1:nz)) = VAL(1:NZ)
  !!           or
  !!               A(IA(1:nz),JA(1:nz)) =  A(IA(1:nz),JA(1:nz))+VAL(1:NZ)
  !!           according to the value of DUPLICATE.
  !!         - Coefficients with (IA(I),JA(I)) outside the ranges specified by
  !!           IMIN:IMAX,JMIN:JMAX will be ignored. 
  !!           
  !!  \param nz    number of triples in input
  !!  \param ia(:)  the input row indices
  !!  \param ja(:)  the input col indices
  !!  \param val(:)  the input coefficients
  !!  \param imin  minimum row index 
  !!  \param imax  maximum row index 
  !!  \param jmin  minimum col index 
  !!  \param jmax  maximum col index 
  !!  \param info  return code
  !!  \param gtl(:) [none] an array to renumber indices   (iren(ia(:)),iren(ja(:))
  !!
  !
  interface 
    subroutine psb_lz_base_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_lpk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_csput_a
  end interface
  
  interface 
    subroutine psb_lz_base_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_vect_type), intent(inout)  :: val
      class(psb_l_base_vect_type), intent(inout)  :: ia, ja
      integer(psb_lpk_), intent(in)             :: nz, imin, imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_csput_v
  end interface
  
  !
  !
  !> Function  csgetrow:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Get a (subset of) row(s)
  !!        
  !!        getrow is the basic method by which the other (getblk, clip) can
  !!        be implemented.
  !!        
  !!        Returns the set
  !!           NZ, IA(1:nz), JA(1:nz), VAL(1:NZ)
  !!         each identifying the position of a nonzero in A
  !!         between row indices IMIN:IMAX; 
  !!         IA,JA are reallocated as necessary.
  !!         
  !!  \param imin  the minimum row index we are interested in 
  !!  \param imax  the minimum row index we are interested in 
  !!  \param nz the number of output coefficients
  !!  \param ia(:)  the output row indices
  !!  \param ja(:)  the output col indices
  !!  \param val(:)  the output coefficients
  !!  \param info  return code
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_lz_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_base_csgetrow
  end interface
  
  !
  !> Function  csgetblk:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Get a (subset of) row(s)
  !!        
  !!        getblk is very similar to getrow, except that the output
  !!        is packaged in a psb_lz_coo_sparse_mat object
  !!         
  !!  \param imin  the minimum row index we are interested in 
  !!  \param imax  the minimum row index we are interested in 
  !!  \param b     the output (sub)matrix
  !!  \param info  return code
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_lz_base_csgetblk(imin,imax,a,b,info,&
         & jmin,jmax,iren,append,rscale,cscale)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_base_csgetblk
  end interface
  
  !
  !
  !> Function  csclip:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Get a submatrix.
  !!        
  !!        csclip is practically identical to getblk.
  !!        One of them has to go away.....
  !!         
  !!  \param b     the output submatrix
  !!  \param info  return code
  !!  \param imin [1] the minimum row index we are interested in 
  !!  \param imax [a%get_nrows()] the minimum row index we are interested in 
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!           
  !
  interface 
    subroutine psb_lz_base_csclip(a,b,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(out) :: b
      integer(psb_ipk_),intent(out)            :: info
      integer(psb_lpk_), intent(in), optional  :: imin,imax,jmin,jmax
      logical, intent(in), optional            :: rscale,cscale
    end subroutine psb_lz_base_csclip
  end interface
  !
  !> Function  tril:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief  Copy the lower triangle, i.e. all entries
  !!         A(I,J) such that J-I <= DIAG
  !!         default value is DIAG=0, i.e. lower triangle up to
  !!         the main diagonal.
  !!         DIAG=-1 means copy the strictly lower triangle
  !!         DIAG= 1 means copy the lower triangle plus the first diagonal
  !!                 of the upper triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!         
  !!  \param l     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in 
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in 
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param u  [none]  copy of the complementary triangle
  !!           
  !
  interface 
    subroutine psb_lz_base_tril(a,l,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,u)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(out) :: l
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_lpk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_lz_coo_sparse_mat), optional, intent(out) :: u
    end subroutine psb_lz_base_tril
  end interface
  
  !
  !> Function  triu:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief  Copy the upper triangle, i.e. all entries
  !!         A(I,J) such that DIAG <= J-I
  !!         default value is DIAG=0, i.e. upper triangle from 
  !!         the main diagonal up.
  !!         DIAG= 1 means copy the strictly upper triangle
  !!         DIAG=-1 means copy the upper triangle plus the first diagonal
  !!                 of the lower triangle.
  !!         Moreover, apply a clipping by copying entries A(I,J) only if
  !!         IMIN<=I<=IMAX
  !!         JMIN<=J<=JMAX
  !!         Optionally copies the lower triangle at the same time
  !!         
  !!  \param u     the output (sub)matrix
  !!  \param info  return code
  !!  \param diag [0] the last diagonal (J-I) to be considered.
  !!  \param imin [1] the minimum row index we are interested in 
  !!  \param imax [a\%get_nrows()] the minimum row index we are interested in 
  !!  \param jmin [1] minimum col index 
  !!  \param jmax [a\%get_ncols()] maximum col index 
  !!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
  !!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
  !!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
  !!          ( iren cannot be specified with rscale/cscale)
  !!  \param append [false] append to ia,ja 
  !!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
  !!  \param l  [none]  copy of the complementary triangle
  !!           
  !
  interface 
    subroutine psb_lz_base_triu(a,u,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,l)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(out) :: u
      integer(psb_ipk_),intent(out)              :: info
      integer(psb_lpk_), intent(in), optional    :: diag,imin,imax,jmin,jmax
      logical, intent(in), optional              :: rscale,cscale
      class(psb_lz_coo_sparse_mat), optional, intent(out) :: l
    end subroutine psb_lz_base_triu
  end interface
  
  
  !
  !> Function  get_diag:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Extract the diagonal of A. 
  !!        
  !!   D(i) = A(i:i), i=1:min(nrows,ncols)
  !!
  !! \param d(:)  The output diagonal
  !! \param info  return code. 
  ! 
  interface 
    subroutine psb_lz_base_get_diag(a,d,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_get_diag
  end interface
  
  !
  !> Function  mold:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Allocate a class(psb_lz_base_sparse_mat) with the
  !!     same dynamic type as the input.
  !!     This is equivalent to allocate(  mold=  ) and is provided
  !!     for those compilers not yet supporting mold.
  !!   \param b The output variable
  !!   \param info return code
  ! 
  interface 
    subroutine psb_lz_base_mold(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(in)                 :: a
      class(psb_lz_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_lz_base_mold
  end interface

  !
  !
  !> Function  clone:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Allocate and clone  a class(psb_lz_base_sparse_mat) with the
  !!     same dynamic type as the input. 
  !!     This is equivalent to allocate( source=  ) except that
  !!     it should guarantee a deep copy wherever needed.
  !!     Should also be equivalent to calling mold and then copy,
  !!     but it can also be implemented by default using cp_to_fmt.
  !!   \param b The output variable
  !!   \param info return code
  ! 
  interface 
    subroutine psb_lz_base_clone(a,b, info)
      import 
      implicit none 
      class(psb_lz_base_sparse_mat), intent(inout)              :: a
      class(psb_lz_base_sparse_mat), allocatable, intent(inout) :: b
      integer(psb_ipk_), intent(out)                           :: info      
    end subroutine psb_lz_base_clone
  end interface


  !
  !
  !> Function  make_nonunit:
  !! \memberof  psb_lz_base_make_nonunit
  !! \brief Given a matrix for which is_unit() is true, explicitly
  !!     store the unit diagonal and set is_unit() to false. 
  !!     This is needed e.g. when scaling
  ! 
  interface 
    subroutine psb_lz_base_make_nonunit(a)
      import 
      implicit none 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
    end subroutine psb_lz_base_make_nonunit
  end interface

  
  !
  !> Function  cp_to_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert to psb_lz_coo_sparse_mat
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_to_coo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_to_coo
  end interface
  
  !
  !> Function  cp_from_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert from psb_lz_coo_sparse_mat
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_from_coo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(in)     :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_from_coo
  end interface
  
  !
  !> Function  cp_to_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert to a class(psb_lz_base_sparse_mat)
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%cp_to_coo(tmp) and then b%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_to_fmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_to_fmt
  end interface
  
  !
  !> Function  cp_from_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert from a class(psb_lz_base_sparse_mat)
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%cp_to_coo(tmp) and then a%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_from_fmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_from_fmt
  end interface
  
  !
  !> Function  mv_to_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert to psb_lz_coo_sparse_mat, freeing the source.
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_to_coo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_to_coo
  end interface
  
  !
  !> Function  mv_from_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert from psb_lz_coo_sparse_mat, freeing the source.
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_from_coo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_from_coo
  end interface
  
  !
  !> Function  mv_to_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert to a class(psb_lz_base_sparse_mat), freeing the source.
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%mv_to_coo(tmp) and then b%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_to_fmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_to_fmt
  end interface
  
  !
  !> Function  mv_from_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert from a class(psb_lz_base_sparse_mat), freeing the source.
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%mv_to_coo(tmp) and then a%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_from_fmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_from_fmt
  end interface

  
  !
  !> Function  cp_to_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert to psb_lz_coo_sparse_mat
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_to_icoo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_to_icoo
  end interface
  
  !
  !> Function  cp_from_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert from psb_lz_coo_sparse_mat
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_from_icoo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)     :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_from_icoo
  end interface
  
  !
  !> Function  cp_to_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert to a class(psb_lz_base_sparse_mat)
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%cp_to_coo(tmp) and then b%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_to_ifmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_to_ifmt
  end interface
  
  !
  !> Function  cp_from_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Copy and convert from a class(psb_lz_base_sparse_mat)
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%cp_to_coo(tmp) and then a%cp_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_cp_from_ifmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_cp_from_ifmt
  end interface
  
  !
  !> Function  mv_to_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert to psb_lz_coo_sparse_mat, freeing the source.
  !!        Invoked from the source object.
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_to_icoo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_to_icoo
  end interface
  
  !
  !> Function  mv_from_coo:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert from psb_lz_coo_sparse_mat, freeing the source.
  !!        Invoked from the target object.
  !!   \param b The input variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_from_icoo(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_from_icoo
  end interface
  
  !
  !> Function  mv_to_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert to a class(psb_lz_base_sparse_mat), freeing the source.
  !!        Invoked from the source object. Can be implemented by
  !!        simply invoking a%mv_to_coo(tmp) and then b%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_to_ifmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_to_ifmt
  end interface
  
  !
  !> Function  mv_from_fmt:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Convert from a class(psb_lz_base_sparse_mat), freeing the source.
  !!        Invoked from the target object. Can be implemented by
  !!        simply invoking b%mv_to_coo(tmp) and then a%mv_from_coo(tmp).
  !!   \param b The output variable
  !!   \param info return code
  !  
  interface 
    subroutine psb_lz_base_mv_from_ifmt(a,b,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      class(psb_z_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_mv_from_ifmt
  end interface



  !
  !> 
  !! \memberof  psb_lz_base_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_clean_zeros
  !
  interface
    subroutine  psb_lz_base_clean_zeros(a, info)
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_lz_base_clean_zeros
  end interface

  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_maxval
  interface 
    function psb_lz_coo_maxval(a) result(res)
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_lz_coo_maxval
  end interface

  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_csnmi
  interface 
    function psb_lz_coo_csnmi(a) result(res)
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_lz_coo_csnmi
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_csnm1
  interface 
    function psb_lz_coo_csnm1(a) result(res)
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_lz_coo_csnm1
  end interface

  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_rowsum
  interface 
    subroutine psb_lz_coo_rowsum(d,a) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_coo_rowsum
  end interface
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_arwsum
  interface 
    subroutine psb_lz_coo_arwsum(d,a) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_coo_arwsum
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_colsum
  interface 
    subroutine psb_lz_coo_colsum(d,a) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_coo_colsum
  end interface

  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_z_base_mat_mod::psb_lz_base_aclsum
  interface 
    subroutine psb_lz_coo_aclsum(d,a) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_coo_aclsum
  end interface
    
  !
  !> Function  base_scals:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Scale a matrix by a single scalar value
  !!
  !! \param d      Scaling factor 
  !! \param info   return code
  !
  interface 
    subroutine psb_lz_base_scals(d,a,info) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_base_scals
  end interface
  
  !
  !> Function  base_scal:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Scale a matrix by a vector
  !!
  !! \param d(:)   Scaling vector
  !! \param info   return code
  !! \param side   [L] Scale on the Left (rows) or on the Right (columns)
  !
  interface 
    subroutine psb_lz_base_scal(d,a,info,side) 
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_lz_base_scal
  end interface
  
  !
  !> Function  base_maxval:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Maximum absolute value of all coefficients;
  !! 
  !
  interface 
    function psb_lz_base_maxval(a) result(res)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_lz_base_maxval
  end interface
  
  !
  !
  !> Function  base_csnmi:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Operator infinity norm
  !! 
  !
  interface 
    function psb_lz_base_csnmi(a) result(res)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_lz_base_csnmi
  end interface

  !
  !
  !> Function  base_csnmi:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Operator 1-norm
  !! 
  !
  interface 
    function psb_lz_base_csnm1(a) result(res)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_lz_base_csnm1
  end interface

  !
  !
  !> Function  base_rowsum:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Sum along the rows
  !! \param d(:) The output row sums
  !! 
  !
  interface 
    subroutine psb_lz_base_rowsum(d,a) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_base_rowsum
  end interface

  !
  !> Function  base_arwsum:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Absolute value sum along the rows
  !! \param d(:) The output row sums
  !! 
  interface 
    subroutine psb_lz_base_arwsum(d,a) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_base_arwsum
  end interface
  
  !
  !
  !> Function  base_colsum:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Sum along the columns
  !! \param d(:) The output col sums
  !! 
  !
  interface 
    subroutine psb_lz_base_colsum(d,a) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_base_colsum
  end interface

  !
  !> Function  base_aclsum:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Absolute value sum along the columns
  !! \param d(:) The output col sums
  !! 
  interface 
    subroutine psb_lz_base_aclsum(d,a) 
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_lz_base_aclsum
  end interface

  
  !
  !> Function  transp:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        Copyout version
  !!   \param b The output variable
  !  
   interface 
    subroutine psb_lz_base_transp_2mat(a,b)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lbase_sparse_mat), intent(out)    :: b
    end subroutine psb_lz_base_transp_2mat
  end interface
  
  !
  !> Function  transc:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Conjugate Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        Copyout version.
  !!   \param b The output variable
  !  
  interface  
    subroutine psb_lz_base_transc_2mat(a,b)
      import 
      class(psb_lz_base_sparse_mat), intent(in) :: a
      class(psb_lbase_sparse_mat), intent(out)    :: b
    end subroutine psb_lz_base_transc_2mat
  end interface
  
  !
  !> Function  transp:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        In-place version.
  !  
  interface 
    subroutine psb_lz_base_transp_1mat(a)
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
    end subroutine psb_lz_base_transp_1mat
  end interface
  
  !
  !> Function  transc:
  !! \memberof  psb_lz_base_sparse_mat
  !! \brief Conjugate Transpose. Can always be implemented by staging through a COO
  !!        temporary for which transpose is very easy. 
  !!        In-place version.
  !  
  interface 
    subroutine psb_lz_base_transc_1mat(a)
      import 
      class(psb_lz_base_sparse_mat), intent(inout) :: a
    end subroutine psb_lz_base_transc_1mat
  end interface
  
  ! == ===============
  !
  ! COO interfaces
  !
  ! == ===============

  !
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_reallocate_nz
  !
  interface
    subroutine  psb_lz_coo_reallocate_nz(nz,a) 
      import 
      integer(psb_lpk_), intent(in) :: nz
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
    end subroutine psb_lz_coo_reallocate_nz
  end interface
  
  !
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_reinit
  !
  interface 
    subroutine psb_lz_coo_reinit(a,clear)
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_lz_coo_reinit
  end interface
  !
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_trim
  !
  interface
    subroutine  psb_lz_coo_trim(a)
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
    end subroutine psb_lz_coo_trim
  end interface
  !
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_clean_zeros
  !
  interface
    subroutine  psb_lz_coo_clean_zeros(a,info)
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_lz_coo_clean_zeros
  end interface
  
  !
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_allocate_mnnz
  !
  interface
    subroutine  psb_lz_coo_allocate_mnnz(m,n,a,nz) 
      import 
      integer(psb_lpk_), intent(in) :: m,n
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      integer(psb_lpk_), intent(in), optional :: nz
    end subroutine psb_lz_coo_allocate_mnnz
  end interface

  
  !> \memberof psb_lz_coo_sparse_mat
  !| \see psb_base_mat_mod::psb_base_mold
  interface 
    subroutine psb_lz_coo_mold(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in)                  :: a
      class(psb_lz_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_lz_coo_mold
  end interface
  
  
  !
  !> Function print.
  !! \memberof  psb_lz_coo_sparse_mat
  !! \brief Print the matrix to file in MatrixMarket format
  !!
  !! \param iout  The unit to write to
  !! \param iv    [none] Renumbering for both rows and columns
  !! \param head  [none] Descriptive header for the file
  !! \param ivr   [none] Row renumbering
  !! \param ivc   [none] Col renumbering
  !!
  !
  interface
    subroutine psb_lz_coo_print(iout,a,iv,head,ivr,ivc)
      import 
      integer(psb_ipk_), intent(in)               :: iout
      class(psb_lz_coo_sparse_mat), intent(in) :: a   
      integer(psb_lpk_), intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_lz_coo_print
  end interface
  
  
  
  !
  !> Function get_nz_row.
  !! \memberof  psb_lz_coo_sparse_mat
  !! \brief How many nonzeros in a row?
  !!
  !! \param idx  The row to search.
  !!
  !
  interface 
    function  psb_lz_coo_get_nz_row(idx,a) result(res)
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: idx
      integer(psb_lpk_) :: res
    end function psb_lz_coo_get_nz_row
  end interface
  
  
  !
  !> Funtion: fix_coo_inner
  !! \brief Make sure the entries are sorted and duplicates are handled.
  !!   Used internally by fix_coo
  !! \param nzin  Number of entries on input to be  handled
  !! \param dupl  What to do with duplicated entries.
  !! \param ia(:) Row indices
  !! \param ja(:) Col indices
  !! \param val(:) Coefficients
  !! \param nzout  Number of entries after sorting/duplicate handling
  !! \param info   return code
  !! \param idir [psb_row_major_] Sort in row major order or col major order 
  !! 
  !
  interface 
    subroutine psb_lz_fix_coo_inner(nr,nc,nzin,dupl,ia,ja,val,nzout,info,idir) 
      import 
      integer(psb_lpk_), intent(in)           :: nr,nc,nzin
      integer(psb_ipk_), intent(in)           :: dupl
      integer(psb_lpk_), intent(inout)        :: ia(:), ja(:)
      complex(psb_dpk_), intent(inout) :: val(:)
      integer(psb_lpk_), intent(out)          :: nzout
      integer(psb_ipk_), intent(out)          :: info
      integer(psb_ipk_), intent(in), optional :: idir
    end subroutine psb_lz_fix_coo_inner
  end interface
  
  !
  !> Function fix_coo
  !! \memberof  psb_lz_coo_sparse_mat
  !! \brief Make sure the entries are sorted and duplicates are handled.
  !! \param info   return code
  !! \param idir [psb_row_major_] Sort in row major order or col major order 
  !!
  !
  interface 
    subroutine psb_lz_fix_coo(a,info,idir) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)                :: info
      integer(psb_ipk_), intent(in), optional :: idir
    end subroutine psb_lz_fix_coo
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_cp_to_coo
  interface 
    subroutine psb_lz_cp_coo_to_coo(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_cp_coo_to_coo
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_cp_from_coo
  interface 
    subroutine psb_lz_cp_coo_from_coo(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_lz_cp_coo_from_coo
  end interface
  
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_cp_to_coo
  interface 
    subroutine psb_lz_cp_coo_to_icoo(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      class(psb_z_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_cp_coo_to_icoo
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_cp_from_coo
  interface 
    subroutine psb_lz_cp_coo_from_icoo(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      class(psb_z_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_lz_cp_coo_from_icoo
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_cp_from_coo
  !! 
  interface 
    subroutine psb_lz_cp_coo_to_fmt(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in)   :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_lz_cp_coo_to_fmt
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_cp_from_fmt
  !! 
   interface 
    subroutine psb_lz_cp_coo_from_fmt(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_lz_cp_coo_from_fmt
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_mv_to_coo
  interface 
    subroutine psb_lz_mv_coo_to_coo(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(inout)   :: b
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_mv_coo_to_coo
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_mv_from_coo
  interface 
    subroutine psb_lz_mv_coo_from_coo(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      class(psb_lz_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_lz_mv_coo_from_coo
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_mv_to_fmt
  interface 
    subroutine psb_lz_mv_coo_to_fmt(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      class(psb_lz_base_sparse_mat), intent(inout)  :: b
      integer(psb_ipk_), intent(out)                        :: info
    end subroutine psb_lz_mv_coo_to_fmt
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_mv_from_fmt
  interface 
    subroutine psb_lz_mv_coo_from_fmt(a,b,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout)  :: a
      class(psb_lz_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)                         :: info
    end subroutine psb_lz_mv_coo_from_fmt
  end interface
  
  interface 
    subroutine psb_lz_coo_cp_from(a,b)
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      type(psb_lz_coo_sparse_mat), intent(in)   :: b
    end subroutine psb_lz_coo_cp_from
  end interface
  
  interface 
    subroutine psb_lz_coo_mv_from(a,b)
      import 
      class(psb_lz_coo_sparse_mat), intent(inout)  :: a
      type(psb_lz_coo_sparse_mat), intent(inout) :: b
    end subroutine psb_lz_coo_mv_from
  end interface
  
  
  !> Function csput
  !! \memberof  psb_lz_coo_sparse_mat
  !! \brief  Add coefficients into the matrix.
  !!
  !! \param nz  Number of entries to be added
  !! \param ia(:) Row indices
  !! \param ja(:) Col indices
  !! \param val(:) Values
  !! \param imin  Minimum row index to accept
  !! \param imax  Maximum row index to accept
  !! \param jmin  Minimum col index to accept
  !! \param jmax  Maximum col index to accept
  !! \param info return code
  !! \param gtl [none] Renumbering for rows/columns
  !!
  !
  interface 
    subroutine psb_lz_coo_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: val(:)
      integer(psb_lpk_), intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_coo_csput_a
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_base_mat_mod::psb_base_csgetptn
  interface 
    subroutine psb_lz_coo_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_coo_csgetptn
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_csgetrow
  interface 
    subroutine psb_lz_coo_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      integer(psb_lpk_), intent(in)                  :: imin,imax
      integer(psb_lpk_), intent(out)                 :: nz
      integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer(psb_ipk_),intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer(psb_lpk_), intent(in), optional        :: iren(:)
      integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_lz_coo_csgetrow
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_get_diag
  interface 
    subroutine psb_lz_coo_get_diag(a,d,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(in) :: a
      complex(psb_dpk_), intent(out)     :: d(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_coo_get_diag
  end interface
  
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_scal
  interface 
    subroutine psb_lz_coo_scal(d,a,info,side) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d(:)
      integer(psb_ipk_), intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_lz_coo_scal
  end interface
  
  !> 
  !! \memberof  psb_lz_coo_sparse_mat
  !! \see psb_lz_base_mat_mod::psb_lz_base_scals
  interface
    subroutine psb_lz_coo_scals(d,a,info) 
      import 
      class(psb_lz_coo_sparse_mat), intent(inout) :: a
      complex(psb_dpk_), intent(in)      :: d
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_lz_coo_scals
  end interface
  
contains 
 
  
  ! == ==================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ==================================
  
  
  
  function z_coo_sizeof(a) result(res)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    res = 3*psb_sizeof_ip
    res = res + (2*psb_sizeof_dp)  * psb_size(a%val)
    res = res + psb_sizeof_ip * psb_size(a%ia)
    res = res + psb_sizeof_ip * psb_size(a%ja)
    
  end function z_coo_sizeof
  
  
  function z_coo_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'COO'
  end function z_coo_get_fmt
  
  
  function z_coo_get_size(a) result(res)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res = -1
    
    if (allocated(a%ia)) res = size(a%ia)
    if (allocated(a%ja)) then 
      if (res >= 0) then 
        res = min(res,size(a%ja))
      else 
        res = size(a%ja)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if
  end function z_coo_get_size
  
  
  function z_coo_get_nzeros(a) result(res)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: res
    res  = a%nnz
  end function z_coo_get_nzeros
  
  function z_coo_is_by_rows(a) result(res)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(in) :: a
    logical :: res
    res  = (a%sort_status == psb_row_major_)
  end function z_coo_is_by_rows
  
  function z_coo_is_by_cols(a) result(res)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(in) :: a
    logical :: res
    res  = (a%sort_status == psb_col_major_)
  end function z_coo_is_by_cols
  
  function z_coo_is_sorted(a) result(res)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(in) :: a
    logical :: res
    res  = (a%sort_status == psb_row_major_) &
         & .or.(a%sort_status == psb_col_major_)
  end function z_coo_is_sorted
  
 
  
  ! == ==================================
  !
  !
  !
  ! Setters 
  !
  !
  !
  !
  !
  !
  ! == ==================================
  
  subroutine  z_coo_set_nzeros(nz,a)
    implicit none 
    integer(psb_ipk_), intent(in) :: nz
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    a%nnz = nz
    
  end subroutine z_coo_set_nzeros
  
  function z_coo_get_sort_status(a) result(res)
    implicit none 
    integer(psb_ipk_)   :: res
    class(psb_z_coo_sparse_mat), intent(in) :: a
    
    res = a%sort_status
  end function z_coo_get_sort_status
  
  subroutine  z_coo_set_sort_status(ist,a)
    implicit none 
    integer(psb_ipk_), intent(in)   :: ist
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    a%sort_status = ist
    call a%set_sorted((a%sort_status == psb_row_major_) &
         & .or.(a%sort_status == psb_col_major_))   
  end subroutine z_coo_set_sort_status
 
  
  subroutine  z_coo_set_by_rows(a)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    a%sort_status = psb_row_major_
    call a%set_sorted()
  end subroutine z_coo_set_by_rows
 
  
  subroutine  z_coo_set_by_cols(a)
    implicit none 
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    a%sort_status = psb_col_major_
    call a%set_sorted()
  end subroutine z_coo_set_by_cols
  
  ! == ==================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ==================================
  
  subroutine  z_coo_free(a) 
    implicit none 
    
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_ipk_)
    call a%set_ncols(0_psb_ipk_)
    call a%set_nzeros(0_psb_ipk_)
    call a%set_sort_status(psb_unsorted_)
    
    return
    
  end subroutine z_coo_free
  
  
  
  ! == ==================================
  !
  !
  !
  ! Computational routines
  !
  !
  !
  !
  !
  !
  ! == ==================================
  subroutine z_coo_transp_1mat(a)
    implicit none 
    
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    integer(psb_ipk_), allocatable :: itemp(:) 
    integer(psb_ipk_) :: info
    
    call a%psb_z_base_sparse_mat%psb_base_sparse_mat%transp()
    call move_alloc(a%ia,itemp)
    call move_alloc(a%ja,a%ia)
    call move_alloc(itemp,a%ja)
        
    call a%set_sorted(.false.)
    call a%set_sort_status(psb_unsorted_)
    
    return
    
  end subroutine z_coo_transp_1mat
  
  subroutine z_coo_transc_1mat(a)
    implicit none 
    
    class(psb_z_coo_sparse_mat), intent(inout) :: a
    
    call a%transp() 
    ! This will morph into conjg() for C and Z
    ! and into a no-op for S and D, so a conditional
    ! on a constant ought to take it out completely. 
    if (psb_z_is_complex_) a%val(:) = conjg(a%val(:))

  end subroutine z_coo_transc_1mat

  
  ! == ==================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ==================================
  
  
  
  function lz_coo_sizeof(a) result(res)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res
    res = 3*psb_sizeof_lp
    res = res + (2*psb_sizeof_dp)  * psb_size(a%val)
    res = res + psb_sizeof_lp * psb_size(a%ia)
    res = res + psb_sizeof_lp * psb_size(a%ja)
    
  end function lz_coo_sizeof
  
  
  function lz_coo_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'COO'
  end function lz_coo_get_fmt
  
  
  function lz_coo_get_size(a) result(res)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: res
    res = -1
    
    if (allocated(a%ia)) res = size(a%ia)
    if (allocated(a%ja)) then 
      if (res >= 0) then 
        res = min(res,size(a%ja))
      else 
        res = size(a%ja)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if
  end function lz_coo_get_size
  
  
  function lz_coo_get_nzeros(a) result(res)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: res
    res  = a%nnz
  end function lz_coo_get_nzeros
  
  function lz_coo_is_by_rows(a) result(res)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    logical :: res
    res  = (a%sort_status == psb_row_major_)
  end function lz_coo_is_by_rows
  
  function lz_coo_is_by_cols(a) result(res)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    logical :: res
    res  = (a%sort_status == psb_col_major_)
  end function lz_coo_is_by_cols
  
  function lz_coo_is_sorted(a) result(res)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    logical :: res
    res  = (a%sort_status == psb_row_major_) &
         & .or.(a%sort_status == psb_col_major_)
  end function lz_coo_is_sorted
  
 
  
  ! == ==================================
  !
  !
  !
  ! Setters 
  !
  !
  !
  !
  !
  !
  ! == ==================================
  
  subroutine  lz_coo_set_nzeros(nz,a)
    implicit none 
    integer(psb_lpk_), intent(in) :: nz
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    a%nnz = nz
    
  end subroutine lz_coo_set_nzeros
  
  function lz_coo_get_sort_status(a) result(res)
    implicit none 
    integer(psb_ipk_)   :: res
    class(psb_lz_coo_sparse_mat), intent(in) :: a
    
    res = a%sort_status
  end function lz_coo_get_sort_status
  
  subroutine  lz_coo_set_sort_status(ist,a)
    implicit none 
    integer(psb_ipk_), intent(in)   :: ist
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    a%sort_status = ist
    call a%set_sorted((a%sort_status == psb_row_major_) &
         & .or.(a%sort_status == psb_col_major_))   
  end subroutine lz_coo_set_sort_status
 
  
  subroutine  lz_coo_set_by_rows(a)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    a%sort_status = psb_row_major_
    call a%set_sorted()
  end subroutine lz_coo_set_by_rows
 
  
  subroutine  lz_coo_set_by_cols(a)
    implicit none 
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    a%sort_status = psb_col_major_
    call a%set_sorted()
  end subroutine lz_coo_set_by_cols
  
  ! == ==================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ==================================
  
  subroutine  lz_coo_free(a) 
    implicit none 
    
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0_psb_lpk_)
    call a%set_ncols(0_psb_lpk_)
    call a%set_nzeros(0_psb_lpk_)
    call a%set_sort_status(psb_unsorted_)
    
    return
    
  end subroutine lz_coo_free
  
  
  
  ! == ==================================
  !
  !
  !
  ! Computational routines
  !
  !
  !
  !
  !
  !
  ! == ==================================
  subroutine lz_coo_transp_1mat(a)
    implicit none 
    
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    integer(psb_lpk_), allocatable :: itemp(:) 
    integer(psb_ipk_) :: info
    
    call a%psb_lz_base_sparse_mat%psb_lbase_sparse_mat%transp()
    call move_alloc(a%ia,itemp)
    call move_alloc(a%ja,a%ia)
    call move_alloc(itemp,a%ja)
        
    call a%set_sorted(.false.)
    call a%set_sort_status(psb_unsorted_)
    
    return
    
  end subroutine lz_coo_transp_1mat
  
  subroutine lz_coo_transc_1mat(a)
    implicit none 
    
    class(psb_lz_coo_sparse_mat), intent(inout) :: a
    
    call a%transp() 
    ! This will morph into conjg() for C and Z
    ! and into a no-op for S and D, so a conditional
    ! on a constant ought to take it out completely. 
    if (psb_lz_is_complex_) a%val(:) = conjg(a%val(:))

  end subroutine lz_coo_transc_1mat


end module psb_z_base_mat_mod



