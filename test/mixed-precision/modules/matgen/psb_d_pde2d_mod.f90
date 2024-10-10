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
! File: psb_d_pde2d.f90
!
! Program: psb_d_pde2d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs.
!
!
! The PDE is a general second order equation in 2d
!
!   a1 dd(u)  a2 dd(u)   b1 d(u)   b2 d(u)
! -   ------ -  ------   -----  +  ------  + c u = f
!      dxdx     dydy        dx       dy
!
! with Dirichlet boundary conditions
!   u = g
!
!  on the unit square  0<=x,y<=1.
!
!
! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
!
! There are three choices available for data distribution:
! 1. A simple BLOCK distribution
! 2. A ditribution based on arbitrary assignment of indices to processes,
!    typically from a graph partitioner
! 3. A 2D distribution in which the unit square is partitioned
!    into rectangles, each one assigned to a process.
!
module psb_d_pde2d_mod

  use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_desc_type,&
       &  psb_dspmat_type, psb_d_vect_type, dzero,&
       &  psb_d_base_sparse_mat, psb_d_base_vect_type, psb_i_base_vect_type

  interface
    function d_func_2d(x,y) result(val)
      import :: psb_dpk_
      real(psb_dpk_), intent(in) :: x,y
      real(psb_dpk_) :: val
    end function d_func_2d
  end interface

contains

  function d_null_func_2d(x,y) result(val)

    real(psb_dpk_), intent(in) :: x,y
    real(psb_dpk_) :: val

    val = dzero

  end function d_null_func_2d

  !
  ! functions parametrizing the differential equation
  !

  !
  ! Note: b1 and b2 are the coefficients of the first
  ! derivative of the unknown function. The default
  ! we apply here is to have them zero, so that the resulting
  ! matrix is symmetric/hermitian and suitable for
  ! testing with CG and FCG.
  ! When testing methods for non-hermitian matrices you can
  ! change the B1/B2 functions to e.g. done/sqrt((2*done))
  !
  function b1(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y
    b1=dzero
  end function b1
  function b2(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y
    b2=dzero
  end function b2
  function c(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y
    c=0.d0
  end function c
  function a1(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none
    real(psb_dpk_) ::  a1
    real(psb_dpk_), intent(in) :: x,y
    a1=done/80
  end function a1
  function a2(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y
    a2=done/80
  end function a2
  function g(x,y)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then
      g = exp(-y**2)
    end if
  end function g


  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs.
  !
  subroutine psb_d_gen_pde2d(ctxt,idim,a,bv,xv,desc_a,afmt,info,&
       & f,amold,vmold,imold,partition,nrl,iv)
    use psb_base_mod
    use psb_util_mod
#if defined(OPENMP)
    use omp_lib
#endif
    !
    !   Discretizes the partial differential equation
    !
    !   a1 dd(u)  a2 dd(u)    b1 d(u)  b2 d(u)
    ! -   ------ -  ------ +  -----  +  ------  + c u = f
    !      dxdx     dydy         dx       dy
    !
    ! with Dirichlet boundary conditions
    !   u = g
    !
    !  on the unit square  0<=x,y<=1.
    !
    !
    ! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
    !
    implicit none
    integer(psb_ipk_)     :: idim
    type(psb_dspmat_type) :: a
    type(psb_d_vect_type) :: xv,bv
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type)   :: ctxt
    integer(psb_ipk_)     :: info
    character(len=*)      :: afmt
    procedure(d_func_2d), optional :: f
    class(psb_d_base_sparse_mat), optional :: amold
    class(psb_d_base_vect_type), optional :: vmold
    class(psb_i_base_vect_type), optional :: imold
    integer(psb_ipk_), optional :: partition, nrl,iv(:)

    ! Local variables.

    integer(psb_ipk_), parameter :: nb=20
    type(psb_d_csc_sparse_mat)  :: acsc
    type(psb_d_coo_sparse_mat)  :: acoo
    type(psb_d_csr_sparse_mat)  :: acsr
    real(psb_dpk_)           :: zt(nb),x,y,z
    integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_, mysz
    integer(psb_lpk_) :: m,n,glob_row,nt
    integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
    ! For 2D partition
    ! Note: integer control variables going directly into an MPI call
    ! must be 4 bytes, i.e. psb_mpk_
    integer(psb_mpk_) :: npdims(2), npp, minfo
    integer(psb_ipk_) :: npx,npy,iamx,iamy,mynx,myny
    integer(psb_ipk_), allocatable :: bndx(:),bndy(:)
    ! Process grid
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: icoeff
    integer(psb_lpk_), allocatable     :: myidx(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_dpk_)            :: deltah, sqdeltah, deltah2
    real(psb_dpk_), parameter :: rhs=dzero,one=done,zero=dzero
    real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    procedure(d_func_2d), pointer :: f_
    character(len=20)  :: name, ch_err,tmpfmt

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)


    if (present(f)) then
      f_ => f
    else
      f_ => d_null_func_2d
    end if

    deltah   = done/(idim+1)
    sqdeltah = deltah*deltah
    deltah2  = (2*done)* deltah

    if (present(partition)) then
      if ((1<= partition).and.(partition <= 3)) then
        partition_ = partition
      else
        write(*,*) 'Invalid partition choice ',partition,' defaulting to 3'
        partition_ = 3
      end if
    else
      partition_ = 3
    end if

    ! initialize array descriptor and sparse matrix storage. provide an
    ! estimate of the number of non zeroes

    m   = (1_psb_lpk_)*idim*idim
    n   = m
    nnz = ((n*7)/(np))
    if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n
    t0 = psb_wtime()
    select case(partition_)
    case(1)
      ! A BLOCK partition
      if (present(nrl)) then
        nr = nrl
      else
        !
        ! Using a simple BLOCK distribution.
        !
        nt = (m+np-1)/np
        nr = max(0,min(nt,m-(iam*nt)))
      end if

      nt = nr
      call psb_sum(ctxt,nt)
      if (nt /= m) then
        write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return
      end if

      !
      ! First example  of use of CDALL: specify for each process a number of
      ! contiguous rows
      !
      call psb_cdall(ctxt,desc_a,info,nl=nr)
      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

    case(2)
      ! A  partition  defined by the user through IV

      if (present(iv)) then
        if (size(iv) /= m) then
          write(psb_err_unit,*) iam, 'Initialization error: wrong IV size',size(iv),m
          info = -1
          call psb_barrier(ctxt)
          call psb_abort(ctxt)
          return
        end if
      else
        write(psb_err_unit,*) iam, 'Initialization error: IV not present'
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return
      end if

      !
      ! Second example  of use of CDALL: specify for each row the
      ! process that owns it
      !
      call psb_cdall(ctxt,desc_a,info,vg=iv)
      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

    case(3)
      ! A 2-dimensional partition

      ! A nifty MPI function will split the process list
      npdims = 0
#if defined(SERIAL_MPI)
      npdims = 1
#else
      call mpi_dims_create(np,2,npdims,info)
#endif
      npx = npdims(1)
      npy = npdims(2)

      allocate(bndx(0:npx),bndy(0:npy))
      ! We can reuse idx2ijk for process indices as well.
      call idx2ijk(iamx,iamy,iam,npx,npy,base=0)
      ! Now let's split the 2D square in rectangles
      call dist1Didx(bndx,idim,npx)
      mynx = bndx(iamx+1)-bndx(iamx)
      call dist1Didx(bndy,idim,npy)
      myny = bndy(iamy+1)-bndy(iamy)

      ! How many indices do I own?
      nlr = mynx*myny
      allocate(myidx(nlr))
      ! Now, let's generate the list of indices I own
      nr = 0
      do i=bndx(iamx),bndx(iamx+1)-1
        do j=bndy(iamy),bndy(iamy+1)-1
          nr = nr + 1
          call ijk2idx(myidx(nr),i,j,idim,idim)
        end do
      end do
      if (nr /= nlr) then
        write(psb_err_unit,*) iam,iamx,iamy, 'Initialization error: NR vs NLR ',&
             & nr,nlr,mynx,myny
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
      end if

      !
      ! Third example  of use of CDALL: specify for each process
      ! the set of global indices it owns.
      !
      call psb_cdall(ctxt,desc_a,info,vl=myidx)

      !
      ! Specify process topology
      !
      block
        !
        ! Use adjcncy methods 
        ! 
        integer(psb_mpk_), allocatable :: neighbours(:)
        integer(psb_mpk_) :: cnt
        logical, parameter :: debug_adj=.true.
        if (debug_adj.and.(np > 1)) then 
          cnt = 0
          allocate(neighbours(np))
          if (iamx < npx-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx+1,iamy,npx,npy,base=0)
          end if
          if (iamy < npy-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy+1,npx,npy,base=0)
          end if
          if (iamx >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx-1,iamy,npx,npy,base=0)
          end if
          if (iamy >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy-1,npx,npy,base=0)
          end if
          call psb_realloc(cnt, neighbours,info)
          call desc_a%set_p_adjcncy(neighbours)
          !write(0,*) iam,' Check on neighbours: ',desc_a%get_p_adjcncy()
        end if
      end block

    case default
      write(psb_err_unit,*) iam, 'Initialization error: should not get here'
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end select


    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz,&
         & bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
    ! define  rhs from boundary conditions; also build initial guess
    if (info == psb_success_) call psb_geall(xv,desc_a,info)
    if (info == psb_success_) call psb_geall(bv,desc_a,info,&
         & bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)

    call psb_barrier(ctxt)
    talc = psb_wtime()-t0

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    call psb_barrier(ctxt)
    t1 = psb_wtime()
    !$omp parallel shared(deltah,myidx,a,desc_a)
    !
    block 
      integer(psb_ipk_) :: i,j,k,ii,ib,icoeff, ix,iy, ith,nth
      integer(psb_lpk_) :: glob_row
      integer(psb_lpk_), allocatable     :: irow(:),icol(:)
      real(psb_dpk_), allocatable :: val(:)
      real(psb_dpk_)    :: x,y, zt(nb)
#if defined(OPENMP)
      nth = omp_get_num_threads()
      ith = omp_get_thread_num()
#else
      nth = 1
      ith = 0
#endif
      allocate(val(20*nb),irow(20*nb),&
           &icol(20*nb),stat=info)
      if (info /= psb_success_ ) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        !goto 9999
      endif
      
      !$omp  do schedule(dynamic)
      !     
      do ii=1, nlr,nb
        if(info /= psb_success_) cycle
        ib = min(nb,nlr-ii+1)
        icoeff = 1
        do k=1,ib
          i=ii+k-1
          ! local matrix pointer
          glob_row=myidx(i)
          ! compute gridpoint coordinates
          call idx2ijk(ix,iy,glob_row,idim,idim)
          ! x, y coordinates
          x = (ix-1)*deltah
          y = (iy-1)*deltah
          
          zt(k) = f_(x,y)
          ! internal point: build discretization
          !
          !  term depending on   (x-1,y)
          !
          val(icoeff) = -a1(x,y)/sqdeltah-b1(x,y)/deltah2
          if (ix == 1) then
            zt(k) = g(dzero,y)*(-val(icoeff)) + zt(k)
          else
            call ijk2idx(icol(icoeff),ix-1,iy,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y-1)
          val(icoeff)  = -a2(x,y)/sqdeltah-b2(x,y)/deltah2
          if (iy == 1) then
            zt(k) = g(x,dzero)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy-1,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          
          !  term depending on     (x,y)
          val(icoeff)=(2*done)*(a1(x,y) + a2(x,y))/sqdeltah + c(x,y)
          call ijk2idx(icol(icoeff),ix,iy,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
          !  term depending on     (x,y+1)
          val(icoeff)=-a2(x,y)/sqdeltah+b2(x,y)/deltah2
          if (iy == idim) then
            zt(k) = g(x,done)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy+1,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x+1,y)
          val(icoeff)=-a1(x,y)/sqdeltah+b1(x,y)/deltah2
          if (ix==idim) then
            zt(k) = g(done,y)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix+1,iy,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          
        end do
#if defined(OPENMP)
!!$        write(0,*) omp_get_thread_num(),' Check insertion ',&
!!$             & irow(1:icoeff-1),':',icol(1:icoeff-1)
#endif
        call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
        if(info /= psb_success_) cycle
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
        if(info /= psb_success_) cycle
        zt(:)=dzero
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
        if(info /= psb_success_) cycle
      end do
      !$omp end do
      deallocate(val,irow,icol)
    end block
    !$omp end parallel
    
    
    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info,mold=imold)
    tcdasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    if (info == psb_success_) then
      if (present(amold)) then
        call psb_spasb(a,desc_a,info,mold=amold)
      else
        call psb_spasb(a,desc_a,info,afmt=afmt)
      end if
    end if
    call psb_barrier(ctxt)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (info == psb_success_) call psb_geasb(xv,desc_a,info,mold=vmold)
    if (info == psb_success_) call psb_geasb(bv,desc_a,info,mold=vmold)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    tasb = psb_wtime()-t1
    call psb_barrier(ctxt)
    ttot = psb_wtime() - t0

    call psb_amx(ctxt,talc)
    call psb_amx(ctxt,tgen)
    call psb_amx(ctxt,tasb)
    call psb_amx(ctxt,ttot)
    if(iam == psb_root_) then
      tmpfmt = a%get_fmt()
      write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")')&
           &   tmpfmt
      write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
      write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
      write(psb_out_unit,'("-desc asbly  time : ",es12.5)') tcdasb
      write(psb_out_unit,'("- mat asbly  time : ",es12.5)') tasb
      write(psb_out_unit,'("-total       time : ",es12.5)') ttot

    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_d_gen_pde2d

end module psb_d_pde2d_mod