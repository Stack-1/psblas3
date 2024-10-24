!    
!                Parallel Sparse BLAS  GPU plugin
!      (C) Copyright 2013
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
!    
! File: spdegenmv.f90
!
! Program: pdegenmv
! This sample program measures the performance of the matrix-vector product.
! The matrix is generated in the same way as for the pdegen test case of
! the main PSBLAS library.
!
!
module psb_s_pde3d_mod

 
    use psb_base_mod, only : psb_spk_, psb_ipk_, psb_lpk_, psb_desc_type,&
         &  psb_sspmat_type, psb_s_vect_type, szero,&
         &  psb_s_base_sparse_mat, psb_s_base_vect_type, &
         &  psb_i_base_vect_type, psb_l_base_vect_type
  
    interface 
      function s_func_3d(x,y,z) result(val)
        import :: psb_spk_
        real(psb_spk_), intent(in) :: x,y,z
        real(psb_spk_) :: val
      end function s_func_3d
    end interface 
  
    
  contains
  
    function s_null_func_3d(x,y,z) result(val)
  
      real(psb_spk_), intent(in) :: x,y,z
      real(psb_spk_) :: val
      
      val = szero
  
    end function s_null_func_3d
    !
    ! functions parametrizing the differential equation 
    !  
    function b1(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) :: b1
      real(psb_spk_), intent(in) :: x,y,z
      b1=sone/sqrt((3*sone))
    end function b1
    function b2(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  b2
      real(psb_spk_), intent(in) :: x,y,z
      b2=sone/sqrt((3*sone))
    end function b2
    function b3(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  b3
      real(psb_spk_), intent(in) :: x,y,z      
      b3=sone/sqrt((3*sone))
    end function b3
    function c(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  c
      real(psb_spk_), intent(in) :: x,y,z      
      c=szero
    end function c
    function a1(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  a1   
      real(psb_spk_), intent(in) :: x,y,z
      a1=sone/80
    end function a1
    function a2(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  a2
      real(psb_spk_), intent(in) :: x,y,z
      a2=sone/80
    end function a2
    function a3(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  a3
      real(psb_spk_), intent(in) :: x,y,z
      a3=sone/80
    end function a3
    function g(x,y,z)
      use psb_base_mod, only : psb_spk_, sone, szero
      implicit none 
      real(psb_spk_) ::  g
      real(psb_spk_), intent(in) :: x,y,z
      g = szero
      if (x == sone) then
        g = sone
      else if (x == szero) then 
        g = exp(y**2-z**2)
      end if
    end function g
  
    
    !
    !  subroutine to allocate and fill in the coefficient matrix and
    !  the rhs. 
    !
    subroutine psb_s_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,info,&
         & f,amold,vmold,imold,partition,nrl,iv)
      use psb_base_mod
      use psb_util_mod
      !
      !   Discretizes the partial differential equation
      ! 
      !   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
      ! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
      !      dxdx     dydy       dzdz        dx       dy         dz   
      !
      ! with Dirichlet boundary conditions
      !   u = g 
      !
      !  on the unit cube  0<=x,y,z<=1.
      !
      !
      ! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
      !
      implicit none
      integer(psb_ipk_)     :: idim
      type(psb_sspmat_type) :: a
      type(psb_s_vect_type) :: xv,bv
      type(psb_desc_type)   :: desc_a
      type(psb_ctxt_type) :: ctxt
      integer(psb_ipk_)     :: info
      character(len=*)      :: afmt
      procedure(s_func_3d), optional :: f
      class(psb_s_base_sparse_mat), optional :: amold
      class(psb_s_base_vect_type), optional :: vmold 
      class(psb_i_base_vect_type), optional :: imold
      integer(psb_ipk_), optional :: partition, nrl,iv(:)
  
      ! Local variables.
  
      integer(psb_ipk_), parameter :: nb=20
      type(psb_s_csc_sparse_mat)  :: acsc
      type(psb_s_coo_sparse_mat)  :: acoo
      type(psb_s_csr_sparse_mat)  :: acsr
      real(psb_spk_)           :: zt(nb),x,y,z
      integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_
      integer(psb_lpk_) :: m,n,glob_row,nt
      integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
      ! For 3D partition
      ! Note: integer control variables going directly into an MPI call
      ! must be 4 bytes, i.e. psb_mpk_
      integer(psb_mpk_) :: npdims(3), npp, minfo
      integer(psb_ipk_) :: npx,npy,npz, iamx,iamy,iamz,mynx,myny,mynz
      integer(psb_ipk_), allocatable :: bndx(:),bndy(:),bndz(:)
      ! Process grid
      integer(psb_ipk_) :: np, iam
      integer(psb_ipk_) :: icoeff
      integer(psb_lpk_), allocatable     :: irow(:),icol(:),myidx(:)
      real(psb_spk_), allocatable :: val(:)
      ! deltah dimension of each grid cell
      ! deltat discretization time
      real(psb_spk_)            :: deltah, sqdeltah, deltah2
      real(psb_spk_), parameter :: rhs=szero,one=sone,zero=szero
      real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
      integer(psb_ipk_) :: err_act
      procedure(s_func_3d), pointer :: f_
      character(len=20)  :: name, ch_err,tmpfmt
  
      info = psb_success_
      name = 'create_matrix'
      call psb_erractionsave(err_act)
  
      call psb_info(ctxt, iam, np)
  
  
      if (present(f)) then 
        f_ => f
      else
        f_ => s_null_func_3d
      end if
  
      deltah   = sone/(idim+2)
      sqdeltah = deltah*deltah
      deltah2  = (2*sone)* deltah
  
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
      
      m   = (1_psb_lpk_*idim)*idim*idim
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
        ! A 3-dimensional partition
  
        ! A nifty MPI function will split the process list
        npdims = 0
        call mpi_dims_create(np,3,npdims,info)
        npx = npdims(1)
        npy = npdims(2)
        npz = npdims(3)
  
        allocate(bndx(0:npx),bndy(0:npy),bndz(0:npz))
        ! We can reuse idx2ijk for process indices as well. 
        call idx2ijk(iamx,iamy,iamz,iam,npx,npy,npz,base=0)
        ! Now let's split the 3D cube in hexahedra
        call dist1Didx(bndx,idim,npx)
        mynx = bndx(iamx+1)-bndx(iamx)
        call dist1Didx(bndy,idim,npy)
        myny = bndy(iamy+1)-bndy(iamy)
        call dist1Didx(bndz,idim,npz)
        mynz = bndz(iamz+1)-bndz(iamz)
  
        ! How many indices do I own? 
        nlr = mynx*myny*mynz
        allocate(myidx(nlr))
        ! Now, let's generate the list of indices I own
        nr = 0
        do i=bndx(iamx),bndx(iamx+1)-1
          do j=bndy(iamy),bndy(iamy+1)-1
            do k=bndz(iamz),bndz(iamz+1)-1
              nr = nr + 1
              call ijk2idx(myidx(nr),i,j,k,idim,idim,idim)
            end do
          end do
        end do
        if (nr /= nlr) then
          write(psb_err_unit,*) iam,iamx,iamy,iamz, 'Initialization error: NR vs NLR ',&
               & nr,nlr,mynx,myny,mynz
          info = -1
          call psb_barrier(ctxt)
          call psb_abort(ctxt)
        end if
  
        !
        ! Third example  of use of CDALL: specify for each process
        ! the set of global indices it owns.
        ! 
        call psb_cdall(ctxt,desc_a,info,vl=myidx)
        
      case default
        write(psb_err_unit,*) iam, 'Initialization error: should not get here'
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return
      end select
  
      
      if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz,&
           & dupl=psb_dupl_err_)
      ! define  rhs from boundary conditions; also build initial guess 
      if (info == psb_success_) call psb_geall(xv,desc_a,info)
      if (info == psb_success_) call psb_geall(bv,desc_a,info)
  
      call psb_barrier(ctxt)
      talc = psb_wtime()-t0
  
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='allocation rout.'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
  
      ! we build an auxiliary matrix consisting of one row at a
      ! time; just a small matrix. might be extended to generate 
      ! a bunch of rows per call. 
      ! 
      allocate(val(20*nb),irow(20*nb),&
           &icol(20*nb),stat=info)
      if (info /= psb_success_ ) then 
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999
      endif
  
  
      ! loop over rows belonging to current process in a block
      ! distribution.
  
      call psb_barrier(ctxt)
      t1 = psb_wtime()
      do ii=1, nlr,nb
        ib = min(nb,nlr-ii+1) 
        icoeff = 1
        do k=1,ib
          i=ii+k-1
          ! local matrix pointer 
          glob_row=myidx(i)
          ! compute gridpoint coordinates
          call idx2ijk(ix,iy,iz,glob_row,idim,idim,idim)
          ! x, y, z coordinates
          x = (ix-1)*deltah
          y = (iy-1)*deltah
          z = (iz-1)*deltah
          zt(k) = f_(x,y,z)
          ! internal point: build discretization
          !   
          !  term depending on   (x-1,y,z)
          !
          val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
          if (ix == 1) then 
            zt(k) = g(szero,y,z)*(-val(icoeff)) + zt(k)
          else
            call ijk2idx(icol(icoeff),ix-1,iy,iz,idim,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y-1,z)
          val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
          if (iy == 1) then 
            zt(k) = g(x,szero,z)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy-1,iz,idim,idim,idim)          
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y,z-1)
          val(icoeff)=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
          if (iz == 1) then 
            zt(k) = g(x,y,szero)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy,iz-1,idim,idim,idim)          
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
  
          !  term depending on     (x,y,z)
          val(icoeff)=(2*sone)*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah &
               & + c(x,y,z)
          call ijk2idx(icol(icoeff),ix,iy,iz,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1                  
          !  term depending on     (x,y,z+1)
          val(icoeff)=-a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
          if (iz == idim) then 
            zt(k) = g(x,y,sone)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy,iz+1,idim,idim,idim)          
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y+1,z)
          val(icoeff)=-a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
          if (iy == idim) then 
            zt(k) = g(x,sone,z)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy+1,iz,idim,idim,idim)          
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x+1,y,z)
          val(icoeff)=-a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
          if (ix==idim) then 
            zt(k) = g(sone,y,z)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix+1,iy,iz,idim,idim,idim)          
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
  
        end do
        call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
        if(info /= psb_success_) exit
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
        if(info /= psb_success_) exit
        zt(:)=szero
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
        if(info /= psb_success_) exit
      end do
  
      tgen = psb_wtime()-t1
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='insert rout.'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
  
      deallocate(val,irow,icol)
  
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
    end subroutine psb_s_gen_pde3d
  
  
  end module psb_s_pde3d_mod