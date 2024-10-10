module psb_ds_gen_mod
  use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_desc_type,&
         &  psb_dspmat_type, psb_d_vect_type, dzero,&
         &  psb_d_base_sparse_mat, psb_d_base_vect_type, psb_i_base_vect_type
  
  contains
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
    subroutine psb_ds_gen_matrix(ctxt,idim,a,b,xv,a_lower_precision, b_lower_precision, x_lower_precision, desc_a,afmt,info)
      use psb_base_mod
      use psb_util_mod
  
      implicit none
      
      integer(psb_ipk_)     :: idim
      type(psb_dspmat_type) :: a
      type(psb_sspmat_type) :: a_lower_precision
      type(psb_d_vect_type) :: xv,b
      type(psb_s_vect_type) :: x_lower_precision, b_lower_precision
      type(psb_desc_type)   :: desc_a
      type(psb_ctxt_type)   :: ctxt
      integer(psb_ipk_)     :: info
      character(len=*)      :: afmt
  
      ! Local variables.
  
      integer(psb_ipk_), parameter :: nb=20
      real(psb_dpk_)           :: zt(nb),x,y,z
      integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_, mysz
      integer(psb_lpk_) :: m,n,glob_row,nt
      integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
      ! For 2D partition
      ! Note: integer control variables going directly into an MPI call
      ! must be 4 bytes, i.e. psb_mpk_
      integer(psb_mpk_) :: npdims(2), npp, minfo
      integer(psb_ipk_) :: npx,npy,my_rankx,my_ranky,mynx,myny
      integer(psb_ipk_), allocatable :: bndx(:),bndy(:)
      ! Process grid
      integer(psb_ipk_) :: np, my_rank
      integer(psb_ipk_) :: icoeff
      integer(psb_lpk_), allocatable     :: myidx(:)
      ! deltah dimension of each grid cell
      ! deltat discretization time
      real(psb_dpk_)            :: deltah, sqdeltah, deltah2
      real(psb_dpk_), parameter :: rhs=dzero,one=done,zero=dzero
      real(psb_dpk_)            :: initial_timer, temporary_timer, allocation_time, generation_time, building_time, total_time
      integer(psb_ipk_) :: err_act
      character(len=20)  :: name, ch_err,tmpfmt
      logical :: debug

      info = psb_success_
      name = 'create_matrix'
      call psb_erractionsave(err_act)
  
      debug = .false.

      call psb_info(ctxt, my_rank, np)
  

  
      deltah   = done/(idim+1)
      sqdeltah = deltah*deltah
      deltah2  = (2*done)* deltah
  
  
      ! initialize array descriptor and sparse matrix storage. provide an
      ! estimate of the number of non zeroes
  
      m   = (1_psb_lpk_)*idim*idim
      n   = m
      nnz = ((n*7)/(np))

      if(my_rank == psb_root_) then 
        print '("=================================================================================")'
        print '("=                        MATRIX GENERATION SECTION                              =")'
        print '("=================================================================================")'
        if(debug .eqv. .true.) print '("[INFO] Step 0/6")'
        print '("[INFO] Generating Matrix ",i0,"x",i0,"...")' , n , n
        print '(" ")'
      end if

      initial_timer = psb_wtime()

      !
      ! Using a simple BLOCK ROWS distribution.
      !
      nt = (m+np-1)/np
      nr = max(0,min(nt,m-(my_rank*nt)))

      nt = nr
      call psb_sum(ctxt,nt)
      if (nt /= m) then
        write(psb_err_unit,*) '[ERROR] Process ', my_rank, ': Initialization error ',nr,nt,m
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
      
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='allocating matrix descriptor' 
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Step 1/6")'
        print '("[INFO] Descriptor succesfully allocated")'
        print '(" ")'
      end if
  
      if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Step 2/6")'
        print '("[INFO] Matrix in double precision succesfully allocated")'
      end if

      if (info == psb_success_) call psb_geall(xv,desc_a,info)
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] x vector in double precision succesfully allocated")'
      end if

      if (info == psb_success_) call psb_geall(b,desc_a,info)
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] b vector in double precision succesfully allocated")'
        print '(" ")'
      end if

      call psb_barrier(ctxt)
      
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='allocating double precision data structures' 
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      ! Allocates lower precision data structures
      if (info == psb_success_) call psb_spall(a_lower_precision,desc_a,info,nnz=nnz)
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Step 3/6")'
        print '("[INFO] Matrix in single precision succesfully allocated")'
      end if

      ! define  rhs from boundary conditions; also build initial guess
      if (info == psb_success_) call psb_geall(x_lower_precision,desc_a,info)
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] x vector in single precision succesfully allocated")'
      end if

      if (info == psb_success_) call psb_geall(b_lower_precision,desc_a,info)
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] x vector in single precision succesfully allocated")'
        print '(" ")'
      end if
      

  
      call psb_barrier(ctxt)      
      
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='allocating single precision data structures' 
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
   
      allocation_time = psb_wtime() - initial_timer

  
  
      call psb_barrier(ctxt)

      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Step 4/6")'
        print '("[INFO] Starting coefficents generation... ")'
      end if


      temporary_timer = psb_wtime()

      !$omp parallel shared(deltah,myidx,a,desc_a)
      !
      block 
        integer(psb_ipk_) :: i,j,k,ii,ib,icoeff, ix,iy, ith,nth
        integer(psb_lpk_) :: glob_row
        integer(psb_lpk_), allocatable     :: irow(:),icol(:)
        real(psb_dpk_), allocatable :: val(:)
        real(psb_dpk_)    :: x,y, zt(nb)
        
        real(psb_spk_), allocatable :: val_lower(:)
        real(psb_spk_)              :: x_lower,y_lower, zt_lower(nb)
        
        nth = 1
        ith = 0
        allocate(val(20*nb),irow(20*nb),&
              &icol(20*nb),stat=info)
        
        if (info /= psb_success_ ) then
          info=psb_err_alloc_dealloc_
          ch_err='allocating val vector' 
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
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
              zt(k) = dzero
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


            val_lower = val 
            zt_lower = zt

    
            call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
            if(info /= psb_success_) cycle
            call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),b,desc_a,info)
            if(info /= psb_success_) cycle
            zt(:) = dzero
            call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
            if(info /= psb_success_) cycle


    
            call psb_spins(icoeff-1,irow,icol,val_lower,a_lower_precision,desc_a,info)
            if(info /= psb_success_) cycle
            call psb_geins(ib,myidx(ii:ii+ib-1),zt_lower(1:ib),b_lower_precision,desc_a,info)
            if(info /= psb_success_) cycle
            zt_lower(:) = szero
            call psb_geins(ib,myidx(ii:ii+ib-1),zt_lower(1:ib),x_lower_precision,desc_a,info)
            if(info /= psb_success_) cycle



        end do
        !$omp end do
        deallocate(val,irow,icol)
      end block
      !$omp end parallel

      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Coefficents generated succesfully")'
        print '(" ")'
        print '("[INFO] Step 5/6")'
        print '("[INFO] Coefficents inserted in double precision data structures succesfully")'
        print '("[INFO] Coefficents inserted in single precision data structures succesfully")'
        print '(" ")'
      end if

      call psb_barrier(ctxt)

      
      if(info /= psb_success_) then
        info = psb_err_from_subroutine_
        ch_err = 'insert rout.'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      generation_time = psb_wtime() - temporary_timer

      temporary_timer = psb_wtime()

      ! Build all entities
      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Step 6/6")'
        print '("[INFO] Building generated data structures...")'
      end if

      call psb_cdasb(desc_a,info)
      if(info == psb_success_) call psb_spasb(a,desc_a,info)
      if(info == psb_success_) call psb_spasb(a_lower_precision,desc_a,info)
      
      if(info /= psb_success_) then
        info = psb_err_from_subroutine_
        ch_err = 'build data structures'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      building_time = psb_wtime() - temporary_timer

      if( (my_rank == psb_root_).and.(debug .eqv. .true.) ) then 
        print '("[INFO] Data structures generated")'
        print '(" ")'
      end if


      total_time = psb_wtime() - initial_timer
    
      ! Needed to pick the timer of the slowest process
      call psb_amx(ctxt,allocation_time)
      call psb_amx(ctxt,generation_time)
      call psb_amx(ctxt,building_time)
      call psb_amx(ctxt,total_time)
        
      if(my_rank == psb_root_) then
        tmpfmt = a%get_fmt()
        print '("[INFO] The matrix has been generated and assembled in ",a3," format.")' , tmpfmt
        print '("       - allocation time             : ",es12.5, " s")' , allocation_time
        print '("       - coefficent generation time  : ",es12.5, " s")' , generation_time
        print '("       - building entities time      : ",es12.5, " s")' , building_time
        print '("       - total time                  : ",es12.5, " s")' , total_time
        print '(" ")'
      end if 
      
      call psb_erractionrestore(err_act)
      
      return
  
  9999 call psb_error_handler(ctxt,err_act)
  
      return
    end subroutine psb_ds_gen_matrix
  
  end module psb_ds_gen_mod