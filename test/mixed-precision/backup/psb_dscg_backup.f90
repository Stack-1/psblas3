!   
!   Parallel Sparse BLAS  version 3.5
!   (C) Copyright 2006-2018
!       Salvatore Filippone    
!       Alfredo Buttari      
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
! File:  psb_dcg.f90
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   C                                                                      C
!   C  References:                                                         C
!   C          [1] Duff, I., Marrone, M., Radicati, G., and Vittoli, C.    C
!   C              Level 3 basic linear algebra subprograms for sparse     C
!   C              matrices: a user level interface                        C
!   C              ACM Trans. Math. Softw., 23(3), 379-401, 1997.          C
!   C                                                                      C
!   C                                                                      C
!   C         [2]  S. Filippone, M. Colajanni                              C
!   C              PSBLAS: A library for parallel linear algebra           C
!   C              computation on sparse matrices                          C
!   C              ACM Trans. on Math. Softw., 26(4), 527-550, Dec. 2000.  C
!   C                                                                      C
!   C         [3] M. Arioli, I. Duff, M. Ruiz                              C
!   C             Stopping criteria for iterative solvers                  C
!   C             SIAM J. Matrix Anal. Appl., Vol. 13, pp. 138-144, 1992   C
!   C                                                                      C
!   C                                                                      C
!   C         [4] R. Barrett et al                                         C
!   C             Templates for the solution of linear systems             C
!   C             SIAM, 1993                                          
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_dscg.f90
! 
! Author: Staccone Simone
!
! Subroutine: psb_dscg
!    This subroutine implements the Conjugate Gradient method in mixed precision.
!
!
! Arguments:
!
!    matrix      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)       Input: preconditioner
!    b(:)   -  real                    Input: vector containing the
!                                         right hand side B
!    x(:)   -  real                    Input/Output: vector containing the
!                                         initial guess and final solution X.
!    eps    -  real                       Input: Stopping tolerance; the iteration is
!                                         stopped when the error estimate |err| <= eps
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                   Output: Return code
!
!    itmax  -  integer(optional)          Input: maximum number of iterations to be
!                                         performed.
!    iter   -  integer(optional)          Output: how many iterations have been
!                                         performed.
!                                         performed.
!    err    -  real   (optional)          Output: error estimate on exit. If the
!                                         denominator of the estimate is exactly
!                                         0, it is changed into 1. 
!    itrace -  integer(optional)          Input: print an informational message
!                                         with the error estimate every itrace
!                                         iterations
!    istop  -  integer(optional)          Input: stopping criterion, or how
!                                         to estimate the error. 
!                                         1: err =  |r|/(|a||x|+|b|);  here the iteration is
!                                            stopped when  |r| <= eps * (|a||x|+|b|)
!                                         2: err =  |r|/|b|; here the iteration is
!                                            stopped when  |r| <= eps * |b|
!                                         where r is the (preconditioned, recursive
!                                         estimate of) residual. 
! 


program psb_dscg 
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod


  implicit none
  
  type(psb_dspmat_type) :: matrix_A
  type(psb_sspmat_type) :: matrix_A_lower_precision

  type(psb_desc_type) :: desc_a, desc_a_lower_precision

  type(psb_ldspmat_type) :: global_A
  type(psb_lsspmat_type) :: global_A_lower_precision

  type(psb_ctxt_type) :: comm_ctxt

  type(psb_d_vect_type)    :: b_col, x_col, r_col, b_col_lower_precision

  ! preconditioner metadata
  type(psb_dprec_type)  :: preconditioner



  real(psb_dpk_), pointer  :: b_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob_lower_precision(:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)


  character(len=40) :: krilov_method, preconditioner_string, partition, matrix_file, rhs_file
  character(len=256) :: specs_file_path
  character(len=5) :: matrix_format
  integer(psb_ipk_) :: iret, istopc,info, max_iteretions,itrace_,irst, my_rank, np, err_act, & 
  & rows, cols, nnz, rows_per_process, nlr, total_rows_check, nr
  real(psb_dpk_) :: eps, condition_number, err
  integer(psb_ipk_) :: iteration_number
  character(len=20) :: name
  real(psb_dpk_) :: t1, t2 ! timers variables
  integer(psb_lpk_), allocatable :: my_index(:)
  integer(psb_ipk_), parameter :: iunit=12
  real(psb_dpk_), allocatable, target ::  aux_b(:,:)
  real(psb_dpk_), allocatable, target ::  aux_b_lower_precision(:,:)


  integer :: i

  name = 'mixed_psb_dscg'
  info = psb_success_
  call psb_erractionsave(err_act)


  call psb_init(comm_ctxt)
  call psb_info(comm_ctxt, my_rank, np)


  if( (my_rank < 0).or.(my_rank >= np) ) then
    write(*,*) '[ERROR] Unexpected wrong rank value'
    call psb_exit(comm_ctxt)
    stop
  end if

  if( my_rank == psb_root_) then
    print '("[INFO] Process ",I0 ,": PSBLAS environment started correctly")', my_rank
    print '("[INFO] Process ", I0, ": Reading parameters from file... "  )', my_rank
  end if


  ! Collect parameters from configuration file
  call get_parms(comm_ctxt,specs_file_path,krilov_method,preconditioner_string,partition, matrix_file,&
  &  matrix_format,rhs_file,istopc,max_iteretions,itrace_,irst,eps)

  if( my_rank == psb_root_) then
    print '("[INFO] Process ",I0 ,": All parameters have been aquired correctly")', my_rank
    print '("")' 
    print '("=================================================================================")'
  end if

  call psb_barrier(comm_ctxt)
  if(my_rank == psb_root_) then 
    t1 = psb_wtime() 
    t2 = 0
  endif

  ! Read data of the linera system Ax = b
  if (my_rank == psb_root_) then
    ! We want to aquire the same matrix in different precision
    print '("[INFO] Process ", I0, ": Reading double precision matrix... "  )', my_rank
    call mm_mat_read(global_A,info,iunit=iunit,filename=matrix_file)
    t2 = t2 + ( psb_wtime() - t1 )
    call check_correct(info,name,err_act)
    print '("[INFO] Process ", I0, ": Double precision matrix succesfully read")', my_rank


    print '("[INFO] Process ", I0, ": Reading single precision matrix... "  )', my_rank
    t1 = psb_wtime()
    call mm_mat_read(global_A_lower_precision,info,iunit=iunit,filename=matrix_file)
    t2 = t2 + ( psb_wtime() - t1 )
    call check_correct(info,name,err_act)
    print '("[INFO] Process ", I0, ": Single precision matrix succesfully read")', my_rank

    ! Read rhs (b) vector from file
    print '("[INFO] Process ", I0, ": Reading double precision rhs vector... "  )', my_rank
    t1 = psb_wtime()
    call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
    t2 = t2 + ( psb_wtime() - t1 )
    call check_correct(info,name,err_act)
    print '("[INFO] Process ", I0, ": Double precision rhs vector succesfully read")', my_rank

  ! Read rhs (b) vector from file
    print '("[INFO] Process ", I0, ": Reading single precision rhs vector... "  )', my_rank
    t1 = psb_wtime()
    call mm_array_read(aux_b_lower_precision,info,iunit=iunit,filename=rhs_file)
    t2 = t2 + ( psb_wtime() - t1 )
    call check_correct(info,name,err_act)
    print '("[INFO] Process ", I0, ": Single single rhs vector succesfully read")', my_rank


    ! Broadcast rhs vector to all process, making the global var point to it
    b_col_glob =>aux_b(:,1)
    b_col_glob_lower_precision => aux_b_lower_precision(:,1)

    ! Debug information to check how the matrix is composed
    rows = global_A%get_nrows()
    cols = global_A%get_ncols()
    nnz = global_A%get_nzeros() 
  end if

  call psb_barrier(comm_ctxt)
  if(my_rank == psb_root_) then 
    print '("[INFO] Process ", I0, ": Global matrix is " , I0 , " x " ,I0)', my_rank, rows, cols
    print '("[INFO] Process ", I0, ": Time to read matrix and rhs vector is " , es12.5 )', my_rank, t2
    print '("")' 
    print '("=================================================================================")'
  endif

  ! Extra communication to check if everyone has the correct information about own rows and columns
  call psb_bcast(comm_ctxt,rows)
  call psb_bcast(comm_ctxt,cols)
  rows_per_process = rows / np
  nr = mod(rows,np)
  if( my_rank < nr) then
      rows_per_process = rows_per_process + 1
  end if
  total_rows_check = rows_per_process

  call psb_sum(comm_ctxt,total_rows_check)
  
  ! Check that the sum of the individual rows of each process gives back the total rows of the global matrix
  if (total_rows_check /= rows) then
    write(psb_err_unit,*) my_rank, 'Initialization error ',total_rows_check,rows_per_process,rows
    info = -1
    call psb_barrier(comm_ctxt)
    call psb_abort(comm_ctxt)
    return
  endif


  ! Distribute matrix A over process using BLOCK distribution
  if(my_rank == psb_root_) then 
    print '("[INFO] Process ", I0, ": Starting double precision matrix distribution over processes... "  )', my_rank
    t1 = psb_wtime() 
  endif
  call psb_matdist(global_A, matrix_A,comm_ctxt,desc_a,info,fmt=matrix_format,parts=part_block)
  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
  endif 
  call check_correct(info,name,err_act)
  if(my_rank == psb_root_) then 
    print '("[INFO] Process ", I0, ": Double precision matrix distribution ended correctly")', my_rank
    print '("[INFO] Process ", I0, ": Time to allocate and distribute double precision matrix is " , es12.5 )', my_rank, t2
    print '("")'  
    print '("=================================================================================")' 
  endif


  if(my_rank == psb_root_) then 
    print '("[INFO] Process ", I0, ": Starting single precision matrix distribution over processes... "  )', my_rank
    t1 = psb_wtime() 
  endif
  call psb_matdist(global_A_lower_precision, matrix_A_lower_precision,&
                  & comm_ctxt,desc_a_lower_precision,info,fmt=matrix_format,parts=part_block)

  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
  endif 
  call check_correct(info,name,err_act)
  if(my_rank == psb_root_) then 
    print '("[INFO] Process ", I0, ": Single precision matrix distribution ended correctly")', my_rank
    print '("[INFO] Process ", I0, ": Time to allocate and distribute single precision matrix is " , es12.5 )', my_rank, t2
    print '("")' 
    print '("=================================================================================")'
  endif

  ! Distribute vector RHS over processes
  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Starting double precision rhs distribution over processes... "  )', my_rank
    t1 = psb_wtime() 
  endif

  call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_)

  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
  endif 

  call check_correct(info,name,err_act)

  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
    print '("[INFO] Process ", I0, ": Double precision rhs vector distribution completed succesfully"  )', my_rank
    print '("[INFO] Process ", I0, ": Time to allocate and distribute vector rhs is " , es12.5 )', my_rank, t2
    print '("")' 
    print '("=================================================================================")'
  endif 


  ! Distribute vector RHS over processes
  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Starting single precision rhs distribution over processes... "  )', my_rank
    t1 = psb_wtime() 
  endif

  call psb_scatter(b_col_glob_lower_precision,b_col_lower_precision,desc_a_lower_precision,info,root=psb_root_)

  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
  endif 

  call check_correct(info,name,err_act)

  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
    print '("[INFO] Process ", I0, ": Single precision rhs vector distribution completed succesfully"  )', my_rank
    print '("[INFO] Process ", I0, ": Time to allocate and distribute vector rhs is " , es12.5 )', my_rank, t2
    print '("")' 
    print '("=================================================================================")'
  endif 


  ! Now that the matrix is distributed and the descriptor has been initialized
  ! we need to build our preconditioner
  if(preconditioner_string /= "NONE") then
    call preconditioner%init(comm_ctxt,preconditioner_string,info)

    ! building the preconditioner
    t1 = psb_wtime()
    call preconditioner%build(matrix_A,desc_a,info)
    t2 = psb_wtime()-t1
    if(my_rank == psb_root_) then
      print '("[INFO] Process ", I0, ": Time to build preconditioner is " , es12.5 )', my_rank, t2
    end if
    call check_correct(info,name,err_act)
  else
    if(my_rank == psb_root_) then
      print '("[INFO] Process ", I0, ": No preconditioner was selected,&
      & solving system without using a preconditioner " )', my_rank
    end if
  endif

  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Time to build preconditioner is " , es12.5 )', my_rank, t2
    print '("")' 
    print '("=================================================================================")'
  end if

  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Allocationg local x vector...")', my_rank
  end if

  call psb_geall(x_col,desc_a,info)
  call x_col%zero()
  call psb_geasb(x_col,desc_a,info)
  call check_correct(info,name,err_act)
  call psb_barrier(comm_ctxt)

  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Local vector x allocated correctly")', my_rank
    print '("")'
    print '("=================================================================================")'
  end if


  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Starting computation...")', my_rank
  end if

  

  call psb_dscg_impl(matrix_A,preconditioner,b_col,x_col,eps,desc_a,info, max_iteretions, &
                    & iteration_number, err, itrace=itrace_, istop=istopc, cond=condition_number)
  ! call check_correct(info,name,err_act)
  call psb_barrier(comm_ctxt)

  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Computation ended correctly")', my_rank
    print '("[INFO] Process ", I0, ": Size of solution " , I0)', my_rank, x_col%get_nrows()
    print '("")'
    print '("=================================================================================")'
  end if


  call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
  ! call check_correct(info,name,err_act)

  call psb_barrier(comm_ctxt)


  if(my_rank == psb_root_) then
    write(20,*) 'matrix: ',matrix_file
    write(20,*) 'computed solution on ',np,' processors.'
    write(20,*) 'iterations to convergence: ',iteration_number
    write(20,*) 'error estimate (infinity norm) on exit:', &
         & ' ||r||/(||a||||x||+||b||) = ',err
    write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
    do i=1,desc_a%get_global_rows()
      write(20,*) i,x_col_glob(i),b_col_glob(i)
    enddo
  endif


  if(my_rank == 0) then 
    print '("[INFO] Process ", I0, ": Computation ended correctly, check file to see the results")', my_rank
    print '("[INFO] Process ", I0, ": PSBLAS environment shutting down")', my_rank
  endif

  call psb_exit(comm_ctxt)

contains

! Devug sobroutine to check if info has risen error
subroutine check_correct(info, name, err_act)
  use psb_base_mod
  implicit none

  integer(psb_ipk_) :: info, err_act
  character(len=20) :: name

  if ( info /= psb_success_ ) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    call psb_error_handler(err_act)
  end if

end subroutine


subroutine  get_parms(comm_ctxt,specs_file_path,krilov_method,preconditioner,partition,&
  & matrix_file,matrix_format,rhs_file,istopc,max_iteretions,itrace,irst,eps)
  
  use psb_base_mod
  implicit none

  type(psb_ctxt_type) :: comm_ctxt
  character(len=40) :: krilov_method, matrix_file, preconditioner, rhs_file
  character(len=256) :: specs_file_path
  character(len=20) :: partition, filefmt
  integer(psb_ipk_) :: iret, istopc,max_iteretions,itrace,irst
  character(len=40) :: charbuf
  real(psb_dpk_) :: eps
  character(len=5)    :: matrix_format
  integer(psb_ipk_) :: np, my_rank, info
  integer(psb_ipk_) :: inparms(40), ip, inp_unit

  call psb_info(comm_ctxt,my_rank,np)

  if (my_rank == 0) then
    if (command_argument_count()>0) then
      call get_command_argument(1,specs_file_path)
      inp_unit = 30
      open(inp_unit,file=specs_file_path,action='read',iostat=info)
      
      if (info /= 0) then
        print '("[INFO] Process ", I0, ": Could not open file ",a,"for input ")', my_rank, trim(specs_file_path)
        call psb_abort(comm_ctxt)
        stop
      else
        print '("[INFO] Process ", I0, ": Opened file ",a,"for input ")', my_rank, trim(specs_file_path)
      end if
    
    else
      inp_unit=inp_unit
    end if
 
    ! Read Input Parameters
    read(inp_unit,*) ip

    if (ip >= 5) then
      read(inp_unit,*) matrix_file
      read(inp_unit,*) rhs_file
      read(inp_unit,*) filefmt
      read(inp_unit,*) krilov_method
      read(inp_unit,*) preconditioner
      read(inp_unit,*) matrix_format
      read(inp_unit,*) partition



      call psb_bcast(comm_ctxt,matrix_file)
      call psb_bcast(comm_ctxt,rhs_file)
      call psb_bcast(comm_ctxt,filefmt)
      call psb_bcast(comm_ctxt,krilov_method)
      call psb_bcast(comm_ctxt,preconditioner)
      call psb_bcast(comm_ctxt,matrix_format)
      call psb_bcast(comm_ctxt,partition)
      
      
      ! Read stop criterion parameters or set default
      if (ip >= 7) then
        read(inp_unit,*) istopc
      else
        istopc=1        
      endif

      ! Read maximum iterations parameters or set default
      if (ip >= 8) then
        read(inp_unit,*) max_iteretions
      else
        max_iteretions=500
      endif
      
      ! Read verbose parameters or set default
      if (ip >= 9) then
        read(inp_unit,*) itrace
      else
        itrace=-1
      endif
      
      ! Read restart parameters or set default
      if (ip >= 10) then
        read(inp_unit,*) irst
      else
        irst  = 1
      endif
      
      ! Read error parameters or set default
      if (ip >= 11) then
        read(inp_unit,*) eps
      else
        eps=1.d-6
      endif

      inparms(1) = istopc
      inparms(2) = max_iteretions
      inparms(3) = itrace
      inparms(4) = irst
      call psb_bcast(comm_ctxt,inparms(1:4))
      call psb_bcast(comm_ctxt,eps)

      write(psb_out_unit,'("Solving matrix       : ",a)')  matrix_file      
      write(psb_out_unit,'("Number of processors : ",i3)') np
      write(psb_out_unit,'("Data distribution    : ",a)') partition
      write(psb_out_unit,'("Iterative method     : ",a)')  krilov_method
      write(psb_out_unit,'("Preconditioner       : ",a)')  preconditioner
      write(psb_out_unit,'("Restart parameter    : ",i2)') irst
      write(psb_out_unit,'("Storage format       : ",a)')  matrix_format
      write(psb_out_unit,'(" ")')
    else
      write(psb_err_unit,*) 'Wrong format for input file'
      call psb_abort(comm_ctxt)
      stop 1
    end if
    if (inp_unit /= psb_inp_unit) then
      close(inp_unit)
    end if
  else
    ! Receive Parameters
    call psb_bcast(comm_ctxt,matrix_file)
    call psb_bcast(comm_ctxt,rhs_file)
    call psb_bcast(comm_ctxt,filefmt)
    call psb_bcast(comm_ctxt,krilov_method)
    call psb_bcast(comm_ctxt,preconditioner)
    call psb_bcast(comm_ctxt,matrix_format)
    call psb_bcast(comm_ctxt,partition)

    call psb_bcast(comm_ctxt,inparms(1:4))
    istopc =  inparms(1) 
    max_iteretions  =  inparms(2) 
    itrace =  inparms(3) 
    irst   =  inparms(4) 
    call psb_bcast(comm_ctxt,eps)

  end if

end subroutine get_parms


!
!
!
!

subroutine psb_dscg_impl(a,prec,b,x,eps,desc_a,info,&
  & itmax,iter,err,itrace,istop,cond)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_krylov_conv_mod
  use psb_krylov_mod
  
  implicit none
  type(psb_dspmat_type), intent(in)              :: a
  Type(psb_desc_type), Intent(in)                :: desc_a
  class(psb_dprec_type), intent(inout)           :: prec
  type(psb_d_vect_type), Intent(inout)           :: b
  type(psb_d_vect_type), Intent(inout)           :: x
  Real(psb_dpk_), Intent(in)                     :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_dpk_), Optional, Intent(out)          :: err,cond
  ! =   Local data
  real(psb_dpk_), allocatable, target            :: aux(:),td(:),tu(:),eig(:),ewrk(:)
  integer(psb_mpk_), allocatable                 :: ibl(:), ispl(:), iwrk(:)
  type(psb_d_vect_type), allocatable, target     :: wwrk(:)
  type(psb_d_vect_type), pointer                 :: q, p, r, z, w
  real(psb_dpk_)                                 :: alpha, beta, rho, rho_old, sigma,alpha_old,beta_old
  integer(psb_ipk_)                              :: itmax_, istop_, naux, it, itx, itrace_,&
                                                    &  n_col, n_row,err_act, ieg,nspl, istebz
  integer(psb_lpk_)                              :: mglob
  integer(psb_ipk_)                              :: debug_level, debug_unit
  type(psb_ctxt_type)                            :: ctxt
  integer(psb_ipk_)                              :: np, me
  real(psb_dpk_)                                 :: derr  
  type(psb_itconv_type)                          :: stopdat
  logical                                        :: do_cond
  character(len=20)                              :: name
  character(len=*), parameter                    :: methdname='CG'
  logical, parameter                             :: do_timings=.true.
  integer(psb_ipk_), save                        :: cg_vect=-1, cg_mv=-1, cg_prec=-1
  
  info = psb_success_
  name = 'mixed_psb_dscg'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  
  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (.not.allocated(b%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if ((do_timings).and.(cg_vect==-1))       &
      & cg_vect = psb_get_timer_idx("CG: vector ops ")
  if ((do_timings).and.(cg_mv==-1))       &
      & cg_mv = psb_get_timer_idx("CG: MV product")
  if ((do_timings).and.(cg_prec==-1))       &
      & cg_prec = psb_get_timer_idx("CG: preconditioner")
  
  
  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()
  
  
  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif
  

  call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
  if (info == psb_success_) then 
    call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
  endif
  if(info /= psb_success_) then
   info=psb_err_from_subroutine_    
   call psb_errpush(info,name,a_err='psb_chkvect on X/B')
   goto 9999
  end if

  naux=4*n_col
  allocate(aux(naux), stat=info)
  if (info == psb_success_) call psb_geall(wwrk,desc_a,info,n=5_psb_ipk_)
  if (info == psb_success_) call psb_geasb(wwrk,desc_a,info,mold=x%v,scratch=.true.)  
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if
  
  p  => wwrk(1)
  q  => wwrk(2)
  r  => wwrk(3)
  z  => wwrk(4) 
  w  => wwrk(5)
  
  
  if (present(itmax)) then 
    itmax_ = itmax
  else
    itmax_ = 1000
  endif
  
  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = 0
  end if
  
  do_cond=present(cond)
  if (do_cond) then 
    istebz = 0
    allocate(td(itmax_),tu(itmax_), eig(itmax_),&
        & ibl(itmax_),ispl(itmax_),iwrk(3*itmax_),ewrk(4*itmax_),&
        & stat=info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
  end if
  itx=0
  alpha = dzero
  

  
  restart: do 
  ! =   
  ! =    r0 = b-Ax0
  ! =   
    if (do_timings) call psb_tic(cg_vect)
    if (itx>= itmax_) exit restart 
    it = 0
    call psb_geaxpby(done,b,dzero,r,desc_a,info)
    if (do_timings) call psb_toc(cg_vect)
    if (do_timings) call psb_tic(cg_mv)
    if (info == psb_success_) call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

    if (do_timings) call psb_toc(cg_mv)
   
    if (do_timings) call psb_tic(cg_vect)
    rho = dzero
   
    call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_non_,name)
      goto 9999
    endif
    if (do_timings) call psb_toc(cg_vect)

    ! Clone of the matrix in lower precision
    ! call matrix_single_clone(a, matrix_A_lower_precision)
  
    iteration:  do 
      ! This should be done all in single precision

      it   = it + 1
      itx = itx + 1
      if (do_timings) call psb_tic(cg_prec)
     
      call prec%apply(r,z,desc_a,info,work=aux)
      if (do_timings) call psb_toc(cg_prec)
      if (do_timings) call psb_tic(cg_vect)
 
      rho_old = rho
      rho     = psb_gedot(r,z,desc_a,info)
  
      if (it == 1) then
        call psb_geaxpby(done,z,dzero,p,desc_a,info)
      else
        if (rho_old == dzero) then
          if (debug_level >= psb_debug_ext_)&
              & write(debug_unit,*) me,' ',trim(name),&
              & ': CG Iteration breakdown rho'
          if (do_timings) call psb_toc(cg_vect)
          exit iteration
        endif
        beta = rho/rho_old
        call psb_geaxpby(done,z,beta,p,desc_a,info)
      end if
      if (do_timings) call psb_toc(cg_vect)
      if (do_timings) call psb_tic(cg_mv)

      call psb_spmm(done,a,p,dzero,q,desc_a,info,work=aux)
      if (do_timings) call psb_toc(cg_mv)
      if (do_timings) call psb_tic(cg_vect)



      sigma = psb_gedot(p,q,desc_a,info)

      if (sigma == dzero) then
        if (debug_level >= psb_debug_ext_) then
          write(debug_unit,*) me,' ',trim(name),': CG Iteration breakdown sigma'
        endif
        exit iteration
      endif
      alpha_old = alpha
      alpha = rho/sigma

      if (do_cond) then 
        istebz = istebz + 1
        if (istebz == 1) then 
          td(istebz) = done/alpha
        else 
          td(istebz) = done/alpha + beta/alpha_old
          tu(istebz-1) = sqrt(beta)/alpha_old
        end if
      end if

  
      call psb_geaxpby(alpha,p,done,x,desc_a,info)
      call psb_geaxpby(-alpha,q,done,r,desc_a,info)



      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_non_,name)
        goto 9999
      end If
 
   end do iteration
  end do restart

  if (do_timings) call psb_toc(cg_vect)
  if (do_cond) then 
   if (me == psb_root_) then 
     cond = dzero
     info=psb_success_
   end if
   call psb_bcast(ctxt,cond)
  end if
  
  
  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr
  
  if (info == psb_success_) call psb_gefree(wwrk,desc_a,info)
  if (info == psb_success_) deallocate(aux,stat=info)
  if (info /= psb_success_) then
   call psb_errpush(info,name)
   goto 9999
  end if
  
  call psb_erractionrestore(err_act)
  return
  
  9999 call psb_error_handler(err_act)
  return

end subroutine psb_dscg_impl


end program
