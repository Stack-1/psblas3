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
! Check https://masuday.github.io/fortran_tutorial/random.html to see random number generation



program psb_dscg 
  use psb_base_mod
  use psb_util_mod
  use psb_prec_mod

  implicit none
  
  type(psb_sspmat_type) :: matrix_A_lower_precision

  type(psb_ctxt_type) :: comm_ctxt

  type(psb_d_vect_type)    ::  r_col, b_col




  character(len=40) :: krilov_method,  partition, matrix_file, rhs_file
  character(len=256) :: specs_file_path
  character(len=5) :: matrix_format
  integer(psb_ipk_) :: iret, istopc, max_iteretions,itrace_,irst, err_act
  real(psb_dpk_) :: eps, condition_number, err
  integer(psb_ipk_) :: iteration_number
  character(len=20) :: name
  real(psb_dpk_) :: t1, t2 ! timers variables
  integer(psb_ipk_), parameter :: iunit=12
  real(psb_dpk_), allocatable, target ::  aux_b(:,:)
  integer(psb_ipk_)                   ::  seed_size
  real(psb_dpk_), allocatable, target ::  aux_b_lower_precision(:,:)


  integer :: i

  ! Communication parameters
  integer(psb_ipk_)               :: my_rank, np

  ! Global matrix parameters
  integer(psb_ipk_)               :: global_rows, global_cols, global_nnz ! 4 byte integer
  
  ! Matrix metadata
  type(psb_dspmat_type)           :: local_A
  type(psb_desc_type)             :: desc_a
  type(psb_ldspmat_type)          :: global_A

  ! Data distribution parameters
  integer(psb_ipk_)               :: rows_per_process, process_extra_row
  integer(psb_lpk_), allocatable  :: my_index(:)
  integer(psb_ipk_)               :: process_grid_rows, process_grid_cols

  ! Vectors
  type(psb_d_vect_type)           :: local_b, local_x
  real(psb_dpk_), allocatable     :: global_x(:), global_b(:), global_r(:)

  ! Check variables
  integer(psb_ipk_)               :: info, total_rows_check
  character(len=20)               :: output_file_string

  ! Fill the matrix variables
  integer(psb_ipk_)               :: local_nnz
  integer(psb_lpk_), allocatable  :: row_indexes(:), col_indexes(:)
  real(psb_dpk_), allocatable     :: random_number_vector(:)


  ! Preconditioner metadata
  type(psb_dprec_type)            :: preconditioner
  character(len=40)               :: preconditioner_string

  name = 'mixed_psb_dscg'
  info = psb_success_
  call psb_erractionsave(err_act)


  preconditioner_string = "NONE"
  eps=1.d-6
  max_iteretions = 500
  itrace_ = 0
  irst  = 1
  istopc=1 


  call psb_init(comm_ctxt)
  call psb_info(comm_ctxt, my_rank, np)


  if( (my_rank < 0).or.(my_rank >= np) ) then
    write(*,*) '[ERROR] Unexpected wrong rank value'
    call psb_exit(comm_ctxt)
    stop
  end if

  if( my_rank == psb_root_) then
    print '("=================================================================================")'
    print '("[INFO] Process ",I0 ,": PSBLAS environment started correctly")', my_rank
  end if


  
  if(my_rank == psb_root_) then
    global_rows = 16
    global_cols = 16
    global_nnz  = 30
    print '("[INFO] Process ",i0, ": Matrix size is ",i0," x ",i0,", with ",i0,"/",i0," non zeros values")',&
    & my_rank, global_rows, global_cols, global_nnz, global_rows*global_cols
    print '("[INFO] Process ",i0, ": Data distribution process starting...")', my_rank
  endif
  
  
  call psb_bcast(comm_ctxt,global_rows)
  call psb_bcast(comm_ctxt,global_cols)
  call psb_bcast(comm_ctxt,global_nnz)

  !
  ! Using a simple BLOCK ROW distribution.
  !
  
  rows_per_process = global_rows / np
  if(my_rank < mod(global_rows , np)) then
    rows_per_process = rows_per_process + 1
  end if

  call psb_cdall(comm_ctxt,desc_a,info,nl=rows_per_process)
  call check_correct(info,name,err_act)
  my_index = desc_a%get_global_indices()


  call psb_barrier(comm_ctxt)

  if(my_rank == psb_root_) then
    print '("[INFO] Process ", i0, ": Data distribution indexes correctly set for block row distribution")', my_rank
    print '("=================================================================================")' 
    print '("[INFO] Process ", i0, ": Allocating double precision matrix and vectors ... ")', my_rank
    t1 = psb_wtime() 
  endif

  ! Allocate sparse global matrix
  if (info == psb_success_) then
    call psb_spall(local_A,desc_a,info,nnz=global_nnz,bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
  endif
  ! Define rhs from boundary conditions; also build initial guess
  if (info == psb_success_) then
    call psb_geall(local_x,desc_a,info)
  endif

  if (info == psb_success_) then
    call psb_geall(local_b,desc_a,info,bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
  endif

  
  call psb_barrier(comm_ctxt)

  if(my_rank == psb_root_) then 
    t2 = psb_wtime() - t1
  endif 
  call check_correct(info,name,err_act)

  if(my_rank == psb_root_) then 
    print '("[INFO] Process ", I0, ": Allocation of double precision matrix and vectors ended correctly")', my_rank
    print '("[INFO] Process ", I0, ": Time to allocate double precision matrix and vectors is " , es12.5 )', my_rank, t2
    print '("=================================================================================")'
  endif

  allocate(random_number_vector(global_nnz))
  
  if(my_rank == psb_root_) then
    print '("[INFO] Process ", i0, ": Generating random numbers in random positions of the array...")', my_rank
    ! Generate random numbers
    call generate_random_numbers(random_number_vector, size(random_number_vector), my_rank)
  endif

  call psb_bcast(comm_ctxt,random_number_vector)

  allocate(row_indexes(size(random_number_vector)))
  allocate(col_indexes(size(random_number_vector)))

  if(my_rank == psb_root_) then
    call generate_random_indexes(row_indexes, size(row_indexes), col_indexes, size(col_indexes), my_rank, global_rows, global_cols)
  end if


  call psb_bcast(comm_ctxt,row_indexes)
  call psb_bcast(comm_ctxt,col_indexes)
  call psb_barrier(comm_ctxt)

  call psb_spins(size(random_number_vector), row_indexes, col_indexes, random_number_vector, local_A, desc_a, info )

  if(my_rank == psb_root_) then
    t2 = psb_wtime() - t1
  end if

  ! Write results on file 
  write (output_file_string, '("global_output_",i0,".txt")') my_rank
  
  call psb_gather(global_A,local_A,desc_a,info,root=0)
  call psb_barrier(comm_ctxt)
  
  if(my_rank == psb_root_) then  
    call global_A%print(fname=output_file_string)
  end if
  
  call psb_barrier(comm_ctxt)

  write (output_file_string, '("output_",i0,".txt")') my_rank
  call local_A%print(fname=output_file_string,iv=my_index)
  call psb_barrier(comm_ctxt)

  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Random number generation ended correctly")', my_rank
    print '("[INFO] Process ", I0, ": Time to generate random numbers and fill the matrix is " , es12.5 )', my_rank, t2
    print '("=================================================================================")'
  end if

  ! Now that the matrix is distributed and the descriptor has been initialized
  ! we need to build our preconditioner
  if(preconditioner_string /= "NONE") then
    call preconditioner%init(comm_ctxt,preconditioner_string,info)

    ! building the preconditioner
    t1 = psb_wtime()
    call preconditioner%build(local_A,desc_a,info)
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
    print '("=================================================================================")'
    print '("[INFO] Process ", I0, ": Starting computation...")', my_rank
    global_b = local_b%get_vect()
    write(*,*) global_b
  end if
  


  ! call psb_dscg_impl(local_A,preconditioner,local_b,local_x,eps,desc_a,info, max_iteretions, &
  !                  & iteration_number, err, itrace=itrace_, istop=istopc, cond=condition_number)
  ! call check_correct(info,name,err_act)
  ! call psb_barrier(comm_ctxt)

  ! call psb_gather(global_x,local_x,desc_a,info,root=psb_root_)
  ! call check_correct(info,name,err_act)

  call psb_barrier(comm_ctxt)

  
  if(my_rank == psb_root_) then
    print '("[INFO] Process ", I0, ": Computation ended correctly")', my_rank
    print '("[INFO] Process ", I0, ": Size of solution " , I0)', my_rank, local_x%get_nrows()
    print '("=================================================================================")'
  end if


  ! if(my_rank == psb_root_) then
  !  write(20,*) 'matrix: ',matrix_file
  !  write(20,*) 'computed solution on ',np,' processors.'
  !  write(20,*) 'iterations to convergence: ',iteration_number
  !  write(20,*) 'error estimate (infinity norm) on exit:', &
  !       & ' ||r||/(||a||||x||+||b||) = ',err
  !  write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
  !  do i=1,desc_a%get_global_rows()
  !    write(20,*) i,global_x(i),global_b(i)
  !  enddo
  ! endif

  if(my_rank == 0) then 
    print '("[INFO] Process ", I0, ": Computation ended correctly, check file to see the results")', my_rank
    print '("[INFO] Process ", I0, ": PSBLAS environment shutting down")', my_rank
  endif


  ! Free dynamically allocated memory
  deallocate(random_number_vector)


  call psb_exit(comm_ctxt)

contains

! Debug sobroutine to check if info has risen error
!
! @param info                 -> Integer to return, 0 if no error 
!                                occured, /= 0 otherwise
! @param name                 -> Name of the current program  
! @param err_act              -> Error code for the error handler
!                                - Check PSBLAS user guide to 
!                                  know more
!
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

! This subroutine is usefull to generate random number
! using the default fortran pseudo-random generator
! (NOTE: Should be converted in some other PRG)
! 
! The subrutine returns an array filling it with
! random real numbers (double precision) between 
! 0 - 100
!
! @param random_number_vector -> Real double precision vector to fill
!
! @param vector_size          -> The size of the vector to fill, 
!                                mandatory to dynamically allocate memory
!
! @param nonce                -> Integer used to randomize the seed 
!                                - Usefull when running the same code with
!                                  different processes.
!                                - HINT: Use process PID
!
subroutine generate_random_numbers(random_number_vector,vetor_size,nonce)
  implicit none 

  integer(psb_ipk_)               :: nonce, vetor_size
  integer(psb_ipk_)               :: seed(8)
  real(psb_dpk_)                  :: random_number_vector(vetor_size)

  ! TODO: Generate a random seed and read it from file
  seed = 123456789 + nonce*842342834    ! putting arbitrary seed to all elements
  ! Set seed
  call random_seed(put=seed)
  ! Generate numbers
  call random_number(random_number_vector)  
  ! Normalize output vector
  random_number_vector = random_number_vector*100

end subroutine


! This subroutine is usefull to generate random indexes for the matrix
! to fill. The subroutine uses the default fortran pseudo-random 
! generator.
! (NOTE: Should be converted in some other PRG)
! 
! The subrutine returns an array filling it with random real numbers 
! (double precision) between 1 - max
! 
!
! @param row_indexes          -> Integer vector to fill with row indexes
!
! @param rows_size            -> The size of the rows vector to fill, 
!                                mandatory to dynamically allocate memory
!
! @param col_indexes          -> Integer vector to fill with col indexes
!
! @param col_size             -> The size of the cols vector to fill, 
!                                mandatory to dynamically allocate memory
!
! @param nonce                -> Integer used to randomize the seed 
!                                - Usefull when running the same code with
!                                  different processes.
!                                - HINT: Use process PID
!
! @param max_rows             -> The maximum number to generate for rows
!
! @param max_cols             -> The maximum number to generate for cols
subroutine generate_random_indexes(row_indexes,rows_size,col_indexes, cols_size, nonce, max_rows, max_cols)
  implicit none 

  integer(psb_ipk_)               :: rows_size, cols_size, max_rows, max_cols, nonce, i, j
  integer(psb_lpk_)               :: row_indexes(rows_size), col_indexes(cols_size)
  integer(psb_ipk_)               :: seed(8)
  real(psb_dpk_)                  :: gen_row, gen_col
  logical                         :: exists

  ! TODO: Generate a random seed and read it from file
  seed = 37215645 + nonce*23723644    ! putting arbitrary seed to all elements
  ! Set seed
  call random_seed(put=seed)

  ! Generate numbers
  do i=1, rows_size
    1111 call random_number(gen_row)
    call random_number(gen_col)

    ! Normalize obtained numbers
    gen_row = ( ( gen_row * 10 * (max_rows - 1) ) / 10 ) + 1
    gen_col = ( ( gen_col * 10 * (max_cols - 1) ) / 10 ) + 1

    exists = .false.
    do j = 0, i
      ! Check if there already exist this couple of indexes
      if( (row_indexes(j) == int(gen_row)).and.(col_indexes(j) == int(gen_col) ) ) then
        exists = .true.
        exit
      end if 
    end do

    if(.not.(exists)) then
      row_indexes(i) = int(gen_row,8)
      col_indexes(i) = int(gen_col,8)
    else
      goto 1111
    end if


  end do

  
  
end subroutine


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
    if(my_rank == psb_root_) then 
      write(*,*) x%get_nrows(), r%get_nrows()
    end if
    if (info == psb_success_) call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
    if(my_rank == psb_root_) then 
      print '("[DEBUG] Error in sparse matrix multiplication")'
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
