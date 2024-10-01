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

program psb_dscg
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use psb_ds_gen_mod
  use psb_d_cg
  use psb_s_cg
  use psb_ds_cg_1
  use psb_ds_cg_2
#ifdef HAVE_CUDA
  use psb_ext_mod
  use psb_cuda_mod
#endif

  implicit none

  ! input parameters
  character(len=20)                     :: krylov_method, prec_type
  character(len=5)                      :: afmt
  integer(psb_ipk_)                     :: idim
  integer(psb_epk_)                     :: system_size

  ! timers
  real(psb_dpk_)                        :: initial_time, temporary_time, generation_time, preconditioning_time, &
                                          & computation_time, residual_computation_time, gpu_convertion_time, total_time

  ! sparse matrix
  type(psb_dspmat_type)                 :: local_a, local_a_gpu
  type(psb_sspmat_type)                 :: local_a_lower_precision, local_a_lower_gpu

  ! preconditioner
  type(psb_dprec_type)                  :: prec
  type(psb_sprec_type)                  :: prec_lower


  ! descriptor
  type(psb_desc_type)                   :: desc_a


  ! dense vectors
  type(psb_d_vect_type)                 :: local_x,local_b, local_r
  type(psb_s_vect_type)                 :: local_x_lower_precision, local_b_lower_precision, local_r_lower_precision

  real(psb_dpk_), allocatable, save     :: global_x(:), global_r(:), global_b(:)

  real(psb_spk_), allocatable, save     :: global_x_lower_precision(:), global_r_lower_precision(:), global_b_lower_precision(: )

  ! parallel environment variables
  type(psb_ctxt_type)                   :: ctxt
  integer(psb_ipk_)                     :: my_rank, np, number_of_threads_per_process ! Attention, this is a dummy argument for now

  ! solver parameters
  integer(psb_ipk_)                     :: iter, itmax,itrace, istopc, irst
  integer(psb_epk_)                     :: matrix_memory_size, prec_memory_size, desc_memory_size ! These are all variables containing the memory occupation in bytes
  real(psb_dpk_)                        :: err, eps
  real(psb_spk_)                        :: err_lower, eps_lower


  ! debug output variables
  character(len=40)                     :: output_file_string

  ! cuda variables
#ifdef HAVE_CUDA
  type(psb_d_vect_cuda)                 :: gpu_vector_format_double
  type(psb_d_cuda_elg_sparse_mat)       :: gpu_matrix_format_double

  type(psb_s_vect_cuda)                 :: gpu_vector_format_single
  type(psb_s_cuda_elg_sparse_mat)       :: gpu_matrix_format_single

  type(psb_i_vect_cuda)                 :: gpu_descriptor_format
#endif

  ! Parameters for solvers in Block-Jacobi preconditioner
  type ainvparms
    character(len=12) :: alg, orth_alg, ilu_alg, ilut_scale
    integer(psb_ipk_) :: fill, inv_fill
    real(psb_dpk_)    :: thresh, inv_thresh
  end type ainvparms
  type(ainvparms)     :: parms

  ! other variables
  integer(psb_ipk_) :: info, i, err_act, precision_mode
  character(len=20) :: name,ch_err
  character(len=40) :: fname
  real(psb_dpk_)    :: r_amax, b_amax, scale,resmx,resmxp, condition_number

  ! Common data
  info = psb_success_
  name = 'psb_dscg'
  number_of_threads_per_process = 1 ! This is a dummy value, every process is launched using only a thread
  call psb_erractionsave(err_act)
  
  ! Standard input for krylov PSBLAS standard benchmark
  krylov_method = "CG"
  prec_type = "NONE"
  afmt = "CSR"
  istopc = 2 ! 1 - Normwise backword error, 2 - Relative residual in the 2-norm , 3 - relative residual reduction in the 2-norm
  itmax = 10000
  itrace = 0 ! Debug option on
  irst = 0 ! Restart parameter, ignored for CG

  call psb_init(ctxt)
  call psb_info(ctxt,my_rank,np)
#ifdef HAVE_CUDA
  call psb_cuda_init(ctxt)
#endif


  if ((my_rank < 0).or.(my_rank > np) ) then
    ! This should not happen, but just in case
    info = psb_err_from_subroutine_
    ch_err = 'wrong rank detected'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  if(my_rank == psb_root_) then 
    print '("=================================================================================")'
    print '("=                                 WELCOME                                       =")'
    print '("=================================================================================")'
    print '("[INFO] Welcome to PSBLAS version: ", a )', psb_version_string_
    print '("[INFO] This is the ", a ," sample program")', trim(name)
#ifdef HAVE_CUDA
    print '("[INFO] This program will be computed on GPU ", a )' , psb_cuda_DeviceName() 
#endif
    print '(" ")'
    print '("=================================================================================")'
    print '("=                              INPUT SECTION                                    =")'
    print '("=================================================================================")'
    print '("[INPUT] Insert squared dimension of the matrix")'
    read(*,*) idim
    print '("[INPUT] Insert flating point mode")'
    print '("        1 - Double precision")'
    print '("        2 - Single precision")'
    print '("        3 - Mixed precision 1")'
    print '("        4 - Mixed precision 2")'
    print '("        5 - PSBLAS standatd CG with double precision")'
    print '("        6 - PSBLAS standatd CG with single precision")'
    read(*,*) precision_mode

    if((precision_mode < 1).or.(precision_mode > 7)) then
      print '("[ERROR] Invalid precision selected")'
      goto 9999
    end if
    print '(" ")'
  end if 

  ! Synchronize metadata
  call psb_bcast(ctxt, precision_mode)
  call psb_bcast(ctxt, idim)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess
  !
  initial_time = psb_wtime()
  call psb_ds_gen_matrix(ctxt,idim,local_a,local_b,local_x,&
      & local_a_lower_precision,local_b_lower_precision,local_x_lower_precision, &
      & desc_a,afmt,info)

  generation_time = psb_wtime() - initial_time

  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_ds_gen_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  
 


  if (my_rank == psb_root_) then
    print '("[INFO] Overall matrix creation time : ",es12.5, " s")' , generation_time
    print '(" ")'   
  end if

  !
  !  prepare the preconditioner.
  !
  if(my_rank == psb_root_) then
    print '("=================================================================================")'
    print '("=                         PRECONDITIONER SECTION                                =")'
    print '("=================================================================================")'
    write(psb_out_unit,'("[INFO] Setting preconditioner to:         ",a)') prec_type
  end if
  
  temporary_time = psb_wtime()

  call prec%init(ctxt,prec_type,info)
  call prec_lower%init(ctxt,prec_type,info)

  !
  ! Set the options for the BJAC preconditioner
  !
  if (psb_toupper(prec_type) == "BJAC") then
      call prec%set('sub_solve',       parms%alg,   info)
      select case (psb_toupper(parms%alg))
      case ("ILU")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('ilu_alg',         parms%ilu_alg,    info)
      case ("ILUT")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('sub_iluthrs',     parms%thresh,     info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
      case ("AINV")
        call prec%set('inv_thresh',      parms%inv_thresh, info)
        call prec%set('inv_fillin',      parms%inv_fill,   info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
        call prec%set('ainv_alg',        parms%orth_alg,   info)
      case ("INVK")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('inv_fillin',      parms%inv_fill,   info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
      case ("INVT")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('inv_fillin',      parms%inv_fill,   info)
        call prec%set('sub_iluthrs',     parms%thresh,     info)
        call prec%set('inv_thresh',      parms%inv_thresh, info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
      case default
        ! Do nothing, use default setting in the init routine
      end select
  else
    ! nothing to set for NONE or DIAG preconditioner
  end if



  call psb_barrier(ctxt)

  call prec%build(local_a,desc_a,info)
  call prec_lower%build(local_a_lower_precision,desc_a,info)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_precbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  preconditioning_time = psb_wtime() - temporary_time

  call psb_amx(ctxt,preconditioning_time)

  if (my_rank == psb_root_) then
    print '("[INFO] Overall preconditioning time:   ",es12.5)' , preconditioning_time
    print '(" ")'
  end if 
  call prec%descr(info)
  call prec_lower%descr(info)
  !
  ! iterative method parameters
  !
  call psb_barrier(ctxt)
  eps         = 1.d-6
  eps_lower   = eps

  if(my_rank == psb_root_) then
    print '("=================================================================================")'
    print '("=                        ITERATIVE SOLVER SECTION                               =")'
    print '("=================================================================================")'
    write(psb_out_unit,'("[INFO] Starting computation... ") ')
  end if


#ifdef HAVE_CUDA
  temporary_time = psb_wtime()

  ! marking data structures to use them in GPU
  call local_a%cscnv(local_a_gpu,info,mold=gpu_matrix_format_double)
  if(info == psb_success_) call local_a_lower_precision%cscnv(local_a_lower_gpu,info,mold=gpu_matrix_format_single) 

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='gpu convert mat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  
  call desc_a%cnv(mold=gpu_descriptor_format)

  call psb_geasb(local_b,desc_a,info,scratch=.false.,mold=gpu_vector_format_double)
  if(info == psb_success_) call psb_geasb(local_x,desc_a,info,scratch=.false.,mold=gpu_vector_format_double)
  if(info == psb_success_) call psb_geasb(local_b_lower_precision,desc_a,info,scratch=.false.,mold=gpu_vector_format_single)
  if(info == psb_success_) call psb_geasb(local_x_lower_precision,desc_a,info,scratch=.false.,mold=gpu_vector_format_single)


  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='gpu convert vectors'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call local_a_gpu%cscnv(info,mold=gpu_matrix_format_double)
  call local_a_lower_gpu%cscnv(info,mold=gpu_matrix_format_single)

  gpu_convertion_time = psb_wtime() - temporary_time

  call psb_geall(local_r,desc_a,info)
  call psb_geall(local_r_lower_precision,desc_a,info)

  call psb_geasb(local_r_lower_precision,desc_a,info,scratch=.false.,mold=gpu_vector_format_single)

  write(*,*) local_b_lower_precision%v%v
  print '(" ")'

  call local_b_lower_precision%set_host()
  call local_r_lower_precision%set_host()
  

  call psb_geaxpby(sone, local_b_lower_precision, szero, local_r_lower_precision, desc_a, info)
  call psb_cuda_DeviceSync()

  write(*,*) local_r_lower_precision%get_vect()

  print '(" ")'

  call psb_geasb(local_r,desc_a,info,scratch=.false.,mold=gpu_vector_format_double)

  write(*,*) local_b%v%v
  print '(" ")'

  call local_b%set_host()
  call local_r%set_host()

  call psb_geaxpby(done, local_b, dzero, local_r, desc_a, info)

  call psb_cuda_DeviceSync()
  write(*,*) local_r%get_vect()

  print '(" ")'

#endif


  temporary_time = psb_wtime()

  if(precision_mode == 1) then
    ! Double precision computation
    if(my_rank == psb_root_) print '("[INFO] Calling double precision Conjugate Gradient")'

    call psb_dcg_impl(local_a_gpu,prec,local_b,local_x,eps,desc_a,info,&
    & itmax=itmax,iter=iter,err=err)
  else if(precision_mode == 2) then
    ! Single precision computation
    if(my_rank == psb_root_) print '("[INFO] Calling single precision Conjugate Gradient")'

    call psb_scg_impl(local_a_lower_gpu,prec,local_b_lower_precision,local_x_lower_precision, &
    & eps,desc_a,info, itmax=itmax,iter=iter,err=err)
  else if(precision_mode == 3) then  
    ! Mixed precision computation
    if(my_rank == psb_root_) print '("[INFO] Calling mixed precision 1 Conjugate Gradient")'
    
    call psb_dscg_1_impl(local_a_lower_gpu,prec,local_b,local_x, &
    & eps,desc_a,info, itmax=itmax,iter=iter,err=err)
  else if(precision_mode == 4) then  
    ! Mixed precision computation
    if(my_rank == psb_root_) print '("[INFO] Calling mixed precision 2 Conjugate Gradient")'
    
    call psb_dscg_2_impl(local_a_lower_gpu,prec,local_b_lower_precision,local_x_lower_precision, &
    & eps,desc_a,info, itmax=itmax,iter=iter,err=err)
  else if(precision_mode == 5) then  
    ! PSBLAS CG double computation computation
    if(my_rank == psb_root_) print '("[INFO] Calling iterative method ",a)', krylov_method

    call psb_krylov(krylov_method,local_a_gpu,prec,local_b,local_x,eps,desc_a,info,&
    & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)

  else if(precision_mode == 6) then  
    ! PSBLAS CG double computation computation
    if(my_rank == psb_root_) print '("[INFO] Calling iterative method ",a)', krylov_method

    call psb_krylov(krylov_method,local_a_lower_gpu,prec_lower,local_b_lower_precision,local_x_lower_precision,&
    & eps_lower,desc_a,info,itmax=itmax,iter=iter,err=err_lower,itrace=itrace,istop=istopc,irst=irst)
  else
    print '("[ERROR] Invalid call to psb_dscg")'
    call psb_exit(ctxt)
    stop
  end if 

#ifdef HAVE_CUDA
  call psb_cuda_DeviceSync()
#endif

  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  computation_time = psb_wtime() - temporary_time
  
  call psb_amx(ctxt,computation_time)

  
  temporary_time = psb_wtime()

  ! Compute the explicit residual
  if((precision_mode == 1).or.(precision_mode == 3).or.(precision_mode == 5)) then
    call psb_geall(local_r, desc_a, info)
    ! r = b
    call psb_geaxpby(done, local_b, dzero, local_r, desc_a, info)

    call psb_geasb(local_r, desc_a, info)
      ! r = r - A * x
    call psb_spmm(-done, local_a, local_x, done, local_r, desc_a, info)
    matrix_memory_size = local_a%sizeof()
  else if((precision_mode == 2).or.(precision_mode == 4).or.(precision_mode == 6)) then
    call psb_geall(local_r_lower_precision, desc_a, info)
    ! r = b
    call psb_geaxpby(sone, local_b_lower_precision, szero, local_r_lower_precision, desc_a, info)
    
    call psb_geasb(local_r_lower_precision, desc_a, info)
      ! r = r - A * x
    call psb_spmm(-sone, local_a_lower_precision, local_x_lower_precision, sone, local_r_lower_precision, desc_a, info)

    matrix_memory_size = local_a_lower_precision%sizeof()
  end if

  residual_computation_time = psb_wtime() - temporary_time

  desc_memory_size = desc_a%sizeof()
  prec_memory_size = prec%sizeof()
  system_size = desc_a%get_global_rows()

  call psb_sum(ctxt,matrix_memory_size)
  call psb_sum(ctxt,desc_memory_size)
  call psb_sum(ctxt,prec_memory_size)

  if (my_rank == psb_root_) then
    write(psb_out_unit,'(" ")')
    print '("=================================================================================")'
    print '("=                         TERMINAL OUTPUT SECTION                               =")'
    print '("=================================================================================")'
    print '("[INFO] Printing output statistics")'   
    print '("       - Number of processes                     : ",i0)'            , np
    print '("       - Number of threads                       : ",i0)'            , number_of_threads_per_process
    print '("       - Total number of tasks                   : ",i0)'            , number_of_threads_per_process*np
    print '("       - Linear system size                      : ",i0," x ",i0)'   , system_size, system_size
    print '("       - Time to compute residual                : ",es0.5)'         , residual_computation_time
#ifdef HAVE_CUDA
    print '("       - Time to initialize GPU data structures  : ",es0.5)'         , gpu_convertion_time
#endif
    print '("       - Time to solve system                    : ",es0.5)'         , computation_time
    print '("       - Time per iteration                      : ",es0.5)'         , computation_time/iter
    print '("       - Number of iterations                    : ",i0)'            , iter
    if(precision_mode == 6) then
      print '("       - Convergence indicator on exit           : ",es0.5)'       , err_lower
    else
      print '("       - Convergence indicator on exit           : ",es0.5)'       , err
    end if
    print '("       - Info on exit                            : ",i0)'            , info
    print '("       - Total memory occupation for A           : ",i0,"  bytes")'  , matrix_memory_size
    print '("       - Total memory occupation for PREC        : ",i0,"  bytes")'  , prec_memory_size 
    print '("       - Total memory occupation for DESC_A      : ",i0,"  bytes")'  , desc_memory_size
    print '("       - Storage format for A                    : ",a)'             , local_a%get_fmt()
    print '("       - Storage format for DESC_A               : ",a)'             , desc_a%get_fmt()
  end if


  call psb_gather(global_x,local_x,desc_a,info,root=psb_root_)
  if(info == psb_success_) call psb_gather(global_x_lower_precision,local_x_lower_precision,desc_a,info,root=psb_root_)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='gathering x')
    goto 9999
   end if

  call psb_barrier(ctxt)



  call psb_gather(global_b,local_b,desc_a,info,root=psb_root_)
  if(info == psb_success_) call psb_gather(global_b_lower_precision,local_b_lower_precision,desc_a,info,root=psb_root_)
  
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='gathering b')
    goto 9999
   end if

  call psb_barrier(ctxt)

  
  if( (precision_mode == 1).or.(precision_mode == 3).or.(precision_mode == 5) ) then 
    call psb_gather(global_r,local_r,desc_a,info,root=psb_root_)
  else if( (precision_mode == 2).or.(precision_mode == 4).or.(precision_mode == 6) ) then
    call psb_gather(global_r_lower_precision,local_r_lower_precision,desc_a,info,root=psb_root_)
  end if

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='gathering r')
    goto 9999
   end if

  call psb_barrier(ctxt)


  
  if(my_rank == psb_root_) then
    write(output_file_string, '("../data/gpu/final_result_",i0,".txt")') precision_mode
    open (unit = 20, file = output_file_string)
    write(20,*) 'computed solution on ',np,' processors.'
    write(20,*) 'iterations to convergence: ',iter
    write(20,*) 'time to convergence: ',computation_time
    write(20,*) 'matrix size: ',idim * idim, ' x ', idim * idim
    if(precision_mode == 6) then
      write(20,*) 'error estimate on exit:', &
      & ' ||r|| / ||b|| = ',err_lower
    else
      write(20,*) 'error estimate on exit:', &
      & ' ||r|| / ||b|| = ',err
    end if
    
    write(20,'(a8,4(2x,a20))') 'I','X(I)','B(I)','R(I)'


    if( (precision_mode == 1).or.(precision_mode == 3).or.(precision_mode == 5) ) then 
      do i=1,100
        write(20,*) i ,global_x(i), global_b(i) , global_r(i)
      end do
    else if( (precision_mode == 2).or.(precision_mode == 4).or.(precision_mode == 6) ) then
      do i=1,100
        write(20,*) i, global_x_lower_precision(i), global_b_lower_precision(i), global_r_lower_precision(i)
      end do
    else
      close(20)
      goto 9999
    end if
    close(20)
  endif


  call psb_barrier(ctxt)

  !
  !  cleanup storage and exit
  !
  call psb_gefree(local_b,desc_a,info)
  call psb_gefree(local_x,desc_a,info)
  call psb_spfree(local_a,desc_a,info)
  if(precision_mode == 4) call psb_gefree(local_r,desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

#ifdef HAVE_CUDA
  call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop

end program psb_dscg

