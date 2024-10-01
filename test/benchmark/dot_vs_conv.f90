program dot_vs_conv
    use psb_base_mod

    ! descriptor
    type(psb_desc_type)                 :: desc_a

    ! dense vectors
    type(psb_d_vect_type)               :: local_vect_1, local_vect_2
    type(psb_s_vect_type)               :: local_vect_1_lower, local_vect_2_lower
  
    ! global vectors memory
    real(psb_dpk_), allocatable         :: temp_local_vect_1(:), temp_local_vect_2(:)
    real(psb_spk_), allocatable         :: temp_local_vect_1_lower(:), temp_local_vect_2_lower(:) 

    ! global vectors memory
    real(psb_dpk_), allocatable         :: global_vect_1(:), global_vect_2(:)
    real(psb_spk_), allocatable         :: global_vect_1_lower(:), global_vect_2_lower(:) 

    ! parallel environment variables
    type(psb_ctxt_type)                 :: ctxt
    integer(psb_ipk_)                   :: my_rank, np, number_of_threads ! Attention, this is a dummy argument for now
  
    ! other variables
    integer(psb_ipk_)                   :: info, i, j, err_act, vector_size, number_of_runs
    character(len=20)                   :: name,ch_err

    ! Result                        
    real(psb_dpk_)                      :: result 
    real(psb_spk_)                      :: result_lower 


    ! Timers
    real(psb_dpk_)                      :: temporary_time 
    real(psb_dpk_), allocatable         :: fill_vectors_time(:), single_scalar_product_timer(:), &
    & double_scalar_product_timer(:), double_to_single(:), single_to_double(:)

    ! Mean values
    real(psb_dpk_)                      :: mean_fill, mean_s_scalar_product, mean_d_scalar_product, mean_s_to_d, mean_d_to_s
    real(psb_dpk_)                      :: var_fill, var_s_scalar_product, var_d_scalar_product, var_s_to_d, var_d_to_s
    real(psb_dpk_)                      :: sum

    ! Number generation
    integer(psb_ipk_)                   :: seed(8)
    real(psb_dpk_)                      :: dummy

    ! Common data
    info              = psb_success_
    name              = 'dot_vs_conv'
    number_of_threads = 1 ! This is a dummy value
    call psb_erractionsave(err_act)

    call psb_init(ctxt)
    call psb_info(ctxt,my_rank,np)

    if ((my_rank < 0).or.(my_rank > np) ) then
        ! This should not happen, but just in case
        info = psb_err_from_subroutine_
        ch_err = 'wrong rank detected'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
    end if

    if( my_rank == psb_root_ ) then
      print '("Insert vector size")'
      read(*,*) vector_size
      print '("Insert number of runs")'
      read(*,*) number_of_runs
    end if

    call psb_bcast(ctxt, vector_size)
    call psb_bcast(ctxt, number_of_runs)
  
    call psb_cdall(ctxt,desc_a,info,nl= ( vector_size / np ))

    if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='allocating descriptor' 
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if


    call psb_geall(local_vect_1,desc_a,info)
    if (info == psb_success_) call psb_geall(local_vect_2,desc_a,info)
    if (info == psb_success_) call psb_geall(local_vect_1_lower,desc_a,info)
    if (info == psb_success_) call psb_geall(local_vect_2_lower,desc_a,info)
      
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocating data structures' 
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call local_vect_1%zero()
    call local_vect_2%zero()
    call local_vect_1_lower%zero()
    call local_vect_2_lower%zero()
    
    call psb_cdasb(desc_a,info)
    if (info == psb_success_) call psb_geasb(local_vect_1,desc_a,info)
    if (info == psb_success_) call psb_geasb(local_vect_2,desc_a,info)
    if (info == psb_success_) call psb_geasb(local_vect_1_lower,desc_a,info)
    if (info == psb_success_) call psb_geasb(local_vect_2_lower,desc_a,info)
      
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='building data structures' 
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(temp_local_vect_1(size(local_vect_1%v%v)))
    allocate(temp_local_vect_2(size(local_vect_2%v%v)))
    allocate(temp_local_vect_1_lower(size(local_vect_1%v%v)))
    allocate(temp_local_vect_2_lower(size(local_vect_2%v%v)))

    allocate(fill_vectors_time(number_of_runs))
    allocate(single_scalar_product_timer(number_of_runs))
    allocate(double_scalar_product_timer(number_of_runs))
    allocate(double_to_single(number_of_runs))
    allocate(single_to_double(number_of_runs))


    ! if(my_rank == 0) write(*,*) size(local_vect_1%v%v)
    ! if(my_rank == 0) write(*,*) size(local_vect_2%v%v)
    ! if(my_rank == 0) write(*,*) size(local_vect_1_lower%v%v)
    ! if(my_rank == 0) write(*,*) size(local_vect_2_lower%v%v)
! 
    ! if(my_rank == 0) write(*,*) size(temp_local_vect_1)
    ! if(my_rank == 0) write(*,*) size(temp_local_vect_2)
    ! if(my_rank == 0) write(*,*) size(temp_local_vect_1_lower)
    ! if(my_rank == 0) write(*,*) size(temp_local_vect_2_lower)
! 
    ! if(my_rank == 0) write(*,*) size(fill_vectors_time)
    ! if(my_rank == 0) write(*,*) size(single_scalar_product_timer)
    ! if(my_rank == 0) write(*,*) size(double_scalar_product_timer)
    ! if(my_rank == 0) write(*,*) size(double_to_single)
    ! if(my_rank == 0) write(*,*) size(single_to_double)

    seed = 21523313
    call random_seed(put=seed)
    
    ! Get data
    do i = 1, number_of_runs
      ! Number generation
      do j = 1, ( (vector_size / np) * my_rank ) * 2
        call random_number(dummy)
      end do
  
      call random_number(temp_local_vect_1)  
      call random_number(temp_local_vect_2)  
  

      ! ---------------------------------------------------

      temporary_time = psb_wtime()

      call local_vect_1%v%bld(temp_local_vect_1)
      call local_vect_2%v%bld(temp_local_vect_2)
      

      fill_vectors_time(i) = psb_wtime() - temporary_time

      
      call psb_amx(ctxt,fill_vectors_time(i))


      ! ---------------------------------------------------

      temporary_time = psb_wtime()

      result = psb_gedot(local_vect_1, local_vect_2, desc_a, info)

      double_scalar_product_timer(i) = psb_wtime() - temporary_time

      call psb_amx(ctxt,double_scalar_product_timer(i))

      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='double scalar product' 
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      ! if(my_rank == 0) then
      !   write(*,'(es20.5)') double_scalar_product_timer
      !   write(*,*) result
      ! end if
      
      ! ---------------------------------------------------

      temporary_time = psb_wtime()

      do j = 1, size(local_vect_1_lower%v%v)
        local_vect_1_lower%v%v(i) = local_vect_1%v%v(i)
      end do
      
      call local_vect_1_lower%v%bld(local_vect_1_lower%v%v)

      double_to_single(i) = psb_wtime() - temporary_time

      call psb_amx(ctxt,double_to_single(i))

      ! if(my_rank == 0) then
      !   write(*,'(es20.5)') double_to_single
      ! end if

      do j = 1, size(local_vect_2_lower%v%v)
        local_vect_2_lower%v%v(i) = local_vect_2%v%v(i)
      end do
      
      call local_vect_2_lower%v%bld(local_vect_2_lower%v%v)


      ! ---------------------------------------------------

      temporary_time = psb_wtime()

      result_lower = psb_gedot(local_vect_1_lower, local_vect_2_lower, desc_a, info)

      single_scalar_product_timer(i) = psb_wtime() - temporary_time

      call psb_amx(ctxt,single_scalar_product_timer(i))



      
      ! ---------------------------------------------------

      temporary_time = psb_wtime()

      do j = 1, size(local_vect_1%v%v)
        local_vect_1%v%v(i) = local_vect_1_lower%v%v(i) * 1.d0
      end do
      
      call local_vect_1%v%bld(local_vect_1%v%v)

      single_to_double(i) = psb_wtime() - temporary_time

      call psb_amx(ctxt,single_to_double(i))

      ! write(40 + my_rank,*) my_rank, i
      ! write(20,*) i, result, result_lower
    end do
    
    call psb_barrier(ctxt)
    

    ! Compute mean
    if(my_rank == psb_root_) then
      sum = 0.d0

      do i = 1, size(fill_vectors_time)
        sum = fill_vectors_time(i) + sum
      end do

      mean_fill = sum / ( size(fill_vectors_time) * 1.d0 )

      sum = 0.d0

      do i = 1, size(double_scalar_product_timer)
        sum = double_scalar_product_timer(i) + sum
      end do

      mean_d_scalar_product = sum / ( size(double_scalar_product_timer) * 1.d0 )

      sum = 0.d0

      do i = 1, size(single_scalar_product_timer)
        sum = single_scalar_product_timer(i) + sum
      end do
      mean_s_scalar_product = sum / ( size(single_scalar_product_timer) * 1.d0 )

      sum = 0.d0

      do i = 1, size(single_to_double)
        sum = single_to_double(i) + sum
      end do

      mean_s_to_d = sum / size(single_to_double)

      sum = 0.d0

      do i = 1, size(double_to_single)
        sum = double_to_single(i) + sum
      end do

      mean_d_to_s = sum / ( size(double_to_single) * 1.d0 )
    end if
    
    call psb_barrier(ctxt)

    if(my_rank == psb_root_) then
      ! Compute variance
      sum = 0.d0

      do i = 1, size(fill_vectors_time)
        sum = sum + ( ( fill_vectors_time(i) - mean_fill ) * ( fill_vectors_time(i) - mean_fill ) )
      end do

      var_fill = sum / size(fill_vectors_time)

      sum = 0.d0

      do i = 1, size(double_scalar_product_timer)
        sum = sum + ( ( double_scalar_product_timer(i) - mean_d_scalar_product ) * & 
        & ( double_scalar_product_timer(i) - mean_d_scalar_product ) ) 
      end do

      var_d_scalar_product = sum / size(double_scalar_product_timer)

      sum = 0.d0

      do i = 1, size(single_scalar_product_timer)
        sum = sum + ( ( single_scalar_product_timer(i) - mean_s_scalar_product ) * & 
        & ( single_scalar_product_timer(i) - mean_s_scalar_product ) )
      end do

      var_s_scalar_product = sum / size(single_scalar_product_timer)

      sum = 0.d0

      do i = 1, size(single_to_double)
        sum = sum + ( ( single_to_double(i) - mean_s_to_d ) * ( single_to_double(i) - mean_s_to_d ) ) 
      end do

      var_s_to_d = sum / size(single_to_double)

      sum = 0.d0

      do i = 1, size(double_to_single)
        sum = sum + (( double_to_single(i) - mean_d_to_s ) * ( double_to_single(i) - mean_d_to_s ) )
      end do

      var_d_to_s = sum / size(double_to_single)


      ! ---------------------------------------------------
      open(unit = 12, file = "results.dat",access = 'append')
      write(12,'("Computation with ",i0, " processes and vector size ", i0, " for ",i0, " runs")') np, vector_size, number_of_runs 
      write(12, '("----------------------------------------------------------------")')
      write(12, '("|  COMPUTATION STEP  |        MEAN        |      VARIANCE      |  ")')
      write(12, '("----------------------------------------------------------------")')
      write(12, '("|  Filling vectors   |", es20.5, "|", es20.5, "|")') mean_fill, var_fill
      write(12, '("|  S scalar prod     |", es20.5, "|", es20.5, "|")') mean_s_scalar_product, var_s_scalar_product
      write(12, '("|  D scalar prod     |", es20.5, "|", es20.5, "|")') mean_d_scalar_product, var_d_scalar_product
      write(12, '("|  D to S conversion |", es20.5, "|", es20.5, "|")') mean_d_to_s, var_d_to_s
      write(12, '("|  S to D conversion |", es20.5, "|", es20.5, "|")') mean_s_to_d, var_s_to_d
      write(12, '("----------------------------------------------------------------")')
      write(12,'(" ")')
      close(12)
    end if

    call psb_barrier(ctxt)

!   
    ! open(unit = 13, file = "double_scalar_product.dat",access = 'append')
    ! write(13,'(es12.5)') double_scalar_product_timer
    ! write(13,'(es12.5)') single_to_double
    ! close(12)
! 
    ! open(unit = 12, file = "single_scalar_product.dat",access = 'append')
    ! write(12,'(es12.5)') single_scalar_product_timer
    ! close(12)
! 
    ! open(unit = 12, file = "double_to_single.dat",access = 'append')
    ! write(12,'(es12.5)') double_to_single
    ! close(12)
! 
    ! open(unit = 12, file = "single_to_double.dat",access = 'append')
    ! write(12,'(es12.5)') double_to_single
    ! close(12)
! 
    call psb_gefree(local_vect_1,desc_a,info)
    call psb_gefree(local_vect_2,desc_a,info)
    call psb_gefree(local_vect_1_lower,desc_a,info)
    call psb_gefree(local_vect_2_lower,desc_a,info)
    
    call psb_cdfree(desc_a,info)

    deallocate(fill_vectors_time)
    deallocate(double_scalar_product_timer)
    deallocate(single_scalar_product_timer)
    deallocate(double_to_single)
    deallocate(single_to_double)

    deallocate(temp_local_vect_1)
    deallocate(temp_local_vect_2)
    deallocate(temp_local_vect_1_lower)
    deallocate(temp_local_vect_2_lower)


    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='free routine'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_exit(ctxt)
    stop

    9999 call psb_error_handler(ctxt,err_act)
    stop

  end program dot_vs_conv