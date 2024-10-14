module psb_s_cg

    contains
        subroutine psb_scg_impl(a,prec,b,x,error_stopping_criterion,desc_a,info,&
            & itmax,iter,err)
            use psb_base_mod
            use psb_prec_mod
            use psb_krylov_mod
#ifdef HAVE_CUDA
            use psb_ext_mod
            use psb_cuda_mod
#endif
            implicit none
        
            class(psb_dprec_type), intent(inout)            :: prec
        
            real(psb_dpk_), Intent(in)                      :: error_stopping_criterion
            integer(psb_ipk_), intent(out)                  :: info
            integer(psb_ipk_), Optional, Intent(in)         :: itmax
            integer(psb_ipk_), Optional, Intent(out)        :: iter
            real(psb_spk_), Optional, Intent(out)           :: err
            ! =   Local data
            integer(psb_ipk_)                               :: itmax_, istop_, it, itx,&
                                                            &  n_col, n_row,err_act, ieg,nspl, istebz
            integer(psb_dpk_)                               :: mglob
            integer(psb_dpk_)                               :: debug_level, debug_unit
            type(psb_ctxt_type)                             :: ctxt
            integer(psb_dpk_)                               :: np, my_rank
            character(len=20)                               :: name, ch_err
            
            ! New variables
            type(psb_sspmat_type), Intent(in)               :: a
            Type(psb_desc_type)                             :: desc_a
            type(psb_s_vect_type), Intent(inout)            :: b, x
            type(psb_d_vect_type)                           :: b_double, r_double

            type(psb_s_vect_type)                           :: d, rho
            type(psb_s_vect_type)                           :: r
        
            real(psb_spk_)                                  :: beta, alpha
            real(psb_spk_)                                  :: r_scalar_product, r_scalar_product_next, partial_result_d_rho
        
            ! Check exit variables
            integer(psb_ipk_)                               :: i
        
            ! Norm variables
            real(psb_dpk_)                                  :: r_norm, x_norm, a_norm, b_norm
            
#ifdef HAVE_CUDA
            type(psb_s_vect_cuda)                           :: gpu_vector_format
#endif


            info = psb_success_
            name = 'mixed_psb_dscg'
            call psb_erractionsave(err_act)

            ! Initialize parameters
            itx                   = 0
            alpha                 = szero
            r_scalar_product_next = szero
        
            ctxt = desc_a%get_context()
            call psb_info(ctxt, my_rank, np)

            if (.not.allocated(b%v)) then 
                info = psb_err_from_subroutine_    
                call psb_errpush(info,name,a_err='b not allocated')
                goto 9999
            endif

            if (.not.allocated(x%v)) then 
                info = psb_err_from_subroutine_    
                call psb_errpush(info,name,a_err='x not allocated')
                goto 9999
            endif

            if (.not.allocated(a%a)) then 
                info = psb_err_from_subroutine_    
                call psb_errpush(info,name,a_err='a not allocated')
                goto 9999
            endif

        
            ! Allocate vectors
            call psb_geall(r,desc_a,info)
            if(info == psb_success_) call psb_geall(d,desc_a,info)
            if(info == psb_success_) call psb_geall(rho,desc_a,info)
            if(info == psb_success_) call psb_geall(b_double,desc_a,info)
            if(info == psb_success_) call psb_geall(r_double,desc_a,info)


            if(info /= psb_success_) then
                info = psb_err_from_subroutine_    
                call psb_errpush(info,name,a_err='allocating local vectors')
                goto 9999
            end if





#ifdef HAVE_CUDA
            ! marking data structures to use them in GPU
            call psb_geasb(r,desc_a,info,scratch=.false.,mold=gpu_vector_format)
            call psb_geasb(d,desc_a,info,scratch=.false.,mold=gpu_vector_format)
            call psb_geasb(rho,desc_a,info,scratch=.false.,mold=gpu_vector_format)

            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='gpu convert vectors'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if

#endif


            restart: do 
                ! =   
                ! =    r_0 = b - A * x_0
                ! =   
                if (itx >= itmax) exit restart 
                it = 0

                call psb_geaxpby(sone,b,szero,r,desc_a,info) ! r_0 = b_0


                if(info /= psb_success_) then
                    info=psb_err_from_subroutine_
                    ch_err='r_0 = b_0'
                    call psb_errpush(info,name,a_err=ch_err)
                    goto 9999
                end if

                call psb_spmm(-sone,a,x,sone,r,desc_a,info) ! r_0 = -A * x_0 + r_0

                if(info /= psb_success_) then
                    info=psb_err_from_subroutine_
                    ch_err='computing r_0 in single'
                    call psb_errpush(info,name,a_err=ch_err)
                    goto 9999
                end if
                ! computed r_0 = b - A * x_0

                call psb_geaxpby(sone,r,szero,d,desc_a,info) ! d_0 = r_0

                ! This is the actual CG method 
                iteration:  do   
                    it   = it + 1
                    itx = itx + 1

                    if(it > itmax) exit restart
                    
                    if(it == 1) then
                        r_scalar_product    = psb_gedot(r,r,desc_a,info)   ! r_i * r_i
                    else
                        r_scalar_product    = r_scalar_product_next
                    end if

                    ! call prec%apply(r,z,desc_a,info,work=aux)
            
            
                    call psb_spmm(sone,a,d,szero,rho,desc_a,info) ! rho_i = A * d_i

                    
                    partial_result_d_rho      = psb_gedot(d,rho,desc_a,info) ! d_i * rho_i
                    alpha                     = r_scalar_product / partial_result_d_rho

                    call psb_geaxpby(alpha,d,sone,x,desc_a,info)        ! x_i+1 = x_i + alpha_i * rho_i

                    call psb_geaxpby(-alpha,rho,sone,r,desc_a,info)     ! r_i+1 = r_i - alpha_i * rho_i

                    ! ||r|| / ||b||
                    r_norm = psb_norm2(r, desc_a, info)
                    b_norm = psb_norm2(b, desc_a, info)
                    err = r_norm / b_norm
                    ! write(*,'(i0," ", es15.9," ", es15.9)') it, err, error_stopping_criterion
                    ! write(*,'(i0," ", es15.9," ", es15.9)') it, r_norm, b_norm


                    if(err < error_stopping_criterion) then
                        exit restart
                    end if


                    r_scalar_product_next = psb_gedot(r,r,desc_a,info)  ! r_i+1 * ri+1
                    beta = r_scalar_product_next / r_scalar_product


                    call psb_geaxpby(sone,r,beta,d,desc_a,info)       ! d_i+1 = r_i+1 + beta_i+1 * d_i




            
                end do iteration
            end do restart
            
            
        
            if(my_rank == psb_root_) then
                iter = it
            end if
            
            call psb_bcast(ctxt, iter)
            
        
            
        
            call psb_erractionrestore(err_act)
            return
            
            9999 call psb_error_handler(err_act)
            return
        
        end subroutine psb_scg_impl

end module psb_s_cg