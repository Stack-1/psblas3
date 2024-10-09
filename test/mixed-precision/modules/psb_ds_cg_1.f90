!
!   Parallel Sparse BLAS  version 3.5 
!       (C) Copyright 2006-2018
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
!   Author: Simone Staccone (Stack-1)
!
!   This module implements a mixed bersion of the conjugate gradient algorithm
!   implementing the following pseudo-code:
!   
!   psb_ds_cg(float A,float b)
!       d_0 = r_0 = b - A * x_0
!       i = 0
!       while{max(r_single) < error_stopping_criterion}
!           rho_i       = A * d_i
!           alpha_i     = (r_i * r_i) / (d_i * rho_i)
!           x_i+1       = x_i + (alpha_i * d_i)
!           r_i+1       = r_i - (alpha_i * rho_i)
!           beta_i+i    = (r_i+1 * r_i+1) / (r_i * r_i)
!           d_i+1       = d_i - (beta_i+1 * d_i)
!       end
!       return x_i
!
!   The main idea is to keep orthogonalization properties of d_double vectors
!   that compose the base of the krilov subspace and therefore keep the vector in
!   double precision to mitigate approxiamtion errors.
!
!   In particular the whole code is devveloped as a part of a process parallel environment
!   exploiting the use of PSBLAS library, so this particular algorithm will be excecuted
!   in parallel by multiple processes on a subportion of the whole linear system.
!
!   See paper of this code on: 
!

module psb_ds_cg_1
    
    contains
    subroutine psb_dscg_1_impl(a_single,prec,b_single,x_single,&
        & error_stopping_criterion,desc_a,info,itmax,iter,err)
        use psb_base_mod
        use psb_prec_mod
        use psb_util_mod
        use psb_krylov_mod

        implicit none
        
        ! Preconditioner variables  
        class(psb_dprec_type), intent(inout)            :: prec
        
        real(psb_dpk_), intent(in)                      :: error_stopping_criterion
        integer(psb_ipk_), intent(out)                  :: info
        integer(psb_ipk_), Optional, Intent(in)         :: itmax
        integer(psb_ipk_), Optional, Intent(out)        :: iter
        real(psb_dpk_), Optional, Intent(out)           :: err

        ! Computation metadata variables    
        type(psb_ctxt_type)                             :: ctxt
        integer(psb_dpk_)                               :: np, my_rank
        character(len=20)                               :: name, ch_err
        
        ! Matrix related variables
        type(psb_sspmat_type), intent(in)               :: a_single
        type(psb_desc_type), Intent(in)                 :: desc_a
        
        ! Single precision vectors
        type(psb_s_vect_type), Intent(inout)            :: b_single, x_single
        type(psb_s_vect_type)                           :: rho_single, r_single, d_single


        ! Single precision variables
        real(psb_spk_)                                  :: beta, alpha
        real(psb_spk_)                                  :: r_scalar_product, r_scalar_product_next, partial_result_d_rho
        
        ! Double precision vectors
        type(psb_d_vect_type)                           :: r_double, d_double
        !type(psb_d_vect_type), Intent(inout)            :: b_double, x_double
        !type(psb_dspmat_type), intent(in)               :: a_double

        ! PSBLAS utility variables
        integer(psb_ipk_)                               :: i, j, err_act, it, itx, m, n, ir, jc, nnz

        ! Norm variables
        real(psb_dpk_)                                  :: r_norm, x_norm, a_norm, b_norm

        logical                                         :: stagnation

        real(psb_spk_), allocatable                     :: global_x(:), global_a(:)
    
        character(len=20) :: output_file_name
        integer(psb_lpk_), allocatable :: ltg(:)

        info = psb_success_
        name = 'mixed_psb_dscg'
        call psb_erractionsave(err_act)
    
        ! Initialize parameters
        itx                         = 0
        alpha                       = szero
        r_scalar_product_next       = szero

        ctxt = desc_a%get_context()
        call psb_info(ctxt, my_rank, np)
    
    
        if (.not.allocated(b_single%v)) then 
            info = psb_err_invalid_vect_state_
            call psb_errpush(info,name)
            goto 9999
        endif
        if (.not.allocated(x_single%v)) then 
            info = psb_err_invalid_vect_state_
            call psb_errpush(info,name)
            goto 9999
        endif
        if (.not.allocated(a_single%a)) then 
            info = psb_err_invalid_vect_state_
            call psb_errpush(info,name)
            goto 9999
        endif

        allocate(global_x(size(b_single%v%v)))

        ! Allocate vectors
        call psb_geall(r_single,desc_a,info)
        call psb_geall(r_double,desc_a,info)
        call psb_geall(d_double,desc_a,info)
        call psb_geall(d_single,desc_a,info)
        call psb_geall(rho_single,desc_a,info)


        if(info /= psb_success_) then
            info=psb_err_from_subroutine_    
            call psb_errpush(info,name,a_err='allocating local vectors')
            goto 9999
        end if

        
        it = 0
        itx = 0
        stagnation = .false.

        b_norm = psb_norm2(b_single, desc_a, info)

        n       = a_single%get_ncols()
        m       = a_single%get_nrows()
        nnz     = a_single%get_nzeros()

        !allocate(ltg(n))
!
        !do i = 1, n
        !    ltg(i) = i
        !end do
!
        !write(output_file_name, '("dasdsd_",i0)') my_rank
!
!
        !call a_single%print(fname=output_file_name, iv=ltg)

        ! Restart function should be implemented to helps stabilize the computation
        restart: do 
            ! =   
            ! =    r_0 = b - A * x_0
            ! =   
            if(itx >= itmax) exit restart 
            it = 0

            ! r_0 = b_0
            do i = 1, size(r_double%v%v) 
                r_double%v%v(i) = real(b_single%v%v(i),8)
            end do


            select type(new_a => a_single%a)
                class is(psb_s_coo_sparse_mat)

                    call psb_gather(global_x, x_single, desc_a, info)


                    ! r_0 = - A * x_0 + r_0
                    do i=1,nnz
                        ir = new_a%ia(i)
                        jc = new_a%ja(i)
                        r_double%v%v(ir) = r_double%v%v(ir) &
                        & - ( real(new_a%val(i),8) * real(global_x(jc),8) )
                    end do

                class default
                    write(*,'("Error in class initialization")') 
            end select

            ! computed r_0 = b - A * x_0
            r_norm = psb_norm2(r_double, desc_a, info)

            err = r_norm / b_norm
            write(13,*) itx, err, error_stopping_criterion, r_norm
            if(err < error_stopping_criterion) then
                exit restart
            end if

            ! Initialize d
            ! d_0 = r_0 

            do i=1,size(d_single%v%v)
                r_single%v%v(i)     = r_double%v%v(i)
                d_double%v%v(i)     = r_double%v%v(i)
                d_single%v%v(i)     = d_double%v%v(i)
            end do
            
            ! This is the actual CG method 
            iteration:  do   
                it   = it + 1
                itx = itx + 1
                if(it > itmax) exit restart
                
                if(it == 1) then
                    r_scalar_product    = psb_gedot(r_single,r_single,desc_a,info)   ! r_i * r_i
                else
                    r_scalar_product    = r_scalar_product_next
                end if
                
                ! call prec%apply(r_single,z,desc_a,info,work=aux)

                call psb_spmm(sone,a_single,d_single,szero,rho_single,desc_a,info) ! rho_i = A * d_i                

                partial_result_d_rho      = psb_gedot(d_single,rho_single,desc_a,info) ! d_i * rho_i
                alpha                     = r_scalar_product / partial_result_d_rho


                call psb_geaxpby(alpha,d_single,sone,x_single,desc_a,info)     ! x_i+1 = x_i + alpha_i * d_i

                ! r_i+1 = r_i - alpha_i * rho_i
                do i=1,size(r_double%v%v) 
                    r_double%v%v(i)            = (r_single%v%v(i) * 1.d0) - ( (alpha * 1.d0) * ( rho_single%v%v(i) *1.d0 ) )
                    r_single%v%v(i)            = r_double%v%V(i)
                end do

                r_scalar_product_next = psb_gedot(r_single,r_single,desc_a,info)      ! r_i+1 * r_i+1
                beta = r_scalar_product_next / r_scalar_product


                 ! d_i+1 = r_i+1 + beta_i+1 * d_i     
                do i=1,size(d_double%v%v)
                    d_double%v%v(i)            = (r_single%v%v(i) * 1.d0) + ( (beta * 1.d0) * d_double%v%v(i) )
                    d_single%v%v(i)            = d_double%v%V(i)
                end do
                

                ! ||r|| / ||b||
                r_norm = psb_norm2(r_double, desc_a, info)

                err = r_norm / b_norm
                write(13,*) "---> ",itx, err, error_stopping_criterion, r_norm

                if(err < error_stopping_criterion) then
                    if(stagnation.eqv..true.) exit restart
                    stagnation = .true.

                    exit iteration
                end if
                stagnation = .false.
        
            end do iteration
        end do restart
        
        
    
        if(my_rank == psb_root_) then
            iter = itx
        end if
        
        call psb_bcast(ctxt, iter)
    
        
    
        call psb_erractionrestore(err_act)


        return
        
    9999 call psb_error_handler(err_act)

        return
        
        end subroutine psb_dscg_1_impl

end module psb_ds_cg_1