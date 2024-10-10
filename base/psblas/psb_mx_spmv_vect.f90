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
! File: psb_mx_spmv_vect.f90
!
!
! Subroutine: psb_mx_spmv_vect
!     Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A * Pc * X  + beta * Y,  or
!
!     Y := alpha * Pr * A' * Pr * X  + beta * Y,
!
!  alpha and beta are scalars, X and Y are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Arguments:   
!    alpha   -  real                The scalar alpha.
!    a       -  type(psb_sspmat_type). The sparse matrix containing A.
!    x       -  type(psb_s_vect_type) The input vector containing the entries of ( X ).
!    beta    -  real                The scalar beta.
!    y       -  type(psb_d_vect_type) The input vector containing the entries of ( Y ).
!    desc_a  -  type(psb_desc_type).   The communication descriptor.
!    info    -  integer.               Return code
!    trans   -  character(optional).   Whether A or A'. Default:  'N' 
!    work(:) -  real,(optional).    Working area.
!    doswap  -  logical(optional).     Whether to performe halo updates.
! 
subroutine psb_mx_spmv_vect(alpha,a,x,beta,y,desc_a,info,&
    & trans, work, doswap)   
    use psb_base_mod, psb_protect_name => psb_mx_spmv_vect
    use psi_mod
    implicit none
    
    ! Function parameters
    real(psb_spk_), intent(in)                          :: alpha, beta
    type(psb_s_vect_type), intent(inout)                :: x
    type(psb_d_vect_type), intent(inout)                :: y
    type(psb_sspmat_type), intent(in)                   :: a
    type(psb_desc_type), intent(in)                     :: desc_a
    integer(psb_ipk_), intent(out)                      :: info
    real(psb_spk_), optional, target, intent(inout)     :: work(:)
    character, intent(in), optional                     :: trans
    logical, intent(in), optional                       :: doswap

    ! local variables
    type(psb_ctxt_type)                                 :: ctxt
    integer(psb_ipk_)                                   :: np, me,&
         & err_act, iix, jjx, iia, jja,  nrow, ncol, lldx, lldy, &
         & liwork, iiy, jjy, ib, ip, idx
    integer(psb_lpk_) :: ix, ijx, iy, ijy, m, n, ia, ja
    integer(psb_ipk_), parameter                        :: nb=4
    real(psb_spk_), pointer                             :: x_iwork(:), xp(:)
    real(psb_dpk_), pointer                             :: y_iwork(:), yp(:)

    real(psb_spk_), allocatable                         :: xvsave(:)
    character                                           :: trans_
    character(len=20)                                   :: name, ch_err
    logical                                             :: aliw, doswap_
    integer(psb_ipk_)                                   :: debug_level, debug_unit
    logical, parameter                                  :: do_timings=.true.
    integer(psb_ipk_), save                             :: mv_phase1=-1, mv_phase2=-1, mv_phase3=-1, mv_phase4=-1
    integer(psb_ipk_), save                             :: mv_phase11=-1, mv_phase12=-1

    name='psb_mx_spmv_vect'
    info=psb_success_
    call psb_erractionsave(err_act)
    if  (psb_errstatus_fatal()) then
      info = psb_err_internal_error_ ;    goto 9999
    end if
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ctxt=desc_a%get_context()
    call psb_info(ctxt, me, np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (.not.allocated(y%v)) then 
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (present(doswap)) then
      doswap_ = doswap
    else
      doswap_ = .true.
    endif

    if (present(trans)) then     
      trans_ = psb_toupper(trans)
    else
      trans_ = 'N'
    endif
    if ( (trans_ == 'N').or.(trans_ == 'T')&
         & .or.(trans_ == 'C')) then
    else
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name)
      goto 9999
    end if
    if ((do_timings).and.(mv_phase1==-1))       &
         & mv_phase1 = psb_get_timer_idx("SPMM: and send ")
    if ((do_timings).and.(mv_phase2==-1))       &
         & mv_phase2 = psb_get_timer_idx("SPMM: and cmp ad")
    if ((do_timings).and.(mv_phase3==-1))       &
         & mv_phase3 = psb_get_timer_idx("SPMM: and rcv")
    if ((do_timings).and.(mv_phase4==-1))       &
         & mv_phase4 = psb_get_timer_idx("SPMM: and cmp and")
    if ((do_timings).and.(mv_phase11==-1))       &
         & mv_phase11 = psb_get_timer_idx("SPMM: noand exch ")
    if ((do_timings).and.(mv_phase12==-1))       &
         & mv_phase12 = psb_get_timer_idx("SPMM: noand cmp")


    m    = desc_a%get_global_rows()
    n    = desc_a%get_global_cols()
    nrow = desc_a%get_local_rows()
    ncol = desc_a%get_local_cols()
    lldx = x%get_nrows()
    lldy = y%get_nrows()
    if ((info == 0).and.(lldx<ncol)) call x%reall(ncol,info)
    if ((info == 0).and.(lldy<ncol)) call y%reall(ncol,info)

    if (psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      ch_err='reall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    x_iwork => null()
    y_iwork => null()

    ! check for presence/size of a work area
    liwork= 2*ncol

    if (present(work)) then
      if (size(work) >= liwork) then
        aliw =.false.
      else
        aliw=.true.
      endif
    else
      aliw=.true.
    end if

    if (aliw) then
        allocate(x_iwork(liwork),stat=info)
        if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='Allocate x_iwork'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
        end if
        allocate(y_iwork(liwork),stat=info)
        if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='Allocate y_iwork'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
        end if
    else
        x_iwork => work
        ! y_iwork => work
    endif

    if (debug_level >= psb_debug_comp_) &
         & write(debug_unit,*) me,' ',trim(name),' Allocated work ', info

    if (trans_ == 'N') then
      !  Matrix is not transposed
    
      if (allocated(a%ad)) then
        block
          logical, parameter :: do_timings=.true.
          real(psb_dpk_) :: t1, t2, t3, t4, t5
          !if (me==0) write(0,*) 'going for overlap ',a%ad%get_fmt(),' ',a%and%get_fmt()
          if (do_timings) call psb_barrier(ctxt)
          if (do_timings) call psb_tic(mv_phase1)
          if (doswap_) call psi_swapdata(psb_swap_send_,&
               & szero,x%v,desc_a,x_iwork,info,data=psb_comm_halo_)
          if (do_timings) call psb_toc(mv_phase1)
          if (do_timings) call psb_tic(mv_phase2)          
          call a%ad%spmm(alpha,x%v,beta,y%v,info)
          if (do_timings) call psb_tic(mv_phase3)
          if (doswap_) call psi_swapdata(psb_swap_recv_,&
               & szero,x%v,desc_a,x_iwork,info,data=psb_comm_halo_)
          if (do_timings) call psb_toc(mv_phase3)
          if (do_timings) call psb_tic(mv_phase4)          
          call a%and%spmm(alpha,x%v,sone,y%v,info)
          if (do_timings) call psb_toc(mv_phase4)
        end block
        
      else
        block
          logical, parameter :: do_timings=.true.
          real(psb_dpk_) :: t1, t2, t3, t4, t5
          if (do_timings) call psb_barrier(ctxt)
        
          if (do_timings) call psb_tic(mv_phase11)          
          if (doswap_) then
            call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
                 & szero,x%v,desc_a,x_iwork,info,data=psb_comm_halo_)
          end if
          if (do_timings) call psb_toc(mv_phase11)
          if (do_timings) call psb_tic(mv_phase12)          
          call psb_csmm(alpha,a,x,beta,y,info)
          if (do_timings) call psb_toc(mv_phase12)
        end block
      end if

      if(info /= psb_success_) then
        info = psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if

    else
      !  Matrix is transposed

      !
      ! Non-empty overlap, need a buffer to hold
      ! the entries updated with average operator.
      ! Why the average? because in this way they will contribute
      ! with a proper scale factor (1/np) to the overall product.
      ! 
      call psi_ovrl_save(x%v,xvsave,desc_a,info)
      if (info == psb_success_) call psi_ovrl_upd(x%v,desc_a,psb_avg_,info)

      if (beta /= szero) call y%set(dzero,nrow+1,ncol)
      !  local Matrix-vector product
      if (info == psb_success_) call psb_csmm(alpha,a,x,beta,y,info,trans=trans_)

      if (debug_level >= psb_debug_comp_) &
           & write(debug_unit,*) me,' ',trim(name),' csmm ', info

      if (info == psb_success_) call psi_ovrl_restore(x%v,xvsave,desc_a,info)
      if (info /= psb_success_) then
        info = psb_err_from_subroutine_
        ch_err='psb_csmm'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      if (doswap_) then
        call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
             & done,y%v,desc_a,y_iwork,info)
        if (info == psb_success_) call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
             & done,y%v,desc_a,y_iwork,info,data=psb_comm_ovr_)

        if (debug_level >= psb_debug_comp_) &
             & write(debug_unit,*) me,' ',trim(name),' swaptran ', info
        if(info /= psb_success_) then
          info = psb_err_from_subroutine_
          ch_err='PSI_SwapTran'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
      end if

    end if

    if (aliw) deallocate(x_iwork,stat=info)
    if (debug_level >= psb_debug_comp_) &
         & write(debug_unit,*) me,' ',trim(name),' deallocat ',aliw, info
    if(info /= psb_success_) then
      info = psb_err_from_subroutine_
      ch_err='Deallocate x_iwork'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    

    nullify(x_iwork)
    
    if (aliw) deallocate(y_iwork,stat=info)
    if (debug_level >= psb_debug_comp_) &
         & write(debug_unit,*) me,' ',trim(name),' deallocat ',aliw, info
    if(info /= psb_success_) then
      info = psb_err_from_subroutine_
      ch_err='Deallocate y_iwork'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    
    
    nullify(y_iwork)


    call psb_erractionrestore(err_act)
    if (debug_level >= psb_debug_comp_) then 
      call psb_barrier(ctxt)
      write(debug_unit,*) me,' ',trim(name),' Returning '
    endif
    return  

9999 call psb_error_handler(ctxt,err_act)

    return
end subroutine psb_mx_spmv_vect