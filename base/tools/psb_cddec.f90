!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine psb_cddec(nloc, ictxt, desc_a, info)

  !  Purpose
  !  =======
  !  
  !  Allocate special descriptor for decoupled index space 
  ! 
  !
  !
  ! INPUT
  !======
  ! NLOC         : (Local Input) Integer 
  !                    Number of local indices 
  !                    required.
  !
  ! ictxt      : (Global Input)Integer BLACS context for an NPx1 grid 
  !                required.
  !
  ! OUTPUT
  !=========
  ! desc_a   : TYPEDESC
  ! desc_a OUTPUT FIELDS:
  !
  ! MATRIX_DATA   : Pointer to integer Array 
  !                contains some
  !               local and global information about matrix:
  !
  !  NOTATION        STORED IN		     EXPLANATION
  !  ------------ ---------------------- -------------------------------------
  !  DEC_TYPE        MATRIX_DATA[DEC_TYPE_]   Decomposition type, temporarly is
  !                      setted to 1( matrix not yet assembled)
  !  M 	             MATRIX_DATA[M_]          Total number of equations
  !  N 	             MATRIX_DATA[N_]          Total number of variables
  !  N_ROW           MATRIX_DATA[N_ROW_]      Number of local equations
  !  N_COL           MATRIX_DATA[N_COL_]      Number of local columns (see below)
  !  CTXT_A          MATRIX_DATA[CTXT_]     The BLACS context handle, 
  !                                         indicating
  !	  			            the global context of the operation
  !					    on the matrix.
  !					    The context itself is global.
  !
  !  GLOB_TO_LOC     Array of dimension equal to number of global 
  !                  rows/cols (MATRIX_DATA[M_]). On exit,
  !                  for all global indices either:
  !                  1. The index belongs to the current process; the entry
  !                     is set to the next free local row index.
  !                  2. The index belongs to process P (0<=P<=NP-1); the entry 
  !                     is set to 
  !                     -(NP+P+1)
  !
  !  LOC_TO_GLOB     An array of dimension equal to number of local cols N_COL
  !                  i.e. all columns of the matrix such that there is at least
  !                  one nonzero entry within the local row range. At the time 
  !                  this routine is called N_COL cannot be know, so we set 
  !                  N_COL=N_ROW, and dimension this vector on N_ROW plus an 
  !                  estimate. On exit the vector elements are set
  !                  to the index of the corresponding entry in GLOB_TO_LOC, or  
  !                  to -1 for indices I>N_ROW.
  !
  !
  !  HALO_INDEX      Not touched here, as it depends on the matrix pattern
  !
  !  OVRLAP_INDEX    On exit from this routine, the overlap indices are stored in
  !                  triples (Proc, 1, Index), similar to the assembled format 
  !                  but neither optimized, nor deadlock free. 
  !                  List is terminated with -1
  !
  !  OVRLAP_ELEM     On exit from this routine, just a list of pairs (index,#p).
  !                  List is terminated with -1.
  !                  
  !
  ! END OF desc_a OUTPUT FIELDS
  !
  !

  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psi_mod
  use psb_penv_mod
  implicit None
  !....Parameters...
  Integer, intent(in)               :: nloc,ictxt
  integer, intent(out)              :: info
  Type(psb_desc_type), intent(out)  :: desc_a

  !locals
  Integer              :: i,j,np,me,err,n,itmpov, k,&
       & l_ov_ix,l_ov_el,idx, err_act,m, ip,glx
  Integer              :: INT_ERR(5), thalo(1), tovr(1) 
  integer, allocatable :: nlv(:)
  logical, parameter   :: debug=.false.
  character(len=20)    :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  err=0
  name = 'psb_cddec'

  call psb_info(ictxt, me, np)
  if (debug) write(*,*) 'psb_cdalll: ',np,me

  n = nloc
  !... check nloc and n parameters....
  if (nloc < 1) then
    info = 10
    int_err(1) = 1
    int_err(2) = nloc
  else if (n < 1) then
    info = 10
    int_err(1) = 2
    int_err(2) = n
  endif

  if (info /= 0) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (debug) write(*,*) 'psb_cdall:  doing global checks'  
  !global check on m and n parameters

  allocate(nlv(0:np-1))
  nlv(:)  = 0
  nlv(me) = nloc

  call psb_sum(ictxt,nlv(1:np))
  m = sum(nlv)



  call psb_nullify_desc(desc_a)


  !count local rows number
  if ( m >psb_cd_get_large_threshold()) then 
    allocate(desc_a%loc_to_glob(nloc), desc_a%lprm(1),&
         & desc_a%ptree(2),desc_a%matrix_data(psb_mdata_size_),stat=info)  
    if (info == 0) call InitPairSearchTree(desc_a%ptree,info)
    if (info /= 0) then
      info=2025
      int_err(1)=nloc
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! set LOC_TO_GLOB array to all "-1" values
    desc_a%lprm(1) = 0
    desc_a%loc_to_glob(:) = -1
    desc_a%matrix_data(psb_n_row_)    = nloc
    desc_a%matrix_data(psb_n_col_)    = nloc
    desc_a%matrix_data(psb_m_)        = m
    desc_a%matrix_data(psb_n_)        = m
    desc_a%matrix_data(psb_dec_type_) = psb_desc_large_bld_
    desc_a%matrix_data(psb_ctxt_)     = ictxt
    call psb_get_mpicomm(ictxt,desc_a%matrix_data(psb_mpi_c_))

    do ip=0, np-1
      if (ip==me) then 
        do i=1, nlv(ip) 
          call SearchInsKeyVal(desc_a%ptree,j,i,glx,info)
          desc_a%loc_to_glob(i) = j
          j = j + 1
        enddo
      else
        do i=1, nlv(ip) 
          j = j + 1
        enddo
      endif
    enddo

    tovr  = -1 
    thalo = -1

    desc_a%lprm(:)         = 0

    call psi_cnv_dsc(thalo,tovr,desc_a,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_bld_cdesc')
      goto 9999
    end if


    desc_a%matrix_data(psb_dec_type_) = psb_desc_large_asb_
  else

    allocate(desc_a%glob_to_loc(m),desc_a%matrix_data(psb_mdata_size_),&
         &   desc_a%loc_to_glob(m),desc_a%lprm(1),stat=info)
    if (info /= 0) then     
      info=2025
      int_err(1)=m
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif


    desc_a%matrix_data(psb_n_row_)    = nloc
    desc_a%matrix_data(psb_n_col_)    = nloc
    desc_a%matrix_data(psb_m_)        = m
    desc_a%matrix_data(psb_n_)        = m
    desc_a%matrix_data(psb_dec_type_) = psb_desc_bld_
    desc_a%matrix_data(psb_ctxt_)     = ictxt
    call psb_get_mpicomm(ictxt,desc_a%matrix_data(psb_mpi_c_))

    j = 1
    do ip=0, np-1
      if (ip==me) then 
        do i=1, nlv(ip) 
          desc_a%glob_to_loc(j) = i
          desc_a%loc_to_glob(i) = j
          j = j + 1
        enddo
      else
        do i=1, nlv(ip) 
          desc_a%glob_to_loc(j) = -(np+ip+1)
          j = j + 1
        enddo
      endif
    enddo

    tovr  = -1 
    thalo = -1

    desc_a%lprm(:)         = 0

    desc_a%matrix_data(psb_ovl_state_) = psb_cd_ovl_bld_

    call psi_cnv_dsc(thalo,tovr,desc_a,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_bld_cdesc')
      goto 9999
    end if

    desc_a%matrix_data(psb_dec_type_) = psb_desc_asb_

  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == act_abort) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cddec