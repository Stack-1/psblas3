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
! File: psb_dspins.f90
!
! Subroutine: psb_dspins
!    Takes a cloud of points and inserts them into a sparse matrix.
! 
! Parameters: 
!    nz       - integer.                          The number of points to insert.
!    ia       - integer,dimension(:).             The row indices of the points.
!    ja       - integer,dimension(:).             The column indices of the points.
!    val      - real,dimension(:).                The values of the points to be inserted.
!    a        - type(<psb_dspmat_type>).          The sparse destination matrix.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Error code
!    rebuild  - logical                           Allows to reopen a matrix under
!                                                 certain circumstances.
!
subroutine psb_dspins(nz,ia,ja,val,a,desc_a,info,rebuild)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout)   :: desc_a
  type(psb_dspmat_type), intent(inout) :: a
  integer, intent(in)                  :: nz,ia(:),ja(:)
  real(kind(1.d0)), intent(in)         :: val(:)
  integer, intent(out)                 :: info
  logical, intent(in), optional        :: rebuild
  !locals.....

  integer :: i,icontxt,nprocs ,glob_row,row,k,start_row,end_row,&
       & first_loc_row,nrow,j, err,locix,locjx,err_act,&
       & dectype,mglob, nnza,m,n, pnt_halo,ncol, nh, ip, spstate
  integer,pointer        :: tia1(:),tia2(:), temp(:)
  integer                :: nprow,npcol, myrow ,mycol, iflag, isize, irlc
  logical, parameter     :: debug=.false.
  integer, parameter     :: relocsz=200
  logical                :: rebuild_

  interface psb_cdins
     subroutine psb_cdins(nz,ia,ja,desc_a,info)
       use psb_descriptor_type
       implicit none
       type(psb_desc_type), intent(inout) ::  desc_a
       integer, intent(in)            ::  nz,ia(:),ja(:)
       integer, intent(out)           :: info
     end subroutine psb_cdins
  end interface

  character(len=20)  :: name, ch_err

  info = 0
  name = 'psb_dspins'
  call psb_erractionsave(err_act)


  icontxt = desc_a%matrix_data(psb_ctxt_)
  dectype = desc_a%matrix_data(psb_dec_type_)
  mglob   = desc_a%matrix_data(psb_m_)
  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (npcol.ne.1) then
     info = 2030
     call psb_errpush(info,name)
     goto 9999
  endif
  if (.not.psb_is_ok_dec(dectype)) then 
     info = 3110
     call psb_errpush(info,name)
     goto 9999
  endif

  if (nz <= 0) then 
     info = 1111
     call psb_errpush(info,name)
     goto 9999
  end if
  if (size(ia) < nz) then 
     info = 1111
     call psb_errpush(info,name)
     goto 9999
  end if

  if (size(ja) < nz) then 
     info = 1111
     call psb_errpush(info,name)
     goto 9999
  end if
  if (size(val) < nz) then 
     info = 1111
     call psb_errpush(info,name)
     goto 9999
  end if

  if (present(rebuild)) then 
    rebuild_ = rebuild
  else
    rebuild_ = .false.
  endif

  spstate = a%infoa(psb_state_)
  if (psb_is_bld_dec(dectype)) then 
     call  psb_cdins(nz,ia,ja,desc_a,info)
     if (info /= 0) then
        info=4010
        ch_err='psb_cdins'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
     nrow = desc_a%matrix_data(psb_n_row_)
     ncol = desc_a%matrix_data(psb_n_col_)

     if (spstate == psb_spmat_bld_) then 
        call psb_coins(nz,ia,ja,val,a,desc_a%glob_to_loc,1,nrow,1,ncol,info)
        if (info /= 0) then
           info=4010
           ch_err='psb_coins'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     else
        info = 1123
        call psb_errpush(info,name)
        goto 9999
     end if
  else if (psb_is_asb_dec(dectype)) then 
     nrow = desc_a%matrix_data(psb_n_row_)
     ncol = desc_a%matrix_data(psb_n_col_)
     call psb_coins(nz,ia,ja,val,a,desc_a%glob_to_loc,1,nrow,1,ncol,&
          & info,rebuild=rebuild_)
     if (info /= 0) then
        info=4010
        ch_err='psb_coins'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
  else
     info = 1122
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dspins

