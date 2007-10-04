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
! File: psb_zgelp.f90
!
!
! Subroutine: psb_zgelp
!             Apply a left permutation to a dense matrix
!
! Arguments:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:,:).
! info     - integer.                 Return code.
subroutine psb_zgelp(trans,iperm,x,info)
  use psb_serial_mod, psb_protect_name => psb_zgelp
  use psb_const_mod
  use psb_error_mod
  implicit none

  complex(kind(1.d0)), intent(inout)      ::  x(:,:)
  integer, intent(in)                  ::  iperm(:)
  integer, intent(out)                 ::  info
  character, intent(in)                :: trans

  ! local variables
  integer                  :: ictxt,np,me,nrow,ncol
  complex(kind(1.d0)),allocatable ::  dtemp(:)
  integer, allocatable     :: itemp(:)
  integer                  :: int_err(5), i1sz, i2sz, i, err_act
  real(kind(1.d0)),parameter    :: one=1
  logical, parameter :: debug=.false.

  interface 
    subroutine zgelp(trans,m,n,p,b,ldb,work,lwork,ierror)
      integer, intent(in)  :: ldb, m, n, lwork
      integer, intent(out) :: ierror
      character, intent(in) :: trans
      complex(kind(1.d0)), intent(inout) ::  b(ldb,*), work(*)
      integer, intent(in)  :: p(*)
    end subroutine zgelp
  end interface

  interface isaperm

    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  character(len=20)   :: name, ch_err
  name = 'psb_zgelp'

  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  i1sz    = size(x,dim=1)
  i2sz    = size(x,dim=2)

  if (debug) write(*,*) 'gelp: ',i1sz,i2sz
  allocate(dtemp(i1sz),itemp(size(iperm)),stat=info)
  if (info /= 0) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if
  itemp(:) = iperm(:) 
  
  if (.not.isaperm(i1sz,itemp)) then
    info = 70
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif


  call zgelp(trans,i1sz,i2sz,itemp,x,i1sz,dtemp,i1sz,info)
  if(info.ne.0) then
    info=4010
    ch_err='zgelp'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  deallocate(dtemp,itemp)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_zgelp



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
!
!  WARNING: This routine should be changed and moved to the serial part 
!           i.e. taking out the descriptor. 
!
!
! Subroutine: psb_zgelpv
!             Apply a left permutation to a dense matrix
!
! Arguments:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:).
! info     - integer.                 Return code.
subroutine psb_zgelpv(trans,iperm,x,info)
  use psb_serial_mod, psb_protect_name => psb_zgelpv
  use psb_const_mod
  use psb_error_mod
  implicit none

  complex(kind(1.d0)), intent(inout) ::  x(:)
  integer, intent(in)                  ::  iperm(:)
  integer, intent(out)                 ::  info
  character, intent(in)              ::  trans

  ! local variables
  integer :: ictxt,np,me
  integer :: int_err(5), i1sz,nrow,ncol, i, err_act
  complex(kind(1.d0)),allocatable  ::  dtemp(:)
  integer, allocatable     :: itemp(:)
  real(kind(1.d0)),parameter    :: one=1
  logical, parameter :: debug=.false.

  interface zgelp
    subroutine zgelp(trans,m,n,p,b,ldb,work,lwork,ierror)
      integer, intent(in)  :: ldb, m, n, lwork
      integer, intent(out) :: ierror
      character, intent(in) :: trans
      complex(kind(1.d0)), intent(inout) ::  b(*), work(*)
      integer, intent(in)  :: p(*)
    end subroutine zgelp
  end interface

  interface isaperm

    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  character(len=20)   :: name, ch_err
  name = 'psb_zgelpv'

  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  i1sz = size(x)


  if (debug) write(*,*) 'gelp: ',i1sz
  allocate(dtemp(i1sz),itemp(size(iperm)),stat=info)
  if (info /= 0) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if
  itemp(:) = iperm(:) 
  
  if (.not.isaperm(i1sz,itemp)) then
    info = 70
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif


  call zgelp(trans,i1sz,1,itemp,x,i1sz,dtemp,i1sz,info)
  if(info.ne.0) then
    info=4010
    ch_err='zgelp'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  deallocate(dtemp)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_zgelpv
