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
! File:  psb_dtransp.f90 
! Subroutine: 
! Arguments:

subroutine psb_dtransp(a,b,c,fmt)
  use psb_spmat_type
  use psb_tools_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_dtransp
  implicit none

  type(psb_dspmat_type), intent(inout) :: a
  type(psb_dspmat_type), intent(out)   :: b
  integer, optional          :: c
  character(len=*), optional :: fmt

  character(len=5)           :: fmt_
  integer  ::c_, info, nz 
  integer, allocatable  :: itmp(:)
  type(psb_dspmat_type)      :: tmp

  if (present(c)) then 
    c_=c
  else
    c_=1
  endif
  if (present(fmt)) then 
    fmt_ = toupper(fmt)
  else 
    fmt_='CSR'
  endif

  call psb_nullify_sp(b)
  
  call psb_spcnv(a,b,info,afmt='coo')
  
  if (info /= 0) then 
    write(0,*) 'transp: info from CSDP ',info
    return
  end if
  call psb_transfer(b%ia1,itmp,info)
  call psb_transfer(b%ia2,b%ia1,info)
  call psb_transfer(itmp,b%ia2,info)

  b%m = a%k 
  b%k = a%m
  call psb_spcnv(b,info,afmt=fmt_)

  return
end subroutine psb_dtransp
