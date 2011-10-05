!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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

module psb_d_prec_mod
  use psb_d_prec_type

  interface psb_precbld
    subroutine psb_dprecbld(a,desc_a,prec,info,upd,mold,afmt)
      use psb_base_mod, only  : psb_desc_type, psb_dspmat_type,&
           & psb_d_base_sparse_mat, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_dprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
      character(len=*), intent(in), optional    :: afmt
      class(psb_d_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_dprecbld
  end interface

  interface psb_precinit
    subroutine psb_dprecinit(prec,ptype,info)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
    end subroutine psb_dprecinit
  end interface

  interface psb_precset
    subroutine psb_dprecseti(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      integer                                :: what, val 
      integer, intent(out)                   :: info
    end subroutine psb_dprecseti
    subroutine psb_dprecsetd(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      integer                                :: what
      real(psb_dpk_)                       :: val 
      integer, intent(out)                   :: info
    end subroutine psb_dprecsetd
  end interface

  interface psb_ilu_fct
    subroutine psb_dilu_fct(a,l,u,d,info,blck)
      use psb_base_mod, only  : psb_desc_type, psb_dspmat_type,&
           & psb_d_csr_sparse_mat, psb_dpk_
      integer, intent(out)                ::     info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_d_csr_sparse_mat),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine psb_dilu_fct
  end interface


end module psb_d_prec_mod