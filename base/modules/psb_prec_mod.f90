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

module psb_prec_mod
  use psb_prec_type


  interface psb_precbld
    subroutine psb_dprecbld(a,desc_a,prec,info,upd)
      use psb_descriptor_type
      use psb_prec_type
      implicit none
      type(psb_dspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_dprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_dprecbld
    subroutine psb_zprecbld(a,desc_a,prec,info,upd)
      use psb_descriptor_type
      use psb_prec_type
      implicit none
      type(psb_zspmat_type), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_zprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_zprecbld
  end interface

  interface psb_precset
    subroutine psb_dprecset(prec,ptype,info,iv,rs,rv,ilev,nlev)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: iv(:)
      integer, optional, intent(in)          :: nlev,ilev
      real(kind(1.d0)), optional, intent(in) :: rs
      real(kind(1.d0)), optional, intent(in) :: rv(:)
    end subroutine psb_dprecset
    subroutine psb_zprecset(prec,ptype,info,iv,rs,rv,ilev,nlev)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      implicit none
      type(psb_zprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
      integer, optional, intent(in)          :: iv(:)
      real(kind(1.d0)), optional, intent(in) :: rs
      real(kind(1.d0)), optional, intent(in) :: rv(:)
      integer, optional, intent(in)          :: nlev,ilev
    end subroutine psb_zprecset
  end interface


  interface psb_precfree
    subroutine psb_dprecfree(p,info)
      use psb_descriptor_type
      use psb_serial_mod
      use psb_const_mod
      use psb_prec_type
      type(psb_dprec_type), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine psb_dprecfree
    subroutine psb_zprecfree(p,info)
      use psb_descriptor_type
      use psb_serial_mod
      use psb_const_mod
      use psb_prec_type
      type(psb_zprec_type), intent(inout) :: p
      integer, intent(out)                :: info
    end subroutine psb_zprecfree
  end interface

  interface psb_prc_aply
    subroutine psb_dprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_dprec_type), intent(in)  :: prec
      real(kind(0.d0)),intent(inout)    :: x(:), y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(kind(0.d0)),intent(inout), optional, target :: work(:)
    end subroutine psb_dprc_aply
    subroutine psb_dprc_aply1(prec,x,desc_data,info,trans)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_dprec_type), intent(in)  :: prec
      real(kind(0.d0)),intent(inout)    :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_dprc_aply1
    subroutine psb_zprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_zprec_type), intent(in)  :: prec
      complex(kind(0.d0)),intent(inout) :: x(:), y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      complex(kind(0.d0)),intent(inout), optional, target :: work(:)
    end subroutine psb_zprc_aply
    subroutine psb_zprc_aply1(prec,x,desc_data,info,trans)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_zprec_type), intent(in)  :: prec
      complex(kind(0.d0)),intent(inout) :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_zprc_aply1
  end interface

  interface psb_baseprc_bld
    subroutine psb_dbaseprc_bld(a,desc_a,p,info,upd)
      Use psb_spmat_type
      use psb_descriptor_type
      use psb_prec_type
      type(psb_dspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_dbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine psb_dbaseprc_bld
    subroutine psb_zbaseprc_bld(a,desc_a,p,info,upd)
      Use psb_spmat_type
      use psb_descriptor_type
      use psb_prec_type
      type(psb_zspmat_type), target              :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_zbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine psb_zbaseprc_bld
  end interface

  interface psb_mlprc_bld
    subroutine psb_dmlprc_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(psb_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                      :: info
    end subroutine psb_dmlprc_bld
    subroutine psb_zmlprc_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      type(psb_zbaseprc_type), intent(inout),target :: p
      integer, intent(out)                      :: info
    end subroutine psb_zmlprc_bld
  end interface


  interface psb_baseprc_aply
    subroutine psb_dbaseprc_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_descriptor_type
      use psb_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(psb_dbaseprc_type), intent(in) :: prec
      real(kind(0.d0)),intent(inout)      :: x(:), y(:)
      real(kind(0.d0)),intent(in)         :: alpha,beta
      character(len=1)                    :: trans
      real(kind(0.d0)),target             :: work(:)
      integer, intent(out)                :: info
    end subroutine psb_dbaseprc_aply

    subroutine psb_zbaseprc_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
      use psb_descriptor_type
      use psb_prec_type
      type(psb_desc_type),intent(in)      :: desc_data
      type(psb_zbaseprc_type), intent(in) :: prec
      complex(kind(1.d0)),intent(inout)   :: x(:), y(:)
      complex(kind(1.d0)),intent(in)      :: alpha,beta
      character(len=1)                    :: trans
      complex(kind(1.d0)),target          :: work(:)
      integer, intent(out)                :: info
    end subroutine psb_zbaseprc_aply
  end interface

  interface psb_mlprc_aply
     subroutine psb_dmlprc_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)      :: desc_data
       type(psb_dbaseprc_type), intent(in) :: baseprecv(:)
       real(kind(0.d0)),intent(in)         :: alpha,beta
       real(kind(0.d0)),intent(inout)      :: x(:), y(:)
       character                           :: trans
       real(kind(0.d0)),target             :: work(:)
       integer, intent(out)                :: info
     end subroutine psb_dmlprc_aply
     subroutine psb_zmlprc_aply(alpha,baseprecv,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type),intent(in)      :: desc_data
       type(psb_zbaseprc_type), intent(in) :: baseprecv(:)
       complex(kind(0.d0)),intent(in)      :: alpha,beta
       complex(kind(0.d0)),intent(inout)   :: x(:), y(:)
       character                           :: trans
       complex(kind(0.d0)),target          :: work(:)
       integer, intent(out)                :: info
     end subroutine psb_zmlprc_aply
  end interface

  interface psb_bjac_aply
     subroutine psb_dbjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type), intent(in)       :: desc_data
       type(psb_dbaseprc_type), intent(in)   :: prec
       real(kind(0.d0)),intent(inout)        :: x(:), y(:)
       real(kind(0.d0)),intent(in)           :: alpha,beta
       character(len=1)                      :: trans
       real(kind(0.d0)),target               :: work(:)
       integer, intent(out)                  :: info
     end subroutine psb_dbjac_aply

     subroutine psb_zbjac_aply(alpha,prec,x,beta,y,desc_data,trans,work,info)
       use psb_descriptor_type
       use psb_prec_type
       type(psb_desc_type), intent(in)       :: desc_data
       type(psb_zbaseprc_type), intent(in)   :: prec
       complex(kind(0.d0)),intent(inout)        :: x(:), y(:)
       complex(kind(0.d0)),intent(in)           :: alpha,beta
       character(len=1)                      :: trans
       complex(kind(0.d0)),target               :: work(:)
       integer, intent(out)                  :: info
     end subroutine psb_zbjac_aply
  end interface
  

  interface psb_diagsc_bld
    subroutine psb_ddiagsc_bld(a,desc_data,p,upd,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(psb_dbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
    end subroutine psb_ddiagsc_bld
    subroutine psb_zdiagsc_bld(a,desc_data,p,upd,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      integer, intent(out) :: info
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(psb_zbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
    end subroutine psb_zdiagsc_bld
  end interface

  interface psb_ilu_bld
    subroutine psb_dilu_bld(a,desc_data,p,upd,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      integer, intent(out) :: info
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(psb_dbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
    end subroutine psb_dilu_bld
    subroutine psb_zilu_bld(a,desc_data,p,upd,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      integer, intent(out) :: info
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type),intent(in)            :: desc_data
      type(psb_zbaseprc_type), intent(inout)    :: p
      character, intent(in)                     :: upd
    end subroutine psb_zilu_bld
  end interface

  interface psb_slu_bld
    subroutine psb_dslu_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_dspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_dslu_bld
    subroutine psb_zslu_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_zspmat_type), intent(inout)   :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_zbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_zslu_bld
  end interface

  interface psb_umf_bld
    subroutine psb_dumf_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_dspmat_type), intent(inout)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_dumf_bld
    subroutine psb_zumf_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_zspmat_type), intent(in)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_zbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_zumf_bld
  end interface


  interface psb_ilu_fct
    subroutine psb_dilu_fct(a,l,u,d,info,blck)
      use psb_spmat_type
      integer, intent(out)                ::     info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine psb_dilu_fct
    subroutine psb_zilu_fct(a,l,u,d,info,blck)
      use psb_spmat_type
      integer, intent(out)                ::     info
      type(psb_zspmat_type),intent(in)    :: a
      type(psb_zspmat_type),intent(inout) :: l,u
      type(psb_zspmat_type),intent(in), optional, target :: blck
      complex(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine psb_zilu_fct
  end interface

  interface psb_as_matbld
    Subroutine psb_dasmatbld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_serial_mod
      Use psb_descriptor_type
      Use psb_prec_type
      integer, intent(in)                  :: ptype,novr
      Type(psb_dspmat_type), Intent(in)    ::  a
      Type(psb_dspmat_type), Intent(inout) ::  blk
      Type(psb_desc_type), Intent(inout)   :: desc_p
      Type(psb_desc_type), Intent(in)      :: desc_data 
      Character, Intent(in)                :: upd
      integer, intent(out)                 :: info
      character(len=5), optional           :: outfmt
    end Subroutine psb_dasmatbld
    Subroutine psb_zasmatbld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_serial_mod
      Use psb_descriptor_type
      Use psb_prec_type
      integer, intent(in)                  :: ptype,novr
      Type(psb_zspmat_type), Intent(in)    ::  a
      Type(psb_zspmat_type), Intent(inout) ::  blk
      Type(psb_desc_type), Intent(inout)   :: desc_p
      Type(psb_desc_type), Intent(in)      :: desc_data 
      Character, Intent(in)                :: upd
      integer, intent(out)                 :: info
      character(len=5), optional           :: outfmt
    end Subroutine psb_zasmatbld
  end interface

  interface psb_sp_renum
    subroutine psb_dsp_renum(a,desc_a,blck,p,atmp,info)
      use psb_prec_type
      use psb_descriptor_type
      use psb_spmat_type
      implicit none

      !     .. array Arguments ..                                                     
      type(psb_dspmat_type), intent(in)      :: a,blck
      type(psb_dspmat_type), intent(inout)   :: atmp
      type(psb_dbaseprc_type), intent(inout) :: p
      type(psb_desc_type), intent(in)        :: desc_a
      integer, intent(out)   :: info
    end subroutine psb_dsp_renum
    subroutine psb_zsp_renum(a,desc_a,blck,p,atmp,info)
      use psb_prec_type
      use psb_descriptor_type
      use psb_spmat_type
      implicit none

      !     .. array Arguments ..                                                     
      type(psb_zspmat_type), intent(in)      :: a,blck
      type(psb_zspmat_type), intent(inout)   :: atmp
      type(psb_zbaseprc_type), intent(inout) :: p
      type(psb_desc_type), intent(in)        :: desc_a
      integer, intent(out)   :: info
    end subroutine psb_zsp_renum
  end interface


  interface psb_genaggrmap
    subroutine psb_dgenaggrmap(aggr_type,a,desc_a,nlaggr,ilaggr,info)
      use psb_spmat_type
      use psb_descriptor_type
      implicit none
      integer, intent(in)               :: aggr_type
      type(psb_dspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable              :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine psb_dgenaggrmap
    subroutine psb_zgenaggrmap(aggr_type,a,desc_a,nlaggr,ilaggr,info)
      use psb_spmat_type
      use psb_descriptor_type
      implicit none
      integer, intent(in)               :: aggr_type
      type(psb_zspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer, allocatable              :: ilaggr(:),nlaggr(:)
      integer, intent(out)              :: info
    end subroutine psb_zgenaggrmap
  end interface

  interface psb_bldaggrmat
    subroutine psb_dbldaggrmat(a,desc_a,ac,desc_ac,p,info)
      use psb_prec_type
      use psb_descriptor_type
      use psb_spmat_type
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in)           :: desc_a
      type(psb_dspmat_type), intent(inout),target :: ac
      type(psb_desc_type), intent(inout)        :: desc_ac
      type(psb_dbaseprc_type), intent(inout), target :: p
      integer, intent(out)                      :: info
    end subroutine psb_dbldaggrmat
    subroutine psb_zbldaggrmat(a,desc_a,ac,desc_ac,p,info)
      use psb_prec_type
      use psb_descriptor_type
      use psb_spmat_type
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_zbaseprc_type), intent(inout),target    :: p
      type(psb_zspmat_type), intent(inout),target :: ac
      type(psb_desc_type), intent(in)           :: desc_a
      type(psb_desc_type), intent(inout)        :: desc_ac
      integer, intent(out)                      :: info
    end subroutine psb_zbldaggrmat
  end interface

end module psb_prec_mod