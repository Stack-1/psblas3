subroutine psb_set_coher(ictxt,isvch)
  integer :: ictxt, isvch
  ! Ensure global coherence for convergence checks.
  Call blacs_get(ictxt,16,isvch)
  Call blacs_set(ictxt,16,1)
end subroutine psb_set_coher
subroutine psb_restore_coher(ictxt,isvch)
  integer :: ictxt, isvch
  ! Ensure global coherence for convergence checks.
  Call blacs_set(ictxt,16,isvch)
end subroutine psb_restore_coher
subroutine psb_get_mpicomm(ictxt,comm)
  integer :: ictxt, comm
  call blacs_get(ictxt,10,comm)
end subroutine psb_get_mpicomm
subroutine psb_get_rank(rank,ictxt,id)
  integer :: rank,ictxt, id
  integer :: blacs_pnum
  rank =  blacs_pnum(ictxt,id,0)
end subroutine psb_get_rank
