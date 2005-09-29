
module psb_const_mod

  integer, parameter   :: psb_nohalo_=0,                 psb_halo_=4
  integer, parameter   :: psb_none_=0,                   psb_sum_=1
  integer, parameter   :: psb_avg_=2,                    psb_square_root_=3
  integer, parameter   :: psb_swap_send_=1,              psb_swap_recv_=2
  integer, parameter   :: psb_swap_sync_=4,              psb_swap_mpi_=8
  integer, parameter   :: psb_deadlock_check_=0,         psb_local_mtrx_check_=1
  integer, parameter   :: psb_local_comm_check_=2,       psb_consistency_check_=3
  integer, parameter   :: psb_global_check_=4,           psb_order_communication_=5
  integer, parameter   :: psb_change_represent_=6,       psb_loc_to_glob_check_=7
  integer, parameter   :: psb_convert_halo_=1,           psb_convert_ovrlap_=2
  integer, parameter   :: psb_act_ret_=0,                psb_act_abort_=1, no_err_=0
  integer, parameter   :: psb_dec_type_=1,               psb_m_=2,psb_n_=3
  integer, parameter   :: psb_n_row_=4,                  psb_n_col_=5,psb_ctxt_=6
  integer, parameter   :: psb_loc_to_glob_=7,            psb_mpi_c_=9,psb_mdata_size_=10
  integer, parameter   :: psb_desc_asb_=3099,            psb_desc_bld_=psb_desc_asb_+1
  integer, parameter   :: psb_desc_repl_=3199
  integer, parameter   :: psb_desc_upd_=psb_desc_bld_+1, psb_desc_upd_asb_=psb_desc_upd_+1
  integer, parameter   :: psb_upd_glb_=998,              psb_upd_loc_=997
  integer, parameter   :: psb_proc_id_=0,                psb_n_elem_recv_=1
  integer, parameter   :: psb_elem_recv_=2,              psb_n_elem_send_=2
  integer, parameter   :: psb_elem_send_=3,              psb_n_ovrlp_elem_=1
  integer, parameter   :: psb_ovrlp_elem_to_=2,          psb_ovrlp_elem_=0, psb_n_dom_ovr_=1
  integer, parameter   :: psb_nnz_=1
  integer, parameter   :: psb_no_comm_=-1
  integer, parameter   :: ione=1,izero=0,mone=-1
  integer, parameter   :: itwo=2, ithree=3,              psb_root_=0
  integer, parameter   :: psb_nztotreq_=1,               psb_nzrowreq_=2, psb_nzsizereq_=3
  integer, parameter   :: psb_del_bnd_=6,                psb_srtd_=7 
  integer, parameter   :: psb_state_=8,                  psb_upd_=9
  integer, parameter   :: psb_upd_pnt_=10,               psb_ifasize_=10
  integer, parameter   :: psb_spmat_null_=0,             psb_spmat_bld_=1
  integer, parameter   :: psb_spmat_asb_=2,              psb_spmat_upd_=4
  integer, parameter   :: psb_ireg_flgs_=10,             psb_ip2_=0
  integer, parameter   :: psb_iflag_=2,                  psb_ichk_=3
  integer, parameter   :: psb_nnzt_=4,                   psb_zero_=5,psb_ipc_=6
  integer, parameter   :: psb_perm_update_=98765,        psb_isrtdcoo_=98764
  integer, parameter   :: psb_maxjdrows_=8,              psb_minjdrows_=4
  integer, parameter   :: psb_dbleint_=2
  integer, parameter   :: psb_comm_halo_=0,              psb_comm_ovr_=1

  real(kind(1.d0)), parameter :: psb_colrow_=0.33, psb_percent_=0.7
  real(kind(1.d0)), parameter :: dzero=0.d0, done=1.d0

  character, parameter :: psb_all_='A',                  psb_topdef_=' '
  character(len=5)     :: psb_fidef_='CSR'
    

                                              
end module psb_const_mod                          
