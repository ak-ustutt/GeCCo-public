      integer ::
     &     irt_sched, irt_contr, len_str_block, len_cnt_block,
     &     force_batching, force_ooc_sort, maxbranch, spinadapt
      ! note: the sv_fix, sv_thresh, etc. options should be abandoned from here
      logical ::
     &     use_tr, sv_fix, sv_old, use_auto_opt, jac_fix
      real(8) ::
     &     sv_thresh, jac_thresh, tikhonov, reo_factor
      common/routes/
     &     irt_sched, irt_contr, len_str_block, len_cnt_block,
     &     force_batching, force_ooc_sort, maxbranch,
     &     use_tr,sv_fix,sv_old,sv_thresh,
     &     jac_thresh,tikhonov,use_auto_opt,
     &     jac_fix,reo_factor,spinadapt
