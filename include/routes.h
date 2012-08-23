      integer ::
     &     irt_sched, irt_contr, len_str_block, len_cnt_block,
     &     force_batching, force_ooc_sort, maxbranch
      logical ::
     &     use_tr, sv_fix
      real(8) ::
     &     sv_thresh, tikhonov
      common/routes/
     &     irt_sched, irt_contr, len_str_block, len_cnt_block,
     &     force_batching, force_ooc_sort, maxbranch,
     &     use_tr,sv_fix,sv_thresh,tikhonov
