      real(8) ::
     &     cnt_op1op2(3), cnt_kernel(3), cnt_dloop(3), 
     &     cnt_wr(3), cnt_rd(3),
     &     cnt_coll1(2), cnt_coll2(2), cnt_dgemm(2), cnt_scatt(2),
     &     cnt_reo(2)
     &     ,cnt_test(12)
      logical ::
     &     cnt_used_reo
      integer(8) ::
     &     mm_call, mm_dim1, mm_dim1sq, mm_dim2, mm_dim2sq,
     &     mm_cnt, mm_cntsq, cnt_maxscr
      common/contr_times/ cnt_op1op2, cnt_kernel, cnt_dloop,
     &     cnt_wr, cnt_rd, cnt_coll1, cnt_coll2, cnt_dgemm, cnt_scatt,
     &     cnt_reo
     &     ,cnt_test,
     &     mm_call, mm_dim1, mm_dim1sq, mm_dim2, mm_dim2sq,
     &     mm_cnt, mm_cntsq, cnt_maxscr,
     &     cnt_used_reo
      integer(8) ::
     &     st_total,st_ex1c_z,st_ex1a_z,st_ex2c_z,st_ex2a_z,
     &     st_cntc_z,st_cnta_z,st_noop1_reo,st_noop2_reo,
     &     st_noop12_reo
      common/contr_stat/
     &     st_total,st_ex1c_z,st_ex1a_z,st_ex2c_z,st_ex2a_z,
     &     st_cntc_z,st_cnta_z,st_noop1_reo,st_noop2_reo,
     &     st_noop12_reo
    
