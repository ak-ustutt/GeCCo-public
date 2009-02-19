      real(8) ::
     &     cnt_op1op2(3), cnt_kernel(3), cnt_dloop(3), 
     &     cnt_wr(3), cnt_rd(3),
     &     cnt_coll1(2), cnt_coll2(2), cnt_dgemm(2), cnt_scatt(2)
      common/contr_times/ cnt_op1op2, cnt_kernel, cnt_dloop,
     &     cnt_wr, cnt_rd, cnt_coll1, cnt_coll2, cnt_dgemm, cnt_scatt
