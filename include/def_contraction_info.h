      type contraction_info
      
      ! occupations in full form
c      integer ::
c     &     njoined_op1, njoined_op2, njoined_op1op2,
c     &     njoined_cnt
c
c      integer ::
c     &     iocc_ex1(:,:,:), iocc_ex2(:,:,:),
c     &     iocc_cnt(:,:,:),
c     &     iocc_op1(:,:,:), iocc_op2(:,:,:), iocc_op1op2(:,:,:),
c     &     irst_op1(:,:,:), irst_op2(:,:,:), irst_op1op2(:,:,:),
c      logical ::
c     &     tra_op1, tra_op2, tra_op1op2
c      integer ::
c     &     mst_op1, mst_op2, mst_op1op2,
c     &     gam_op1, gam_op2, gam_op1op2,
c      integer ::
c     &     merge_op1, merge_op2, merge_op1op2, merge_op2op1

      ! condensed occupation version
      integer ::
     &     ncblk_op1, nablk_op1,
     &     ncblk_op2, nablk_op2,
     &     ncblk_op1op2, nablk_op1op2,
     &     ncblk_op1op2tmp, nablk_op1op2tmp,
     &     ncblk_ex1, nablk_ex1,
     &     ncblk_ex2, nablk_ex2,
     &     ncblk_cnt, nablk_cnt,
     &     diag_type1, diag_type2, diag_type12

      integer, pointer ::
     &     cinfo_op1c(:,:), cinfo_op1a(:,:),
     &     cinfo_op2c(:,:), cinfo_op2a(:,:),
     &     cinfo_op1op2c(:,:), cinfo_op1op2a(:,:),
     &     cinfo_op1op2tmpc(:,:),cinfo_op1op2tmpa(:,:),
     &     cinfo_ex1c(:,:), cinfo_ex1a(:,:),
     &     cinfo_ex2c(:,:), cinfo_ex2a(:,:),
     &     cinfo_cntc(:,:), cinfo_cnta(:,:),
     &     map_info_1c(:), map_info_1a(:),
     &     map_info_2c(:), map_info_2a(:),
     &     map_info_12c(:), map_info_12a(:)

      end type contraction_info
