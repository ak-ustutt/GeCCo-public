      type contraction_info

      integer ::
     &     ncblk_op1, nablk_op1,
     &     ncblk_op2, nablk_op2,
     &     ncblk_op1op2, nablk_op1op2,
     &     ncblk_op1op2tmp, nablk_op1op2tmp,
     &     ncblk_ex1, nablk_ex1,
     &     ncblk_ex2, nablk_ex2,
     &     ncblk_cnt, nablk_cnt

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
