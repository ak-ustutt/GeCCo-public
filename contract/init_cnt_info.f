*----------------------------------------------------------------------*
      subroutine init_cnt_info(cnt_info,
     &     iocc_op1,iocc_ex1,nj_op1,iocc_op2,iocc_ex2,nj_op2,
     &     iocc_cnt,nj_cnt,
     &     iocc_op1op2,nj_op1op2,iocc_op1op2tmp,nj_op1op2tmp)
*----------------------------------------------------------------------*
*     initialize dimensions and arrays
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction_info.h'

      type(contraction_info), intent(inout) ::
     &     cnt_info

      integer, intent(in) ::
     &     nj_op1,nj_op2,nj_cnt,
     &     nj_op1op2,nj_op1op2tmp
      integer, intent(in) ::
     &     iocc_op1(ngastp,2,nj_op1),iocc_ex1(ngastp,2,nj_op1),
     &     iocc_op2(ngastp,2,nj_op2),iocc_ex2(ngastp,2,nj_op2),
     &     iocc_cnt(ngastp,2,nj_cnt),
     &     iocc_op1op2(ngastp,2,nj_op1op2),
     &     iocc_op1op2tmp(ngastp,2,nj_op1op2tmp)

      call get_num_subblk(cnt_info%ncblk_op1,cnt_info%nablk_op1,
     &                    iocc_op1,nj_op1)
      call get_num_subblk(cnt_info%ncblk_ex1,cnt_info%nablk_ex1,
     &                    iocc_ex1,nj_op1)
      call get_num_subblk(cnt_info%ncblk_op2,cnt_info%nablk_op2,
     &                    iocc_op2,nj_op2)
      call get_num_subblk(cnt_info%ncblk_ex2,cnt_info%nablk_ex2,
     &                    iocc_ex2,nj_op2)
      call get_num_subblk(cnt_info%ncblk_op1op2,cnt_info%nablk_op1op2,
     &                    iocc_op1op2,nj_op1op2)
      call get_num_subblk(cnt_info%ncblk_op1op2tmp,
     &                    cnt_info%nablk_op1op2tmp,
     &                    iocc_op1op2tmp,nj_op1op2tmp)
      call get_num_subblk(cnt_info%ncblk_cnt,cnt_info%nablk_cnt,
     &                    iocc_cnt,nj_cnt)
      
      allocate(
     &     cnt_info%cinfo_op1c(cnt_info%ncblk_op1,3),
     &     cnt_info%cinfo_op1a(cnt_info%nablk_op1,3),
     &     cnt_info%cinfo_ex1c(cnt_info%ncblk_ex1,3),
     &     cnt_info%cinfo_ex1a(cnt_info%nablk_ex1,3),
     &     cnt_info%cinfo_op2c(cnt_info%ncblk_op2,3),
     &     cnt_info%cinfo_op2a(cnt_info%nablk_op2,3),
     &     cnt_info%cinfo_ex2c(cnt_info%ncblk_ex2,3),
     &     cnt_info%cinfo_ex2a(cnt_info%nablk_ex2,3),
     &     cnt_info%cinfo_cntc(cnt_info%ncblk_cnt,3),
     &     cnt_info%cinfo_cnta(cnt_info%nablk_cnt,3),
     &     cnt_info%cinfo_op1op2c(cnt_info%ncblk_op1op2,3),
     &     cnt_info%cinfo_op1op2a(cnt_info%nablk_op1op2,3),
     &     cnt_info%cinfo_op1op2tmpc(cnt_info%ncblk_op1op2tmp,3),
     &     cnt_info%cinfo_op1op2tmpa(cnt_info%nablk_op1op2tmp,3)
     &     )
      allocate(
     &     cnt_info%map_info_1c(
     &                     max(1,cnt_info%ncblk_op1*2*(nj_op1+nj_cnt))),
     &     cnt_info%map_info_1a(
     &                     max(1,cnt_info%nablk_op1*2*(nj_op1+nj_cnt))),
     &     cnt_info%map_info_2c(
     &                     max(1,cnt_info%ncblk_op2*2*(nj_op2+nj_cnt))),
     &     cnt_info%map_info_2a(
     &                     max(1,cnt_info%nablk_op2*2*(nj_op2+nj_cnt))),
     &     cnt_info%map_info_12c(
     &               max(1,cnt_info%ncblk_op1op2tmp*2*(nj_op1+nj_op2))),
     &     cnt_info%map_info_12a(
     &               max(1,cnt_info%nablk_op1op2tmp*2*(nj_op1+nj_op2)))
     &     )

      return
      end
