*----------------------------------------------------------------------*
      subroutine dealloc_cnt_info(cnt_info)
*----------------------------------------------------------------------*
*     initialize all pointers and set all mxvtx, mxarc, mxfac variables
*     to zero (describing the current size of allocated space for 
*     vertices, arcs, and factorization info)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction_info.h'

      type(contraction_info), intent(inout) ::
     &     cnt_info

      if (associated(cnt_info%cinfo_op1c)) then
        deallocate(cnt_info%cinfo_op1c)
        cnt_info%cinfo_op1c => null()
      end if
      if (associated(cnt_info%cinfo_op1a)) then
        deallocate(cnt_info%cinfo_op1a)
        cnt_info%cinfo_op1a => null()
      end if

      if (associated(cnt_info%cinfo_ex1c)) then
        deallocate(cnt_info%cinfo_ex1c)
        cnt_info%cinfo_ex1c => null()
      end if
      if (associated(cnt_info%cinfo_ex1a)) then
        deallocate(cnt_info%cinfo_ex1a)
        cnt_info%cinfo_ex1a => null()
      end if

      if (associated(cnt_info%cinfo_op2c)) then
        deallocate(cnt_info%cinfo_op2c)
        cnt_info%cinfo_op2c => null()
      end if
      if (associated(cnt_info%cinfo_op2a)) then
        deallocate(cnt_info%cinfo_op2a)
        cnt_info%cinfo_op2a => null()
      end if

      if (associated(cnt_info%cinfo_ex2c)) then
        deallocate(cnt_info%cinfo_ex2c)
        cnt_info%cinfo_ex2c => null()
      end if
      if (associated(cnt_info%cinfo_ex2a)) then
        deallocate(cnt_info%cinfo_ex2a)
        cnt_info%cinfo_ex2a => null()
      end if

      if (associated(cnt_info%cinfo_cntc)) then
        deallocate(cnt_info%cinfo_cntc)
        cnt_info%cinfo_cntc => null()
      end if
      if (associated(cnt_info%cinfo_cnta)) then
        deallocate(cnt_info%cinfo_cnta)
        cnt_info%cinfo_cnta => null()
      end if

      if (associated(cnt_info%cinfo_op1op2c)) then
        deallocate(cnt_info%cinfo_op1op2c)
        cnt_info%cinfo_op1op2c => null()
      end if
      if (associated(cnt_info%cinfo_op1op2a)) then
        deallocate(cnt_info%cinfo_op1op2a)
        cnt_info%cinfo_op1op2a => null()
      end if

      if (associated(cnt_info%cinfo_op1op2tmpc)) then
        deallocate(cnt_info%cinfo_op1op2tmpc)
        cnt_info%cinfo_op1op2tmpc => null()
      end if
      if (associated(cnt_info%cinfo_op1op2tmpa)) then
        deallocate(cnt_info%cinfo_op1op2tmpa)
        cnt_info%cinfo_op1op2tmpa => null()
      end if

      if (associated(cnt_info%map_info_1c)) then
        deallocate(cnt_info%map_info_1c)
        cnt_info%map_info_1c => null()
      end if
      if (associated(cnt_info%map_info_1a)) then
        deallocate(cnt_info%map_info_1a)
        cnt_info%map_info_1a => null()
      end if

      if (associated(cnt_info%map_info_2c)) then
        deallocate(cnt_info%map_info_2c)
        cnt_info%map_info_2c => null()
      end if
      if (associated(cnt_info%map_info_2a)) then
        deallocate(cnt_info%map_info_2a)
        cnt_info%map_info_2a => null()
      end if

      if (associated(cnt_info%map_info_12c)) then
        deallocate(cnt_info%map_info_12c)
        cnt_info%map_info_12c => null()
      end if
      if (associated(cnt_info%map_info_12a)) then
        deallocate(cnt_info%map_info_12a)
        cnt_info%map_info_12a => null()
      end if

      return
      end
