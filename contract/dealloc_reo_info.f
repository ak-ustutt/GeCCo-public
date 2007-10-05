*----------------------------------------------------------------------*
      subroutine dealloc_reo_info(reo_info)
*----------------------------------------------------------------------*
*     initialize all pointers and set all mxvtx, mxarc, mxfac variables
*     to zero (describing the current size of allocated space for 
*     vertices, arcs, and factorization info)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_reorder_info.h'

      type(reorder_info), intent(inout) ::
     &     reo_info

      reo_info%nreo = 0
      reo_info%n_op_reo = 0

      if (associated(reo_info%reo)) then
        deallocate(reo_info%reo)
        reo_info%reo => null()
      end if
      
      if (associated(reo_info%merge_stp1)) then
        deallocate(reo_info%merge_stp1)
        reo_info%merge_stp1 => null()
      end if

      if (associated(reo_info%merge_stp2)) then
        deallocate(reo_info%merge_stp2)
        reo_info%merge_stp2 => null()
      end if

      if (associated(reo_info%iocc_reo)) then
        deallocate(reo_info%iocc_reo)
        reo_info%iocc_reo => null()
      end if

      if (associated(reo_info%iocc_opreo0)) then
        deallocate(reo_info%iocc_opreo0)
        reo_info%iocc_opreo0 => null()
      end if

      if (associated(reo_info%map_reo1c)) then
        deallocate(reo_info%map_reo1c)
        reo_info%map_reo1c => null()
      end if

      if (associated(reo_info%map_reo1a)) then
        deallocate(reo_info%map_reo1a)
        reo_info%map_reo1a => null()
      end if

      if (associated(reo_info%map_reo2c)) then
        deallocate(reo_info%map_reo2c)
        reo_info%map_reo2c => null()
      end if

      if (associated(reo_info%map_reo2a)) then
        deallocate(reo_info%map_reo2a)
        reo_info%map_reo2a => null()
      end if

      if (associated(reo_info%cinfo_opreo0c)) then
        deallocate(reo_info%cinfo_opreo0c)
        reo_info%cinfo_opreo0c => null()
      end if

      if (associated(reo_info%cinfo_opreo0a)) then
        deallocate(reo_info%cinfo_opreo0a)
        reo_info%cinfo_opreo0a => null()
      end if

      if (associated(reo_info%cinfo_reo_c)) then
        deallocate(reo_info%cinfo_reo_c)
        reo_info%cinfo_reo_c => null()
      end if

      if (associated(reo_info%cinfo_reo_a)) then
        deallocate(reo_info%cinfo_reo_a)
        reo_info%cinfo_reo_a => null()
      end if

      return
      end
