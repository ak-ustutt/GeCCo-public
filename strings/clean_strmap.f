*----------------------------------------------------------------------*
      subroutine clean_strmap(strmap_info)
*----------------------------------------------------------------------*
*     clean up everything
*----------------------------------------------------------------------*

      implicit none
      
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'

      type(strmapinf), intent(inout) ::
     &     strmap_info

      deallocate(strmap_info%offsets)
      deallocate(strmap_info%maxlen_blk)
      call file_delete(strmap_info%ffstrmap)

      call mem_clean_vbuffer('bfstrmap')

      return
      end

