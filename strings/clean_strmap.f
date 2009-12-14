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

      integer ::
     &     idx

c      do idx = 1, strmap_info%mxgraph*strmap_info%mxgraph
c        if (associated(strmap_info%offsets(idx)%msms))
c     &      deallocate(strmap_info%offsets(idx)%msms)
c        if (associated(strmap_info%offsets(idx)%msmsgmgm))
c     &      deallocate(strmap_info%offsets(idx)%msmsgmgm)
c      end do
      deallocate(strmap_info%idx_strmap)
      deallocate(strmap_info%offsets)
      deallocate(strmap_info%maxlen_blk)
c      do idx = 1, strmap_info%mxgraph*strmap_info%mxgraph
c        if (associated(strmap_info%offsets_fc(idx)%ms))
c     &      deallocate(strmap_info%offsets_fc(idx)%ms)
c        if (associated(strmap_info%offsets_fc(idx)%msgm))
c     &      deallocate(strmap_info%offsets_fc(idx)%msgm)
c      end do
      deallocate(strmap_info%idx_fcmap)
      deallocate(strmap_info%offsets_fc)
      deallocate(strmap_info%maxlen_blk_fc)
c      do idx = 1, strmap_info%mxgraph
c        if (associated(strmap_info%offsets_flip(idx)%ms))
c     &      deallocate(strmap_info%offsets_flip(idx)%ms)
c        if (associated(strmap_info%offsets_flip(idx)%msgm))
c     &      deallocate(strmap_info%offsets_flip(idx)%msgm)
c      end do
      deallocate(strmap_info%idx_flipmap)
      deallocate(strmap_info%offsets_flip)
      deallocate(strmap_info%maxlen_blk_flip)
      call file_delete(strmap_info%ffstrmap)

      call mem_clean_vbuffer('bfstrmap')

      return
      end

