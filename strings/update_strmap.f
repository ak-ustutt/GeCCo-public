*----------------------------------------------------------------------*
      subroutine update_strmap(str_info,strmap_info)
*----------------------------------------------------------------------*
*     update string mapping info (in case that we have additional
*     graphs)
*----------------------------------------------------------------------*

      implicit none
      
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'

      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info

      integer ::
     &     ngraph, ngraph_old, idx, igraph, jgraph, idx_old, idx_new

      integer, pointer ::
     &     new_idx_strmap(:), new_maxlen_blk(:)
      type(strmap_offsets), pointer ::
     &     new_offsets(:)

      ! everything OK?
      if (str_info%ngraph.le.strmap_info%mxgraph) return

c dbg
      print *,'!!UPDATING STRMAP_INFO!!'
c dbg

      ngraph = str_info%ngraph
      ngraph_old = strmap_info%mxgraph
      strmap_info%mxgraph = ngraph
      allocate(new_idx_strmap(ngraph*ngraph))

      new_idx_strmap(1:ngraph*ngraph) = -1

      allocate(new_offsets(ngraph*ngraph))
      allocate(new_maxlen_blk(ngraph*ngraph))
      do idx = 1, ngraph*ngraph
        nullify(new_offsets(idx)%msms)
        nullify(new_offsets(idx)%msmsgmgm)
      end do

      do igraph = 1, ngraph_old
        idx_old = (igraph-1)*ngraph_old
        idx_new = (igraph-1)*ngraph
        do jgraph = 1, ngraph_old
          idx_new = idx_new+1
          idx_old = idx_old+1
          new_idx_strmap(idx_new) = strmap_info%idx_strmap(idx_old)
          new_maxlen_blk(idx_new) = strmap_info%maxlen_blk(idx_old)
          new_offsets(idx_new)%msms =>
     &         strmap_info%offsets(idx_old)%msms
          new_offsets(idx_new)%msmsgmgm =>
     &         strmap_info%offsets(idx_old)%msmsgmgm
        end do
      end do

      deallocate(strmap_info%idx_strmap,strmap_info%offsets)

      strmap_info%idx_strmap => new_idx_strmap
      strmap_info%offsets    => new_offsets

      return
      end

