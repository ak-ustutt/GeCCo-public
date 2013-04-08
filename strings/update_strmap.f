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
     &     new_idx_strmap(:), new_maxlen_blk(:),
     &     new_idx_flipmap(:), new_maxlen_blk_flip(:),
     &     new_idx_fcmap(:), new_maxlen_blk_fc(:),
     &     new_idx_spprjmap(:), new_maxlen_blk_spprj(:)
      type(strmap_offsets), pointer ::
     &     new_offsets(:)
      type(flpmap_offsets), pointer ::
     &     new_offsets_flip(:),
     &     new_offsets_spprj(:),
     &     new_offsets_fc(:)

      ! everything OK?
      if (str_info%ngraph.le.strmap_info%mxgraph) return

c dbg
c      print *,'!!UPDATING STRMAP_INFO!!'
c dbg

      ngraph = str_info%ngraph
      ngraph_old = strmap_info%mxgraph
      strmap_info%mxgraph = ngraph
      allocate(new_idx_strmap(ngraph*ngraph))
      allocate(new_idx_fcmap(ngraph*ngraph))
      allocate(new_idx_flipmap(ngraph))
      allocate(new_idx_spprjmap(ngraph))

      new_idx_strmap(1:ngraph*ngraph) = -1
      new_idx_fcmap (1:ngraph*ngraph) = -1
      new_idx_flipmap(1:ngraph) = -1
      new_idx_spprjmap(1:ngraph) = -1

      allocate(new_offsets(ngraph*ngraph))
      allocate(new_maxlen_blk(ngraph*ngraph))
      do idx = 1, ngraph*ngraph
        nullify(new_offsets(idx)%msms)
        nullify(new_offsets(idx)%msmsgmgm)
      end do

      allocate(new_offsets_fc(ngraph*ngraph))
      allocate(new_maxlen_blk_fc(ngraph*ngraph))
      do idx = 1, ngraph*ngraph
        nullify(new_offsets_fc(idx)%ms)
        nullify(new_offsets_fc(idx)%msgm)
      end do

      allocate(new_offsets_flip(ngraph))
      allocate(new_maxlen_blk_flip(ngraph))
      do idx = 1, ngraph
        nullify(new_offsets_flip(idx)%ms)
        nullify(new_offsets_flip(idx)%msgm)
      end do

      allocate(new_offsets_spprj(ngraph))
      allocate(new_maxlen_blk_spprj(ngraph))
      do idx = 1, ngraph
        nullify(new_offsets_spprj(idx)%ms)
        nullify(new_offsets_spprj(idx)%msgm)
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

      do igraph = 1, ngraph_old
        idx_old = (igraph-1)*ngraph_old
        idx_new = (igraph-1)*ngraph
        do jgraph = 1, ngraph_old
          idx_new = idx_new+1
          idx_old = idx_old+1
          new_idx_fcmap(idx_new) = strmap_info%idx_fcmap(idx_old)
          new_maxlen_blk_fc(idx_new) =
     &                            strmap_info%maxlen_blk_fc(idx_old)
          new_offsets_fc(idx_new)%ms =>
     &         strmap_info%offsets_fc(idx_old)%ms
          new_offsets_fc(idx_new)%msgm =>
     &         strmap_info%offsets_fc(idx_old)%msgm
        end do
      end do

      do igraph = 1, ngraph_old
        new_idx_flipmap(igraph) = strmap_info%idx_flipmap(igraph)
        new_maxlen_blk_flip(igraph) =
     &                            strmap_info%maxlen_blk_flip(igraph)
        new_offsets_flip(igraph)%ms =>
     &       strmap_info%offsets_flip(igraph)%ms
        new_offsets_flip(igraph)%msgm =>
     &       strmap_info%offsets_flip(igraph)%msgm
      end do

      do igraph = 1, ngraph_old
        new_idx_spprjmap(igraph) = strmap_info%idx_spprjmap(igraph)
        new_maxlen_blk_spprj(igraph) =
     &                            strmap_info%maxlen_blk_spprj(igraph)
        new_offsets_spprj(igraph)%ms =>
     &       strmap_info%offsets_spprj(igraph)%ms
        new_offsets_spprj(igraph)%msgm =>
     &       strmap_info%offsets_spprj(igraph)%msgm
      end do

      deallocate(strmap_info%idx_strmap,strmap_info%maxlen_blk,
     &           strmap_info%offsets)

      strmap_info%idx_strmap => new_idx_strmap
      strmap_info%maxlen_blk => new_maxlen_blk
      strmap_info%offsets    => new_offsets

      deallocate(strmap_info%idx_fcmap,strmap_info%maxlen_blk_fc,
     &           strmap_info%offsets_fc)

      strmap_info%idx_fcmap  => new_idx_fcmap
      strmap_info%maxlen_blk_fc => new_maxlen_blk_fc
      strmap_info%offsets_fc => new_offsets_fc

      deallocate(strmap_info%idx_flipmap,strmap_info%maxlen_blk_flip,
     &           strmap_info%offsets_flip)

      strmap_info%idx_flipmap  => new_idx_flipmap
      strmap_info%maxlen_blk_flip => new_maxlen_blk_flip
      strmap_info%offsets_flip => new_offsets_flip

      deallocate(strmap_info%idx_spprjmap,strmap_info%maxlen_blk_spprj,
     &           strmap_info%offsets_spprj)

      strmap_info%idx_spprjmap  => new_idx_spprjmap
      strmap_info%maxlen_blk_spprj => new_maxlen_blk_spprj
      strmap_info%offsets_spprj => new_offsets_spprj

      return
      end

