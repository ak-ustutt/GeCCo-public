*----------------------------------------------------------------------*
      subroutine init_strmap(str_info,strmap_info)
*----------------------------------------------------------------------*
*     allocate string mapping info and associated file
*----------------------------------------------------------------------*

      implicit none
      
      include 'ioparam.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     initial_ngraph = 10

      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info

      integer ::
     &     ngraph, ifree, idx

      ngraph = max(initial_ngraph,str_info%ngraph)
      strmap_info%mxgraph = ngraph
      allocate(strmap_info%idx_strmap(ngraph*ngraph))

      strmap_info%idx_strmap(1:ngraph*ngraph) = -1
      strmap_info%idx_last = 0
      allocate(strmap_info%offsets(ngraph*ngraph))
      allocate(strmap_info%maxlen_blk(ngraph*ngraph))
      do idx = 1, ngraph*ngraph
        nullify(strmap_info%offsets(idx)%msms)
        nullify(strmap_info%offsets(idx)%msmsgmgm)
      end do

      call file_init(strmap_info%ffstrmap,
     &     'strmaps.da',ftyp_da_unf,lblk_da)

      call file_open(strmap_info%ffstrmap)

      ! currently: hard coded dimensioning
      call mem_init_vbuffer(strmap_info%ffstrmap,
     &     'bfstrmap',256*lblk_da,1024)

      return
      end

