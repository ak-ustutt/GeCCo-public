      integer, parameter ::
     &     ld_gtab = 10
      type strinf
        integer ::
     &     ngraph               ! number of graphs
        integer, pointer ::
     &     gtab(:,:)      ! little "hash"-table for finding
     &                    ! graphs if ispc_occ, and ispc_typ are given
        integer, pointer ::
     &       ispc_typ(:),       ! hole/particle/valence space
     &       ispc_occ(:),       ! occupation of space
     &       igas_restr(:,:,:,:)! associated restrictions per GAS
        type(graph), pointer ::
     &       g(:)               ! info on graph (see def_graph)
      end type strinf
