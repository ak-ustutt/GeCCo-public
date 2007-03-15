      type graph
        integer ::
     &       ndis,              ! number of subspace distributions (total)
     &       ndis_a             ! number of allowed distributions
        integer, pointer ::
     &       y4sg(:),           ! arc weights of graphs in subspace
     &       yinf(:),           ! info array (offsets for graphs)
     &       yssg(:),           ! arc weights for subspace graphs
     &       wssg(:),           ! vertex weights for subspace graphs
     &       idis_m(:),         ! masked distributions (indicated by 0)
     &       lenstr_dgm(:),     ! number of strings per distribution, 
     &                          ! IRREP and MS
     &       ioffstr_dgm(:),    ! offsets per distribution within IRREP, MS
     &       lenstr_gm(:,:)     ! number of strings per IRREP and MS
      end type graph
