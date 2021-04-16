      logical function graph_connected(topo,nn)
      !
      ! analyze graph for connectedness by a depth first search (DFS)
      ! we start at the first vertex and travel along all possible paths
      ! at the end, all vertices must be visited, otherwise the graph is not connected
      ! input: topo -- a topology matrix topo (should be symmetrically filled) 
      !                !=0 - connection  =0 - no connection
      !        nn   -- its dimension
      ! returns: .true. if graph is connected, .false. otherwise
      !        nn = 0 or 1 --> returns .true.
      !
      ! Andreas Koehn, Dec. 2020
      ! 
      implicit none
      
      integer, intent(in) :: nn, topo(nn,nn)

      logical, allocatable :: visited(:)
      
      if (nn.eq.0.or.nn.eq.1) then
         graph_connected = .true.
         return
      end if

      allocate(visited(nn))
      visited(1:nn) = .false.
      
      call dfs(1)
      
      graph_connected = all(visited)
      
      deallocate(visited)

      return

      contains
      
        recursive subroutine dfs(idx)
        
        integer, intent(in) :: idx
        
        integer :: jdx
        
        visited(idx) = .true.
        do jdx = 1, nn
          if (topo(idx,jdx).ne.0 .and. .not.visited(jdx)) then
            call dfs(jdx)
          end if
        end do
        
        end subroutine  
      
      
      end function  
      
