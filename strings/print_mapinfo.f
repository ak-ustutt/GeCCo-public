      subroutine print_mapinfo(lulog,map_info,n12)
      ! 
      ! print the content of the mysterious mapping information
      !
      implicit none
      
      integer, intent(in) ::
     &     lulog, n12, map_info(*)

      integer ::
     &     idx_minf, idx12, idx1, idx2, nsplit

      idx_minf = 0
      idx12 = 0
      do while(idx12.lt.n12)
        idx12 = idx12+1
        idx1 = 0
        idx2 = 0
        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        write(lulog,'("vtx",i4,": vtxs from op1: ",10i4)')
     &       idx12,map_info(idx_minf+1:idx_minf+nsplit)
     &         
        idx_minf = idx_minf+nsplit

        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        write(lulog,'(3x,4x,2x,"vtxs from op2: ",10i4)')
     &       idx12,map_info(idx_minf+1:idx_minf+nsplit)
        idx_minf = idx_minf+nsplit
      end do

      end
