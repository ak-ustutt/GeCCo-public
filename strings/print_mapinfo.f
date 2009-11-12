      subroutine print_mapinfo(luout,map_info,n12)
      ! 
      ! print the content of the mysterious mapping information
      !
      implicit none
      
      integer, intent(in) ::
     &     luout, n12, map_info(*)

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
        write(luout,*) 'idx12, nsplit1: ',idx12,nsplit
        write(luout,*) 'indices: ',
     &         map_info(idx_minf+1:idx_minf+nsplit)
        idx_minf = idx_minf+nsplit

        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        write(luout,*) 'idx12, nsplit2: ',idx12,nsplit
        write(luout,*) 'indices: ',
     &         map_info(idx_minf+1:idx_minf+nsplit)
        idx_minf = idx_minf+nsplit
      end do

      end
