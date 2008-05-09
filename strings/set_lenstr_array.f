      subroutine set_lenstr_array(lenstr,nsym,maxidxms,str_info)

      implicit none

      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      integer, intent(in) ::
     &     nsym, maxidxms
      integer, intent(out) ::
     &     lenstr(nsym,maxidxms,*)
      type(strinf), intent(in) ::
     &     str_info

      integer ::
     &     ngraph, igraph, idxms, isym

      integer, pointer ::
     &     lengm(:,:)

      ngraph = str_info%ngraph
      lenstr(1:nsym,1:maxidxms,1:ngraph) = 0
      do igraph = 1, ngraph
        lengm => str_info%g(igraph)%lenstr_gm
        do idxms = 1, maxidxms
          do isym = 1, nsym
            lenstr(isym,idxms,igraph) = lengm(isym,idxms)
          end do
        end do
      end do

      return
      end
