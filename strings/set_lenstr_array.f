      subroutine set_lenstr_array(lenstr,nsym,maxidxms,str_info)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      integer, parameter :: ntest = 00

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

      if (ntest.ge.100) then
        write(lulog,*) 'the lenstr array:'
        do igraph = 1, ngraph
          write(lulog,*) 'graph ',igraph
          do idxms = 1, maxidxms
            write(lulog,'(1x,i4,2x,8i8)') 
     &           idxms,lenstr(1:nsym,idxms,igraph)
          end do
        end do
      end if

      return
      end
