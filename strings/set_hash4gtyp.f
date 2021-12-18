*----------------------------------------------------------------------*
      subroutine set_hash4gtyp(str_info,max_igtyp)
*----------------------------------------------------------------------*
*     set up a little hash-table for finding the matching graph
*     if occupation, HPV-type, and restriction are given
*     the hash value is calculated from occupation and HPV-type
*     gtab(1) contains the number of graphs with this hash-value
*     (i.e. the same occupation and HPV), but different restriction
*     gtab(2:1+gtab(1)) contains the indices of possibly matching graphs 
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      type(strinf), intent(inout) ::
     &     str_info
      integer, intent(in) ::
     &     max_igtyp

      integer ::
     &     igraph, nocc, ihpv, igtyp

      str_info%max_igtyp = max_igtyp
      str_info%gtab(1:ld_gtab,1:max_igtyp) = 0

      ! loop over occupation classes of op
      do igraph = 1, str_info%ngraph

        nocc = str_info%ispc_occ(igraph)
        ihpv = str_info%ispc_typ(igraph)

        ! "hash"-value of that graph
        igtyp = ngastp*(nocc-1) + ihpv
        if (igtyp.gt.max_igtyp) then
          write(lulog,*) 'dimension problem ',igtyp,max_igtyp
          call quit(1,'set_hash4gtyp',
     &       'max_igtyp was not set correctly')
        end if
        ! set up hashing information
        str_info%gtab(1,igtyp) = str_info%gtab(1,igtyp)+1
        if (str_info%gtab(1,igtyp).gt.ld_gtab-1) then
          write(lulog,*) 'dimension problem',str_info%gtab(1,igtyp),
     &             ld_gtab-1
          call quit(1,'set_hash4gtyp',
     &       'parameter ld_gtab (def_strinf.h) is too small')
        end if
        str_info%gtab(1+str_info%gtab(1,igtyp),igtyp) = igraph

      end do

      return
      end
