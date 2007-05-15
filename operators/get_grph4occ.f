*----------------------------------------------------------------------*
      subroutine get_grph4occ(idx_gr,iocc,irst,
     &     str_info,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     get graphs for each HPV/CA from occupation + restriction
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      type(strinf), intent(in) ::
     &     str_info
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas),
     &     iocc(ngastp,2), irst(2,ngas,2,2)

      integer, intent(out) ::
     &     idx_gr(ngastp,2)

      integer ::
     &     ica, igastp, igtyp, ngr4typ, igr4typ,
     &     idxgraph
      logical, external ::
     &     restr_cmp

      do ica = 1, 2
        do igastp = 1, ngastp
          idx_gr(igastp,ica) = 0
          if (iocc(igastp,ica).eq.0) cycle
          ! get graph type
          igtyp = 4*(iocc(igastp,ica)-1) + igastp
          ! number of graphs with same type
          ngr4typ = str_info%gtab(1,igtyp)
          ! check restrictions
          gr4typ: do igr4typ = 1, ngr4typ
            ! actual index of graph
            idxgraph = str_info%gtab(1+igr4typ,igtyp)
            if (restr_cmp(irst,str_info%igas_restr(1,1,1,idxgraph),
     &                       ica,igastp,ihpvgas,ngas)) exit gr4typ
            idxgraph = -idxgraph ! indicate that this was not what
                                 ! we wanted
          end do gr4typ
          if (idxgraph.le.0) then
            write(luout,*) 'ERROR: string not in list'
            write(luout,*) 'Operator was'
            call wrt_occ(luout,iocc)
            call wrt_rstr(luout,irst,ngas)
            write(luout,*) 'C/A, GAS-TYP: ',ica,igastp
            call quit(1,'get_grph4occ','string not in list')
          end if

          idx_gr(igastp,ica) = idxgraph
        end do
      end do
      return
      end
