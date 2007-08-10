*----------------------------------------------------------------------*
      subroutine get_grph4occ(idx_gr,iocc,irst,
     &     str_info,ihpvgas,ngas,njoined,error_exit)
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
     &     ngas, ihpvgas(ngas), njoined,
     &     iocc(ngastp,2,njoined), irst(2,ngas,2,2,njoined)
      logical, intent(in) ::
     &     error_exit

      integer, intent(out) ::
     &     idx_gr(ngastp,2,njoined)

      integer ::
     &     ica, igastp, igtyp, ngr4typ, igr4typ,
     &     idxgraph, ijoin, ii
      logical, external ::
     &     restr_cmp

      outer_loop: do ijoin = 1, njoined
       do ica = 1, 2
        do igastp = 1, ngastp
          idx_gr(igastp,ica,ijoin) = 0
          if (iocc(igastp,ica,ijoin).eq.0) cycle
          ! get graph type
          igtyp = 4*(iocc(igastp,ica,ijoin)-1) + igastp
          ! number of graphs with same type
          if (igtyp.le.str_info%max_igtyp) then
            ngr4typ = str_info%gtab(1,igtyp)
          else
            ngr4typ = 0
          end if
          idxgraph = -1
          ! check restrictions
          gr4typ: do igr4typ = 1, ngr4typ
            ! actual index of graph
            idxgraph = str_info%gtab(1+igr4typ,igtyp)
            if (restr_cmp(irst(1,1,1,1,ijoin),
     &                    str_info%igas_restr(1,1,1,idxgraph),
     &                    ica,igastp,ihpvgas,ngas)) exit gr4typ
            idxgraph = -idxgraph ! indicate that this was not what
                                 ! we wanted
          end do gr4typ
          if (idxgraph.le.0) then
            if (error_exit) then
              write(luout,*) 'ERROR: string not in list'
              write(luout,*) 'Operator was'
              call wrt_occ_n(luout,iocc,njoined)
              do ii = 1, njoined
                call wrt_rstr(luout,irst(1,1,1,1,ii),ngas)
              end do
              write(luout,*) 'C/A, GAS-TYP, VTX: ',ica,igastp,ijoin
              call quit(1,'get_grph4occ','string not in list')
            else
              idx_gr(1:ngastp,1:2,ijoin) = -1
              exit outer_loop
            end if
          end if
            
          idx_gr(igastp,ica,ijoin) = idxgraph
        end do
       end do
      end do outer_loop

      return
      end
