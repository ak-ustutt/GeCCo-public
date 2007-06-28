*----------------------------------------------------------------------*
      subroutine set_ps_op(oper,name,iocc,irst,mst,igamt,
     &     ngas,ihpvgas,str_info)
*----------------------------------------------------------------------*
*     set up a pseudo-operator (only single block) 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      type(operator), intent(out) ::
     &     oper
      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     ngas,ihpvgas(ngas),
     &     iocc(ngastp,2), irst(2,ngas,2,2), mst, igamt
      type(strinf), intent(in) ::
     &     str_info

      integer ::
     &     ica, ihpv, igtyp, igr4typ, ngr4typ, idxgraph

      integer, external ::
     &     ielsum
      logical, external ::
     &     restr_cmp
      
      oper%dagger = .false.
      oper%gamt = igamt
      oper%mst  = mst

      oper%name = name

      oper%n_occ_cls = 1

      oper%ihpvca_occ(1:ngastp,1:2,1) = iocc(1:ngastp,1:2)
      oper%formal = .false.
      oper%formal_blk(1) = .false.

      oper%ica_occ(1,1) = ielsum(iocc(1,1),ngastp)
      oper%ica_occ(2,1) = ielsum(iocc(1,2),ngastp)

      oper%igasca_restr(1:2,1:ngas,1:2,1:2,1) =
     &             irst(1:2,1:ngas,1:2,1:2)

      do ica = 1, 2
        do ihpv = 1, ngastp

          oper%idx_graph(ihpv,ica,1) = 0
          if (iocc(ihpv,ica).eq.0) cycle
          ! get graph type
          igtyp = 4*(iocc(ihpv,ica)-1) + ihpv
          ! number of graphs with same type
          ngr4typ = str_info%gtab(1,igtyp)
          ! check restrictions
          gr4typ: do igr4typ = 1, ngr4typ
            ! actual index of graph
            idxgraph = str_info%gtab(1+igr4typ,igtyp)
            if (restr_cmp(irst,str_info%igas_restr(1,1,1,idxgraph),
     &                       ica,ihpv,ihpvgas,ngas))
     &           exit gr4typ
            idxgraph = -idxgraph ! indicate that this was not what
                                 ! we wanted
          end do gr4typ
          if (idxgraph.le.0) then
            write(luout,*) 'ERROR: string not in list'
            write(luout,*) 'Operator was'
            call wrt_occ(luout,iocc)
            call wrt_rstr(luout,irst,ngas)
            write(luout,*) 'C/A, GAS-TYP: ',ica,ihpv
            call quit(1,'set_ps_op','string not in list')
          end if

          oper%idx_graph(ihpv,ica,1) = idxgraph
          
        end do
      end do

      return
      end
