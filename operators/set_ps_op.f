*----------------------------------------------------------------------*
      subroutine set_ps_op(oper,name,iocc,irst,njoined,mst,igamt,
     &     orb_info,str_info)
*----------------------------------------------------------------------*
*     set up a pseudo-operator (only single block) 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_orbinf.h'
      include 'def_strinf.h'

      type(operator), intent(out) ::
     &     oper
      character, intent(in) ::
     &     name*(*)
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined), irst(2,orb_info%ngas,2,2,njoined),
     &     mst, igamt
      type(strinf), intent(in) ::
     &     str_info

      integer ::
     &     ica, ihpv, igtyp, igr4typ, ngr4typ, idxgraph, idx, ijoin

      integer, pointer ::
     &     ngas, ihpvgas(:)

      integer, external ::
     &     ielsum
      logical, external ::
     &     restr_cmp
      
      ngas => orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      oper%dagger = .false.
      oper%gamt = igamt
      oper%mst  = mst

      oper%name = name

      oper%n_occ_cls = 1
      oper%njoined = njoined

      ! initial allocation
      call init_operator(0,oper,orb_info)

      oper%ihpvca_occ(1:ngastp,1:2,1:njoined) =
     &           iocc(1:ngastp,1:2,1:njoined)

      oper%ica_occ(1:2,1) = 0
      do ijoin = 1, njoined
        oper%ica_occ(1,1) = oper%ica_occ(1,1)+
     &       ielsum(iocc(1,1,ijoin),ngastp)
        oper%ica_occ(2,1) = oper%ica_occ(2,1)+
     &       ielsum(iocc(1,2,ijoin),ngastp)
      end do

      oper%igasca_restr(1:2,1:ngas,1:2,1:2,1:njoined) =
     &             irst(1:2,1:ngas,1:2,1:2,1:njoined)

      do ijoin = 1, njoined
        do ica = 1, 2
          do ihpv = 1, ngastp

            oper%idx_graph(ihpv,ica,ijoin) = 0
            if (iocc(ihpv,ica,ijoin).eq.0) cycle
            ! get graph type
            igtyp = 4*(iocc(ihpv,ica,ijoin)-1) + ihpv
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
              call wrt_occ_n(luout,iocc,njoined)
              do idx = 1, njoined
                call wrt_rstr(luout,irst(1,1,1,1,ijoin),ngas)
              end do
              write(luout,*) 'C/A, GAS-TYP, vtx: ',ica,ihpv,ijoin
              call quit(1,'set_ps_op','string not in list')
            end if

            oper%idx_graph(ihpv,ica,ijoin) = idxgraph
          
          end do
        end do
      end do

      return
      end
