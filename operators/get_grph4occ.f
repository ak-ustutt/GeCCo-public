*----------------------------------------------------------------------*
      subroutine get_grph4occ(idx_gr,iocc,irst,njoined,
     &     str_info,orb_info,allow_new_g)
! allow_new_g -> error_handling
*----------------------------------------------------------------------*
*     get graphs for each HPV/CA from occupation + restriction
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined),
     &     irst(2,orb_info%ngas,2,2,orb_info%nspin,njoined)
      logical, intent(in) ::
     &     allow_new_g

      integer, intent(out) ::
     &     idx_gr(ngastp,2,njoined)

      integer ::
     &     nspin, ngas, ica, igastp, igtyp, ngr4typ, igr4typ,
     &     idxgraph, ijoin, ii
      integer, pointer ::
     &     ihpvgas(:,:) 
      logical, external ::
     &     restr_cmp

      ngas = orb_info%ngas
      nspin = orb_info%nspin
      ihpvgas => orb_info%ihpvgas
c dbg
c      print *,'get_grph4occ: '
c      call wrt_rstr(6,irst,ngas)
c dbg

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
c dbg
c          print *,'ngr4typ = ',ngr4typ
c dbg
          gr4typ: do igr4typ = 1, ngr4typ
            ! actual index of graph
            idxgraph = str_info%gtab(1+igr4typ,igtyp)
c dbg
c            print *,'comparing to ',idxgraph
c            print '(10(i2,x,i2,2x))',
c     &           str_info%igas_restr(1:2,1:ngas,1,1,idxgraph)
c dbg
            if (restr_cmp(irst(1,1,1,1,1,ijoin),
     &                    str_info%igas_restr(1,1,1,1,idxgraph),
     &                    ica,igastp,ihpvgas,ngas,nspin)) exit gr4typ
            idxgraph = -idxgraph ! indicate that this was not what
                                 ! we wanted
          end do gr4typ
          if (idxgraph.le.0) then
            if (.not.allow_new_g) then
              idx_gr(1:ngastp,1:2,1:njoined) = -1
              exit outer_loop
            end if
           ! add required graph (needed e.g. in case of CASSCF)
           call add_graph(iocc(igastp,ica,ijoin),igastp,ica,
     &              irst(1,1,1,1,1,ijoin),str_info,orb_info)
           idxgraph = str_info%ngraph

          end if
            
          idx_gr(igastp,ica,ijoin) = idxgraph
        end do
       end do
      end do outer_loop
c dbg
c      print *,'result: '
c      call wrt_occ_n(6,idx_gr,njoined)
c dbg

      return
      end
