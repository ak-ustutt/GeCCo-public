*----------------------------------------------------------------------*
      subroutine check_graph(igraph,str_info,orb_info)
*----------------------------------------------------------------------*
*     check graph addressing and string sequence generation routine
*     for consistency
*
*     andreas, august 2008
*      
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, intent(in) ::
     &     igraph
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ihpv

      ihpv = str_info%ispc_typ(igraph)
      call check_graph_kernel(igraph,str_info%ispc_occ(igraph),
     &     str_info%igas_restr(1,1,1,1,igraph),
     &     str_info%g(igraph)%y4sg,
     &     str_info%g(igraph)%yinf,
     &     str_info%g(igraph)%yssg,
     &     str_info%g(igraph)%wssg,
     &     str_info%g(igraph)%ioffstr_dgm,
     &     str_info%g(igraph)%ndis,
     &     orb_info%mostnd(1,1,orb_info%idx_gas(ihpv)),
     &     orb_info%nsym,orb_info%ngas_hpv(ihpv),orb_info%igamorb
     &     )

      return
      end
*----------------------------------------------------------------------*
      subroutine check_graph_kernel(
     &     igraph,iocc,irestr,
     &     g_y4sg,g_yinf,
     &     g_yssg,g_wssg,
     &     g_off_dgm,g_ndis,
     &     mostnd_cur,nsym,ngas_cur,igamorb)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     igraph, iocc, irestr(2,ngas_cur,2),
     &     nsym, ngas_cur,
     &     mostnd_cur(2,nsym),igamorb(*),
     &     g_y4sg(*),g_yinf(*),
     &     g_yssg(*),g_wssg(*),
     &     g_off_dgm(*), g_ndis

      logical ::
     &     first, error
      integer ::
     &     idorb(iocc),idspn(iocc),idgam(iocc),idss(iocc),
     &     idxms, igam,
     &     istr, nstr, idx, ms, idxmap

      integer, external ::
     &     idx4sg, std_spsign
      logical, external ::
     &     next_string

      error = .false.
      do ms = iocc, -iocc, -2
        do igam = 1, nsym

          idxmap = 0

          first = .true.
          istr = 0
          ! loop over strings
          do while(next_string(idorb,idspn,idss,
     &         iocc,ms,igam,first,
     &         irestr,
     &         mostnd_cur,igamorb,
     &         nsym,ngas_cur))

            istr = istr+1
            first = .false.
        
            do idx = 1, iocc
              idgam(idx) = igamorb(idorb(idx))
            end do

            ! get index from graphs
            idx = idx4sg(iocc,idss,idorb,idspn,idgam,
     &           g_y4sg,g_yinf,
     &           g_yssg,g_wssg,g_off_dgm,g_ndis,mostnd_cur,
     &           iocc,nsym,ngas_cur)+1

            ! consistent??
            if (idx.ne.istr) then
              error = .true.
              write(luout,*) 'ERROR for GRAPH/MS/GAMMA: ',igraph,ms,igam
              write(luout,*) 'idorb:   ',idorb(1:iocc)
              write(luout,*) 'idspn:   ',idspn(1:iocc)
              write(luout,*) 'idgam:   ',idgam(1:iocc)
              write(luout,*) 'idss:    ',idss (1:iocc)
              write(luout,*) 'istr, idx = ',istr, idx
            end if

          end do

        end do
      end do

      if (error) then
        call quit(1,'check_graph','inconsistency')
      end if

      return
      end
