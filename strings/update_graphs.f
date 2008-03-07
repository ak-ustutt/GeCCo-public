*----------------------------------------------------------------------*
      subroutine update_graphs(str_info,mel,orb_info)
*----------------------------------------------------------------------*
*
*     successor of unique_graph:
*
*     scan through occupation classes of operator mel%op assiciated
*     with ME-list mel and look for new graphs (i.e. occupation class 
*     for H/P/V alone); set up the list on string info structure str_info
*     
*     the idx_graph array of mel%op is set up accordingly
*
*     andreas, autumn 2006
*              december 2007
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(strinf), intent(inout) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      type(me_list), intent(in) ::
     &     mel

      type(operator), pointer ::
     &     op
      logical ::
     &     unique, same
      integer ::
     &     ngas, nspin,
     &     iocc_cls, ihpv, ica, igas, jgas, idx, idxgr, len(ngastp),
     &     nocc, njoined, ijoin, ioff, igraph
      integer, pointer ::
     &     ihpvgas(:,:)

      logical, external ::
     &     restr_cmp

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'update_graphs')
        write(luout,*) ' unique graphs on input: ',str_info%ngraph
      end if

      nspin = orb_info%nspin
      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      op => mel%op

      njoined = op%njoined

      ! loop over occupation classes of op
      do iocc_cls = 1, op%n_occ_cls
       if(op%formal_blk(iocc_cls))cycle

       ioff = (iocc_cls-1)*njoined
       do ijoin = 1, njoined
        idx = ioff+ijoin

        ! do not forget to init
        mel%idx_graph(1:ngastp,1:2,idx) = 0

        ! loop over h/p/v/x and c/a
        do ihpv = 1,ngastp
          ica_loop: do ica = 1,2

            nocc = op%ihpvca_occ(ihpv,ica,idx)
            if (nocc.eq.0) cycle ica_loop

            unique = .true.
            ! loop over all graphs that are already defined
            igraph_loop: do igraph = 1, str_info%ngraph
c dbg
c            print *,'look at graph # ',igraph,
c     &           str_info%ispc_typ(igraph),
c     &           str_info%ispc_occ(igraph)
c dbg
          
              ! we might compare number of orbitals instead of ihpv
              if (ihpv.eq.str_info%ispc_typ(igraph).and.
     &             nocc.eq.str_info%ispc_occ(igraph)) then

c dbg
c                print *,'comparing: ',ihpv, nocc
c dbg
                ! ... but are restrictions (=shape of graph) the same?
                same = restr_cmp(op%igasca_restr(1,1,1,1,1,idx),
     &                       str_info%igas_restr(1,1,1,1,igraph),
     &                       ica,ihpv,
     &                       ihpvgas,ngas,nspin)
c dbg
c                print *,'same = ',same
c dbg

c                jgas = 1
c                same = .true.
c                cmp_loop: do igas = 1, ngas
c                  if (ihpvgas(igas,1).ne.ihpv) cycle cmp_loop
c                  same = same.and.
c     &                 (op%igasca_restr(1,igas,ica,1,idx).eq.
c     &                  str_info%igas_restr(1,jgas,1,igraph)).and.
c     &                 (op%igasca_restr(2,igas,ica,1,idx).eq.
c     &                  str_info%igas_restr(2,jgas,1,igraph)) .and.
c     &                 (op%igasca_restr(1,igas,ica,2,idx).eq.
c     &                  str_info%igas_restr(1,jgas,2,igraph)).and.
c     &                 (op%igasca_restr(2,igas,ica,2,idx).eq.
c     &                  str_info%igas_restr(2,jgas,2,igraph))
c                  if (.not.same) exit cmp_loop
c                  jgas = jgas+1
c                end do cmp_loop
                
                unique = .not.same

              end if
c dbg
c              print *,'still unique? ',unique
c dbg              

              ! .not.unique == we have this graph already
              ! let op%idx_graph point to this graph ...
              if (.not.unique) mel%idx_graph(ihpv,ica,idx) = igraph
              ! ... and go to next space
              if (.not.unique) cycle ica_loop

            end do igraph_loop
c dbg
c              print *,'final: unique = ',unique
c dbg              

            ! a new unique graph was found:
            call add_graph(ihpv,nocc,ica,
     &           op%igasca_restr(1,1,1,1,1,idx),
     &           str_info,orb_info)

            ! set idxgr to that graph
            idxgr = str_info%ngraph
            mel%idx_graph(ihpv,ica,idx) = idxgr

          end do ica_loop
        end do

       end do ! ijoin
      end do ! iocc_cls

      if (ntest.eq.100) then
        write(luout,*) 'unique graphs on output: ', str_info%ngraph
        len(1:ngastp) = 0
        do igas = 1, ngas
          len(ihpvgas(igas,1)) = len(ihpvgas(igas,1))+1
        end do
        write(luout,*) '    #    typ occ'
        write(luout,*) '                 restrictions'
        do igraph = 1, str_info%ngraph
          write(luout,'(2x,i3,2x,2i4)') igraph,
     &         str_info%ispc_typ(igraph),
     &         str_info%ispc_occ(igraph)
          write(luout,'(15x,10(2i3,x))') 
     &         str_info%igas_restr(1:2,
     &                  1:len(str_info%ispc_typ(igraph)),1,1,igraph)
          if (nspin.eq.2)
     &         write(luout,'(15x,10(2i3,x))') 
     &         str_info%igas_restr(1:2,
     &                  1:len(str_info%ispc_typ(igraph)),1,2,igraph)
        end do
        write(luout,*) 'operator->graph assignments:'
        do iocc_cls = 1, op%n_occ_cls
          if(op%formal_blk(iocc_cls))cycle
          ioff = (iocc_cls-1)*njoined
          write(luout,'(2x,i3,4x,5(4i3,2x))') iocc_cls,
     &         mel%idx_graph(1:ngastp,1,ioff+1:ioff+njoined)
          write(luout,'(5x,4x,5(4i3))') 
     &         mel%idx_graph(1:ngastp,2,ioff+1:ioff+njoined)
        end do
      end if

      return
      end
