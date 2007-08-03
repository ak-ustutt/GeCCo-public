      subroutine unique_graph(str_info,max_igtyp,op,ihpvgas,ngas)
*----------------------------------------------------------------------*
*
*     scan through occupation classes of operator op and find the
*     unique graphs (i.e. occupation class for H/P/V alone) and
*     set up the list on string info structure str_info
*     
*     the idx_graph array of op is set accordingly
*
*     andreas, autumn 2006
*
*----------------------------------------------------------------------*

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      integer, parameter ::
     &     ntest = 100

      type(strinf), intent(inout) ::
     &     str_info
      type(operator), intent(in) ::
     &     op
      integer, intent(out) ::
     &     max_igtyp
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas)

      logical ::
     &     unique, same
      integer ::
     &     iocc_cls, ihpv, ica, igas, jgas, idx, idxgr, len(ngastp),
     &     nocc, njoined, ijoin

      if (ntest.ge.100) then
        write(luout,*) '--------------------'
        write(luout,*) 'this is unique_graph'
        write(luout,*) '--------------------'
        write(luout,*) ' unique graphs on input: ',str_info%ngraph
      end if

      max_igtyp = 0

      njoined = op%njoined

      ! loop over occupation classes of op
      do iocc_cls = 1, op%n_occ_cls

       ioff = (iocc_cls-1)*njoined
       do ijoin = 1, njoined
        idx = ioff+ijoin
c dbg
        print *,'idx = ',idx
c dbg

        ! do not forget to init
        op%idx_graph(1:ngastp,1:2,idx) = 0

        ! loop over h/p/v/x and c/a
        do ihpv = 1,ngastp
          ica_loop: do ica = 1,2

            nocc = op%ihpvca_occ(ihpv,ica,idx)
            if (nocc.eq.0) cycle ica_loop

            unique = .true.
            ! loop over all graphs that are already defined
            igraph_loop: do igraph = 1, str_info%ngraph
          
              ! we might compare number of orbitals instead of ihpv
              if (ihpv.eq.str_info%ispc_typ(igraph).and.
     &             nocc.eq.str_info%ispc_occ(igraph)) then

                ! ... but are restrictions (=shape of graph) the same?
                jgas = 1
                same = .true.
                cmp_loop: do igas = 1, ngas
                  if (ihpvgas(igas).ne.ihpv) cycle cmp_loop
                  same = same.and.
     &                 (op%igasca_restr(1,igas,ica,1,idx).eq.
     &                  str_info%igas_restr(1,jgas,1,igraph)).and.
     &                 (op%igasca_restr(2,igas,ica,1,idx).eq.
     &                  str_info%igas_restr(2,jgas,1,igraph)) .and.
     &                 (op%igasca_restr(1,igas,ica,2,idx).eq.
     &                  str_info%igas_restr(1,jgas,2,igraph)).and.
     &                 (op%igasca_restr(2,igas,ica,2,idx).eq.
     &                  str_info%igas_restr(2,jgas,2,igraph))
                  if (.not.same) exit cmp_loop
                  jgas = jgas+1
                end do cmp_loop
                
                unique = .not.same

              end if
              
              ! .not.unique == we have this graph already
              ! let op%idx_graph point to this graph ...
              if (.not.unique) op%idx_graph(ihpv,ica,idx) = igraph
              ! ... and go to next space
              if (.not.unique) cycle ica_loop

            end do igraph_loop

            ! a new unique graph was found:
            str_info%ngraph = str_info%ngraph + 1

            ! set pointer to that graph
            idxgr = str_info%ngraph
            op%idx_graph(ihpv,ica,idx) = idxgr

            str_info%ispc_typ(idxgr) = ihpv
            str_info%ispc_occ(idxgr) = nocc

            ! dimension needed for hash table
            max_igtyp = max(max_igtyp,ngastp*(nocc-1) + ihpv)
            
            jgas = 1
            igas_loop: do igas = 1, ngas
              if (ihpvgas(igas).ne.ihpv) cycle igas_loop
              str_info%igas_restr(1:2,jgas,1:2,idxgr)
     &             = op%igasca_restr(1:2,igas,ica,1:2,idx)
              jgas = jgas+1
            end do igas_loop

          end do ica_loop
        end do

       end do ! ijoin
      end do ! iocc_cls

      if (ntest.eq.100) then
        write(luout,*) 'unique graphs on output: ', str_info%ngraph
        len(1:ngastp) = 0
        do igas = 1, ngas
          len(ihpvgas(igas)) = len(ihpvgas(igas))+1
        end do
        write(luout,*) '    #    typ occ'
        write(luout,*) '                 restrictions'
        do igraph = 1, str_info%ngraph
          write(luout,'(2x,i3,2x,2i4)') igraph,
     &         str_info%ispc_typ(igraph),
     &         str_info%ispc_occ(igraph)
          write(luout,'(15x,10(2i3,x))') 
     &         str_info%igas_restr(1:2,
     &                  1:len(str_info%ispc_typ(igraph)),1,igraph)
        end do
        write(luout,*) 'operator->graph assignments:'
        do iocc_cls = 1, op%n_occ_cls
          ioff = (iocc_cls-1)*njoined
          write(luout,'(2x,i3,4x,5(4i3,2x))') iocc_cls,
     &         op%idx_graph(1:ngastp,1,ioff+1:ioff+njoined)
          write(luout,'(5x,4x,5(4i3))') 
     &         op%idx_graph(1:ngastp,2,ioff+1:ioff+njoined)
        end do
      end if

      return
      end
