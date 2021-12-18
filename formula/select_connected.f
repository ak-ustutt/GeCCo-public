      subroutine select_connected(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     go through terms in flist and make sure that the operators on
*     array labels are fully connected.
*
*     Contractions that are not connected in this sense are deleted
*     A warning is issued, if the operators on the list do not at all 
*     appear in a contraction
*
*     mode is currently unused
*
*     Andreas Koehn, Dec. 2020
*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest =  000
      character(len=16), parameter ::
     &     i_am = 'select_connected'

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     nlabels
      character(len=*), intent(in) ::
     &     labels(nlabels), mode

      logical ::
     &     delete, error, connected
      integer ::
     &     nop, ii, nvtx, narc, nj, 
     &     ij, ivtx, iarc, idxsvtx, ilink, jlink, iop, jop
      integer ::
     &     idxop(nlabels)
     
      integer, allocatable ::
     &     vtxlist(:), vtxmap(:), topo(:,:)

      integer, pointer ::
     &     svertex(:), joined(:,:)

      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     graph_connected 

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,i_am)
        write(lulog,*) 'mode = ',trim(mode)
        do ii = 1, nlabels
          write(lulog,*) ' label: ',trim(labels(ii))
        end do
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.le.0
      
      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.le.0)
     &       call quit(1,i_am,'no operator label given?')
        call quit(1,i_am,'Labels not on list!')
      end if

      form_pnt => flist
      do ! loop over formula entries
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case (form_pnt%command)
        case (command_end_of_formula)
          if (ntest.ge.1000) write(lulog,*) '[END]'
        case (command_set_target_init)
          if (ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case (command_add_contribution)

          contr => form_pnt%contr

          nvtx = contr%nvtx
          narc = contr%narc
          vertex  => contr%vertex
          arc     => contr%arc
          svertex => contr%svertex
          joined  => contr%joined

          if (contr%nxarc.gt.0) call quit(1,i_am,
     &         'not yet adapted for open diagrams')
 
          allocate(vtxlist(nvtx),vtxmap(nvtx))
          vtxmap(1:nvtx) = 0

          delete = .false.
          if(ntest.ge.1000)then
            write(lulog,*) 'Current formula item:'
            call prt_contr2(lulog,contr,op_info)
          end if
  
          ! count relevant operators in this contraction
          nop = 0
          idxsvtx = 0
          do ivtx = 1, nvtx
            do ii = 1, nlabels
              ! check if operator on list and if this is our first visit to this supervertex
              ! new supervertex numbers always come in ascending order
              if (vertex(ivtx)%idx_op.eq.idxop(ii)
     &           .and.svertex(ivtx).gt.idxsvtx) then
                idxsvtx = svertex(ivtx)  ! get corresponding supervertex
                nop = nop+1
                vtxlist(nop) = ivtx     ! collect this vertex on our vertex list
                vtxmap(ivtx) = nop       ! also create a mapping from the vertex list to the operator list
                nj = joined(0,idxsvtx)   ! for the mapping, we also need the additional vertices (in case this is a multivertex op.)
                do ij = 2, nj
                  vtxmap(joined(ij,idxsvtx)) = nop
                end do
                exit
              end if
            end do
          end do
          
          if (ntest.ge.1000) then
            write(lulog,'(1x,"nop = ",i4)') nop
            write(lulog,'(1x,"vtxlist:  ",10i4)') vtxlist(1:nop)
            write(lulog,'(1x,"vtxmap:   ",10i4)') vtxmap(1:nvtx)
          end if
 
          if (nop.eq.0) then
            call warn(i_am,'no target operator in contraction')
          else
   
            ! allocate topo
            allocate(topo(nop,nop))

            ! set up topo
            topo = 0
            do iop = 1, nop ! loop over relevant operators
              ivtx = vtxlist(iop)      ! get vertex number
              idxsvtx = svertex(ivtx)  ! and super vertex number
              nj = joined(0,idxsvtx)   ! check if operator has further joined vertices
              do ij = 1, nj            ! loop over all joined vertices (first one is equal to ivtx)
                ivtx = joined(ij,idxsvtx)
                do iarc = 1, narc        ! loop over connection arcs
                  do ilink = 1, 2        ! and the two entries on the link
                    jlink = 3-ilink
                    if (arc(iarc)%link(ilink).eq.ivtx) then  ! current link on our list?
                      jop = vtxmap(arc(iarc)%link(jlink))     ! check via map that other link is also on list
                      if (jop.gt.0) then
                        topo(iop,jop) = 1                     ! set topo map (symmetrically)
                        topo(jop,iop) = 1
                      end if
                    end if
                  end do
                end do
              end do
            end do
            
            if (ntest.ge.1000) then
              write(lulog,'(1x,"topo map")')
              do iop = 1, nop
                write(lulog,'(1x,20i4)') topo(iop,1:nop)
              end do
            end if
  
              ! analyze topo by a graph search:
            delete = .not.graph_connected(topo,nop)
        
            ! deallocate topo etc.
            deallocate(topo,vtxmap,vtxlist)
  
            if (delete) then
              ! Print the deleted contraction.
              if(ntest.ge.1000)then
                write(lulog,*) 'Deleted formula item:'
                call prt_contr2(lulog,form_pnt%contr,op_info)
              end if
  
              ! Delete the node.
              call delete_fl_node(form_pnt)
              deallocate(form_pnt)
            end if

          end if    

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,i_am,'command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      end do

      return
      end
      
      
