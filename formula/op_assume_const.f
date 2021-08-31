      subroutine op_assume_const(idx_op,val,form,op_info)
*-----------------------------------------------------------------------
*     loop over formula and replace all (scalar) blocks of operator idx_op
*     and modify factor by value val
*-----------------------------------------------------------------------
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     idx_op
      real(8), intent(in) ::
     &     val
      type(formula_item), target, intent(inout) ::
     &     form
      type(operator_info), intent(in) ::
     &     op_info

      type(operator), pointer ::
     &     op_pnt
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next
      integer ::
     &     idx, idx_form_op, idx_form_blk, idx_blk_out, ivtx,jvtx, iarc,
     &     ieqvfac, nvtx, narc, nxarc, nfac, njoined, idx_join, ij, jdx,
     &     isuper, nvtx_rem, nsvtx_rem, narc_rem
      integer, allocatable ::
     &     occ_temp(:,:,:), vtx_chng_idx(:)
      integer, pointer ::
     &     vtx_map(:), vtx_map_rev(:), occ_vtx(:,:,:), vtx_where(:)
      type(contraction), pointer ::
     &     contr
      logical ::
     &     change, all_zero
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     idx_oplist2 !, iocc_zero

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'op_assume_const')
        write(lulog,*) 'idx_op, val: ',idx_op, val
      end if
        
      op_pnt => op_info%op_arr(idx_op)%op
      ! Allocate a temporary array.
      njoined = op_pnt%njoined
      allocate(occ_temp(ngastp,2,njoined),vtx_chng_idx(njoined))

      form_pnt => form

      do
        ! save pointer to next node
        form_pnt_next => form_pnt%next

        ! Navigate to the correct parts of the formula.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.100) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.100)
     &         write(lulog,*) '[INIT TARGET]',form_pnt%target
        case(command_add_contribution)
          ! The necessary contractions are here.
c          write(lulog,*) '[ADD]'

          change = .false.

          contr => form_pnt%contr
          ! Loop over the contraction's vertices.
          nvtx = contr%nvtx
          nvtx_rem = 0
          nsvtx_rem = 0
          
          do idx = 1,nvtx
            ! Check the operator index of the vertex.
            idx_form_op = contr%vertex(idx)%idx_op

            ! If the index of the operator vertex equals that of the
            ! intermediate operator...
            if (idx_form_op.eq.idx_op) then
              idx_form_blk = contr%vertex(idx)%iblk_op

              occ_temp(1:ngastp,1:2,1)=
     &             op_pnt%ihpvca_occ(1:ngastp,1:2,idx_form_blk)
              
              all_zero = iocc_zero(occ_temp(1,1,1))
              
              ! if necessary, look for further vertices of super-vertex
              if (njoined.gt.1) then
                isuper = contr%svertex(idx)
                if (contr%joined(0,isuper).ne.njoined) then
                  write(lulog,*) 'njoined,joined(0):',
     &                 njoined,contr%joined(0,isuper)
                  call quit(1,'op_assume_const','inconsistency')
                end if
                vtx_where => contr%joined(1:njoined,isuper)
                do ij = 2, njoined
                  occ_temp(1:ngastp,1:2,ij) =
     &                 op_pnt%ihpvca_occ(1:ngastp,1:2,
     &                 idx_form_blk-1+ij)

                  all_zero = all_zero.and.iocc_zero(occ_temp(1,1,ij))

                end do
                idx_form_blk = (idx_form_blk-1)/njoined + 1
              end if

              if (all_zero) then
                if(ntest.ge.100)then
                  write(lulog,*) 'Found target in this contraction'
                  call prt_contr2(lulog,contr,op_info)
                endif

                change = .true.
                ! apply the change directly:
                contr%fac = contr%fac * val

                nvtx_rem = nvtx_rem+njoined
                nsvtx_rem = nsvtx_rem+1
                ! remove vertices
                do ij = 1, njoined
                  jvtx = vtx_where(ij)
                  contr%vertex(jvtx)%idx_op = 0 ! mark as deleted
                end do
                ! check for any (zero) arcs
                narc = contr%narc
                do iarc = 1, narc
                  do ij = 1, njoined
                    if (contr%arc(iarc)%link(1).eq.vtx_where(ij).or.
     &                  contr%arc(iarc)%link(2).eq.vtx_where(ij)) then
                      contr%arc(iarc)%link(1:2) = 0
                    end if
                  end do
                end do

                if(ntest.ge.100.and.change)then
                  write(lulog,*) 'After initial removal:'
                  call prt_contr2(lulog,contr,op_info)
                endif
                
              end if
            end if

          end do


c          if (remove) then
c            call delete_fl_node(form_pnt)
c            deallocate(form_pnt)

          if (change) then
! generate a map for new vertex numbers
            allocate(vtx_map(nvtx),vtx_map_rev(nvtx-nvtx_rem))
            jvtx = 0
            vtx_map = 0
            do ivtx = 1, nvtx
              if (contr%vertex(ivtx)%idx_op.ne.0) then
                jvtx = jvtx+1
                vtx_map_rev(jvtx)=ivtx
                vtx_map(ivtx)=jvtx
              end if
            end do

            if (any(vtx_map(1:nvtx-nvtx_rem).eq.0)) then
              write(lulog,*) 'Implemented but untested case!'
              write(luout,*) '-- see source code ---'
              ! check that the update of the vertices and arcs works!
              ! so far, we only had the trivial case of removing the last vertices
              call quit(1,'op_assume_const','check implementation !')
            end if
            
            do ivtx = 1, nvtx-nvtx_rem
              contr%vertex(ivtx) = contr%vertex(vtx_map_rev(ivtx))
            end do
            narc = contr%narc
            narc_rem = 0
            do iarc = 1, narc
              if (contr%arc(iarc)%link(1).ne.0) then
                contr%arc(iarc)%link(1) =
     &               vtx_map(contr%arc(iarc)%link(1))
                contr%arc(iarc)%link(2) =
     &               vtx_map(contr%arc(iarc)%link(2))
              else
                narc_rem = narc_rem + 1
              end if
            end do

            deallocate(vtx_map,vtx_map_rev)

            contr%nvtx = nvtx - nvtx_rem
            contr%nsupvtx = contr%nsupvtx - nsvtx_rem
            
            if(ntest.ge.100)then
              write(lulog,*) 'Updated contraction (0)'
              call prt_contr2(lulog,contr,op_info)
            endif
            
            
            if (narc_rem>0) then
              write(luout,*) 'Sorry, but this is a case never tested.'
              write(luout,*) '-- see source code ---'
              ! so far, it seems that no dummy contractions (with zeros)
              ! ever appeared, so we never could check this route (where we
              ! have to remove arcs)
              call quit(1,'op_assume_const','not fully implemented')
            end if
            
            ! Ensure everything is properly set up.
            narc = contr%narc
            nfac = contr%nfac
            nxarc = contr%nxarc
            call resize_contr(contr,nvtx,narc,nxarc,nfac)

            call update_svtx4contr(contr)
          
            allocate(vtx_map(nvtx)) ! dummy
            call canon_contr(contr,.false.,vtx_map)
            deallocate(vtx_map)
          
            if(ntest.ge.100.and.change)then
              write(lulog,*) 'Updated contraction'
              call prt_contr2(lulog,contr,op_info)
            endif

          end if 

        case default
          ! just ignore unknown commands
        end select

        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
