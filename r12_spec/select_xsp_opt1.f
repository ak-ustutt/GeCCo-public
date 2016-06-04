      subroutine select_xsp_opt1(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     some quick fix modifications for XSPopt with optimized singles only
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest =  00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
!      include 'def_orbinf.h'
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
     &     error
      integer ::
     &     idx_tbar, idx_top, idx_cbar, idx_cex, idx_cbarf, idx_cexf
      integer ::
     &     ii, nvtx, ivtx, idx_op, iblk_op,iblk_new, iblk_cf1, iblk_cbf1
      integer ::
     &     idxop(nlabels), occ(ngastp,2)


      integer, pointer ::
     &     occ_cbarf(:,:,:), occ_cexf(:,:,:)
      type(operator), pointer ::
     &     op_top, op_tbar, op_cex, op_cbar, op_cexf, op_cbarf
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2, iblk_occ

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_xsp_opt1')
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.6
      
      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.4)
     &       call quit(1,'select_xsp_opt1','need exactly 4 labels')
        call quit(1,'select_xsp_opt1','Labels not on list!')
      end if

      ! mode is unused currently

      ! operator indices
      idx_tbar  = idxop(1)
      idx_cbarf = idxop(2)
      idx_cbar  = idxop(3)
      idx_top   = idxop(4)
      idx_cexf  = idxop(5)
      idx_cex   = idxop(6)

      ! some useful pointers
      op_tbar => op_info%op_arr(idx_tbar)%op
      op_cbarf=> op_info%op_arr(idx_cbarf)%op
      op_cbar => op_info%op_arr(idx_cbar)%op
      op_top  => op_info%op_arr(idx_top)%op
      op_cexf => op_info%op_arr(idx_cexf)%op
      op_cex  => op_info%op_arr(idx_cex)%op

      occ_cbarf => op_cbarf%ihpvca_occ
      occ_cexf  => op_cexf%ihpvca_occ

      ! obtain blocks of CBAR1 and CEX1
      occ = 0
      occ(IHOLE,1) = 1
      occ(IPART,2) = 1
      iblk_cbf1 = iblk_occ(occ,.false.,op_cbarf,1)
      occ = 0
      occ(IPART,1) = 1
      occ(IHOLE,2) = 1
      iblk_cf1  = iblk_occ(occ,.false.,op_cexf,1)

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          ! replace all non_single cexbar/cex by tbar/top
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idx_cexf) then
              if (iblk_op.ne.iblk_cf1 ) then ! replace by T2, T3, ...
                
                iblk_new =  ! I set it hard-coded to blk_version = 1
     &               iblk_occ(occ_cexf(1:,1:,iblk_op),.false.,op_top,1)

                vertex(ivtx)%idx_op = idx_top
                vertex(ivtx)%iblk_op = iblk_new

              else ! replace by T1tilde (nomenclature of the papers)
                   ! here it is op_cex              
                iblk_new =  ! I set it hard-coded to blk_version = 1
     &               iblk_occ(occ_cexf(1:,1:,iblk_op),.false.,op_cex,1)

                vertex(ivtx)%idx_op = idx_cex
                vertex(ivtx)%iblk_op = iblk_new

              end if
            else if (idx_op.eq.idx_cbarf) then
              if (iblk_op.ne.iblk_cbf1 ) then

                iblk_new =  ! I set it hard-coded to blk_version = 1
     &              iblk_occ(occ_cbarf(1:,1:,iblk_op),.false.,op_tbar,1)

                vertex(ivtx)%idx_op = idx_tbar
                vertex(ivtx)%iblk_op = iblk_new

              else

                iblk_new =  ! I set it hard-coded to blk_version = 1
     &            iblk_occ(occ_cbarf(1:,1:,iblk_op),.false.,op_cbar,1)

                vertex(ivtx)%idx_op = idx_cbar
                vertex(ivtx)%iblk_op = iblk_new

              end if
            end if
          end do

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_xsp_opt1','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
