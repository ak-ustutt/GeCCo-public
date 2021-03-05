*---------------------------------------------------------------------*
      subroutine itf_set_code_targets(itf_targets,target_list,nlist,
     &                                op_info)
*---------------------------------------------------------------------*
*
*     Parse the list "target_list"; Every entry in <> is interpreted
*     as a ITF code name, all other entries should be defined GeCCo
*     operator names that are listed as "targets" in the submitted
*     formula. E.g.
*     '<Compute_INTpppp>','INTpppp','<Residual>','O1','O2','EMRCC'
*     generates 2 code sections, one in which the target INTpppp is
*     computed and one in which the residuals O1, O2, and the energy
*     are computed. The sequence of the targets does not need to be
*     the same as on the formula list, but targets for one "code"
*     should be contiguous on the formula list
*
*     andreas, march '21
*---------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'
      include 'def_target.h'

      character(len=20), parameter ::
     &     i_am = 'itf_set_code_targets'
      integer, parameter ::
     &     ntest = 100
      
      integer, intent(in) ::
     &     nlist
      character(len=len_command_par), intent(in) ::
     &     target_list(nlist)
      type(code_targets), intent(inout) ::
     &     itf_targets
      type(operator_info) ::
     &     op_info

      logical ::
     &     is_code
      integer ::
     &     ii, ncodes, ntargets, idxcode, idxtarget, litem, idxop
      character(len=len_command_par) ::
     &     item

      integer, external ::
     &     idx_oplist2
      
      ! check the list for consistency and count
      ncodes = 0
      ntargets = 0
      do ii = 1, nlist
        item = trim(adjustl(target_list(ii)))
        litem = len_trim(item)
        is_code = item(1:1)=='<'.and.item(litem:litem)=='>'
        if (is_code) then
          ncodes = ncodes+1
        else
          if (ii==1) then
            write(lulog,*) 'Found:', item
            call quit(0,i_am,'First entry must be <code name>')
          end if
          ntargets = ntargets+1
        end if
      end do

      itf_targets%ncodes = ncodes
      itf_targets%ntargets = ntargets
      allocate(itf_targets%code_name(ncodes))
      allocate(itf_targets%idx_target(ntargets))
      allocate(itf_targets%idx_code(ntargets))

      itf_targets%idx_target = 0
      itf_targets%idx_code = 0
      
      idxcode = 0
      idxtarget = 0
      do ii = 1, nlist
        item = trim(adjustl(target_list(ii)))
        litem = len_trim(item)
        is_code = item(1:1)=='<'.and.item(litem:litem)=='>'       
        if (is_code) then
          idxcode = idxcode+1
          itf_targets%code_name(idxcode) = item(2:litem-1)
        else
          idxtarget = idxtarget+1
          idxop = idx_oplist2(item,op_info)
          if (idxop.le.0) then
            write(lulog,*) 'Found:', item
            call quit(0,i_am,'Operator not on list')
          end if
          itf_targets%idx_target(idxtarget) = idxop
          itf_targets%idx_code(idxtarget) = idxcode
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'Entries of itf_targets:'
        write(lulog,*) 'ncodes   = ',itf_targets%ncodes
        write(lulog,*) 'ntargets = ',itf_targets%ntargets
        write(lulog,*) 'code_name:'
        do ii = 1, ncodes
          write(lulog,'(1x,i4,2x,a)') ii,trim(itf_targets%code_name(ii))
        end do
        write(lulog,'(1x,a,20i3)') 'idx_target: ',
     &       itf_targets%idx_target(1:ntargets)
        write(lulog,'(1x,a,20i3)') 'idx_code  : ',
     &       itf_targets%idx_code(1:ntargets)
      end if
      
      end subroutine
