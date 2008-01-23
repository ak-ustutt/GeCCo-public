*----------------------------------------------------------------------*
      subroutine select_formula_target(idxselect,nselect,
     &     label,depend,op_info)
*----------------------------------------------------------------------*
*     select target "label" to be evaluated considering all targets
*     named in "depend" which need to be evaluated first
*
*     sets up idxselect and nselect
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_dependency_info.h'

      integer, intent(inout) ::
     &     idxselect(*), nselect

      character(*), intent(in) ::
     &     label
      type(dependency_info), intent(in) ::
     &     depend
      type(operator_info), intent(in) ::
     &     op_info

      integer, parameter ::
     &     ntest = 00
      integer ::
     &     idxmel, idxres, idxdep, idepend

      integer, external ::
     &     idxlist, idx_mel_list

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_formula_target')
        write(luout,*) 'label: ',trim(label)
        write(luout,*) 'nselect on input:   ',nselect
        write(luout,*) 'idxselect on input: ',idxselect(1:nselect)
      end if

      idxmel = idx_mel_list(label,op_info)
      if (idxmel.le.0)
     &     call quit(1,'select_formula_target',
     &     'label not on list: '//trim(label))
      
      idxres = idxlist(idxmel,depend%idxlist,depend%ntargets,1)

      if (idxres.le.0)
     &     call quit(1,'select_formula_target',
     &     'label not on dependency list: '//trim(label))

      ! already selected?
      if (nselect.gt.0.and.idxlist(idxres,idxselect,nselect,1).gt.0)
     &     return

      ! put on list
      nselect = nselect+1
      idxselect(nselect) = idxres
      
      ! check for dependencies:
      do idepend = 1, depend%ndepend
        idxmel = depend%depends_on_idxlist(idepend,idxres)
        if (idxmel.le.0) exit
        ! dependent on any other target named in "depend"?
        idxdep = idxlist(idxmel,depend%idxlist,depend%ntargets,1)
        if (idxdep.le.0) cycle

        ! already on list?
        if (idxlist(idxdep,idxselect,nselect,1).gt.0) cycle

        nselect = nselect+1
        idxselect(nselect) = idxdep

      end do

      if (ntest.ge.100) then
        write(luout,*) 'nselect on output:   ',nselect
        write(luout,*) 'idxselect on output: ',idxselect(1:nselect)
      end if

      return
      end
