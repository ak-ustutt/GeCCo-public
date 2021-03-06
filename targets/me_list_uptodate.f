*----------------------------------------------------------------------*
      logical function me_list_uptodate(idx_res,depend_info,op_info)
*----------------------------------------------------------------------*
*     run set_formula_dependencies() first
*     check, whether the ME-list referenced as entry #idx_res
*     in depend_info is newer than the ME-list on which it depends
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_dependency_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     idx_res
      type(dependency_info), intent(in) ::
     &     depend_info
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     uptodate
      integer ::
     &     irec, idx, jdx, mxidx
      integer(8) ::
     &     time_res, time_int

      idx = depend_info%idxlist(idx_res)
      irec = op_info%mel_arr(idx)%mel%fhand%current_record
      time_res = op_info%mel_arr(idx)%mel%fhand%last_mod(irec)

      if (ntest.ge.100)
     &     write(lulog,*) 'me_list_uptodate: considering: ',
     &     trim(op_info%mel_arr(idx)%mel%label),' time mark',time_res

      if (time_res.lt.0) then
        me_list_uptodate = .false.
        return
      end if

      uptodate = .true.
      mxidx = depend_info%ndepend

      do idx = 1, mxidx
        jdx = depend_info%depends_on_idxlist(idx,idx_res)
        if (jdx.le.0) exit
        irec = op_info%mel_arr(jdx)%mel%fhand%current_record

        if (ntest.ge.100)
     &     write(lulog,*) ' comparing: ',
     &     trim(op_info%mel_arr(jdx)%mel%label),' time mark',
     &       op_info%mel_arr(jdx)%mel%fhand%last_mod(irec)

        time_int = op_info%mel_arr(jdx)%mel%fhand%last_mod(irec)
cmh   following should be commented in, but currently not possible
cmh   due to dirty tricks in evpc_core
c        if (time_int.lt.0) call quit(1,'me_list_uptodate',
c     &        'trying to use list never created: '//
c     &        trim(op_info%mel_arr(jdx)%mel%label))

        uptodate = uptodate.and.time_int.le.time_res
        if (.not.uptodate) exit
      end do

      if (ntest.ge.100) write(lulog,*) ' -> uptodate: ',uptodate

      me_list_uptodate = uptodate

      return
      end
