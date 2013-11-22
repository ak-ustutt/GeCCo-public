*----------------------------------------------------------------------*
      integer function idx_target(tgtname,tgt_info)
*----------------------------------------------------------------------*
*     given an target name, search tgt_info and return index of
*     of corresponding target
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     tgtname*(*)
      type(target_info), intent(in) ::
     &     tgt_info

      integer ::
     &     itgt

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'this is idx_target')
        write(lulog,*) ' looking for: "',trim(tgtname),'"'
      end if

      idx_target = -1
      do itgt = 1, tgt_info%ntargets
        if (trim(tgtname).eq.trim(tgt_info%array(itgt)%tgt%name)) then
          idx_target = itgt
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'result: ',idx_target
      end if

      return
      end
