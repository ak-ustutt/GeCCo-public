*----------------------------------------------------------------------*
      subroutine mel_export(ffarch,mel,orb_info)
*----------------------------------------------------------------------*
*     store contents of ME list in archive format, together with
*     info on orbital space etc.
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_me_list.h'
      include 'def_orbinf.h'

      type(filinf), intent(in) ::
     &     ffarch
      type(me_list), intent(in) ::
     &     mel
      type(orbinf), intent(in) ::
     &     orb_info

      ! open,rewind
      call file_open(ffarch)
      rewind(ffarch%unit)

      ! store orbital info
      call store_orbinf(ffarch,orb_info)

      ! store definition of operator
      call store_opdef(ffarch,mel%op,orb_info)

      ! store definition of ME-list
      call store_meldef(ffarch,mel,orb_info)

      ! store matrix elements
      call store_mel_val(ffarch,mel,orb_info)

      call file_close(ffarch)

      return
      end
