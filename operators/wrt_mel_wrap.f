*----------------------------------------------------------------------*
      subroutine wrt_mel_buf(luout,level,bufmel,mel,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for wrt_mel: operator elements on incore buffer
*     see below for further info
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, level, iblkst, iblknd
      real(8), intent(in) ::
     &     bufmel(*)
      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(filinf) ::
     &     ffdum

      call wrt_mel(luout,level,.true.,bufmel,mel,iblkst,iblknd,
     &     str_info,orb_info)

      return
      end

*----------------------------------------------------------------------*
      subroutine wrt_mel_file(luout,level,mel,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for wrt_mel: matrix-elements on file
*     see below for further info
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, level, iblkst, iblknd
      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      real(8) ::
     &     xdum

      call wrt_mel(luout,level,.false.,xdum,mel,iblkst,iblknd,
     &     str_info,orb_info)

      return
      end

