*----------------------------------------------------------------------*
      subroutine add_unity_drv(label,fac,init,minblk_in,maxblk_in,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      character(len=*), intent(in) ::
     &     label
      real(8), intent(in) ::
     &     fac
      logical, intent(in) ::
     &     init
      integer, intent(in) ::
     &     minblk_in, maxblk_in
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(orbinf) ::
     &     orb_info

      integer ::
     &     idx, jdx, minblk, maxblk

      integer, external ::
     &     idx_mel_list

      type(me_list), pointer ::
     &     mel_pnt

      idx = idx_mel_list(label,op_info)
      if (idx.lt.0)
     &       call quit(1,'process_me_lists','Label not on list: "'//
     &       trim(label)//'"')

      mel_pnt => op_info%mel_arr(idx)%mel

      minblk = 1
      if (minblk_in.gt.0) minblk = minblk_in
      maxblk = mel_pnt%op%n_occ_cls
      if (maxblk_in.gt.0) maxblk = maxblk_in

      if (init) call zeroop(mel_pnt)

      do idx = minblk, maxblk
C???        jdx = (idx-1)*mel_pnt%op%njoined+1
C???        call add_unity(fac,mel_pnt,jdx,orb_info,str_info)
        call add_unity(fac,mel_pnt,idx,orb_info,str_info)
      end do

      end
