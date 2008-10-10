*----------------------------------------------------------------------*
      subroutine set_op_scratch(len_buffer,buffering_type,
     &                            me_op,iblk,maxscr,
     &                            orb_info)
*----------------------------------------------------------------------*
*     return the recommended buffering type and the necessary buffer 
*     length for ME list me_op, when a length of maxscr is intended
*     
*     types:
*         0 -- full list
*         1 -- Ms blocks
*         2 -- Ms/gamma blocks
*         3 -- Ms/gamma/dist blocks   DISABLED
*         4 -- even smaller batches   DISABLED
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_me_list.h'
      include 'def_orbinf.h'

      integer, intent(in) ::
     &     iblk, maxscr
      integer, intent(out) ::
     &     len_buffer, buffering_type
      type(me_list), intent(in) ::
     &     me_op
      type(orbinf), intent(in) ::
     &     orb_info

      integer, external ::
     &     max_dis_blk

      buffering_type=0
      len_buffer = me_op%len_op_occ(iblk)
      if (len_buffer.le.maxscr) return

      buffering_type=1
      len_buffer = max_dis_blk(-1,me_op,iblk,orb_info)
      if (len_buffer.le.maxscr) return

      buffering_type=2
      len_buffer =  max_dis_blk(0,me_op,iblk,orb_info)

      return
      ! disabled:

      if (len_buffer.le.maxscr) return

      buffering_type=3
      len_buffer = max_dis_blk(1,me_op,iblk,orb_info)
      if (len_buffer.le.maxscr) return

      buffering_type=4

      return
      end
