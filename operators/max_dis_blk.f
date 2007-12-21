*----------------------------------------------------------------------*
      integer function max_dis_blk(mode,mel,iblkop,orb_info)
*----------------------------------------------------------------------*
*     get longest distribution block of operator
*     if iblk>0, restrict to block iblk
*     mode=-1:  msa block
*     mode=0 :  msa,gma block
*     mode=1 :  dis,msa,gma block
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'

      type(me_list), intent(in) ::
     &     mel
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkop, mode

      integer ::
     &     iblk, iblk_min, iblk_max, lenblk, lenmsblk, ndis, idis,
     &     idxmsa, igama, msamax, ngam

      type(operator), pointer ::
     &     op

      op => mel%op

      if (iblkop.ge.0) then
        iblk_min = iblkop
        iblk_max = iblkop
      else
        iblk_min = 1
        iblk_max = op%n_occ_cls
      end if
        
      ngam = orb_info%nsym

      max_dis_blk = 0
      do iblk = iblk_min, iblk_max
        if (op%formal_blk(iblk)) cycle

        ! indeed, the number of MS(A) blocks is determined
        ! by the lower of the two occupations:
        msamax = min(op%ica_occ(1,iblk),op%ica_occ(2,iblk))
        
        do idxmsa = 1, msamax+1

          lenmsblk = 0
          do igama = 1, ngam
          
            lenblk = mel%len_op_gmo(iblk)%gam_ms(igama,idxmsa)
            if (lenblk.eq.0) cycle
            lenmsblk = lenmsblk + lenblk

            if (mode.eq.-1) cycle

            ndis = mel%off_op_gmox(iblk)%ndis(igama,idxmsa)

            if (mode.eq.0.or.ndis.eq.1) then
              max_dis_blk = max(max_dis_blk,lenblk)
              cycle
            end if
            
            do idis = 1, ndis
              max_dis_blk = max(max_dis_blk,mel%len_op_gmox(iblk)%
     &             d_gam_ms(idis,igama,idxmsa))
            end do

          end do

          if (mode.eq.-1)
     &         max_dis_blk = max(max_dis_blk,lenmsblk)

        end do
          
      end do

      return
      end

