*----------------------------------------------------------------------*
      integer function ndisblk_mel2(mel,iblkst,iblknd)
*----------------------------------------------------------------------*
*     return number of distribution blocks in ME list
*     i.e. for each distribution,IRREP,MS,op_block
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'

      type(me_list), intent(in) ::
     &     mel
      integer, intent(in) ::
     &     iblkst, iblknd

      integer ::
     &     nblk, nsym, nmsblk, iblk, ms, isym, ist, ind

      integer, pointer ::
     &     ca_occ(:,:), ndis(:,:)

      nblk = mel%op%n_occ_cls
      nsym = mel%nsym
      ca_occ => mel%op%ica_occ

      if(iblkst.lt.1.or.iblknd.gt.nblk)
     &     call quit(1,'ndisblk_mel2','wrong number of blocks')

      if(iblkst.eq.-1)then
        ist = 1
      else
        ist = iblkst
      endif
      if(iblknd.eq.-1)then
        ind = nblk
      else
        ind = iblknd
      endif

      ndisblk_mel2 = 0
      do iblk = ist, ind
        nmsblk = min(ca_occ(1,iblk),ca_occ(2,iblk))+1
        ndis => mel%off_op_gmox(iblk)%ndis
        do ms = 1, nmsblk
          do isym = 1, nsym
            ndisblk_mel2 = ndisblk_mel2 +
     &           ndis(isym,ms)
          end do
        end do
      end do

      return
      end
