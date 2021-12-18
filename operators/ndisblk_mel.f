*----------------------------------------------------------------------*
      integer function ndisblk_mel(mel)
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

      integer ::
     &     nblk, nsym, nmsblk, iblk, ms, isym

      integer, pointer ::
     &     ca_occ(:,:), ndis(:,:)

      nblk = mel%op%n_occ_cls
      nsym = mel%nsym
      ca_occ => mel%op%ica_occ

      ndisblk_mel = 0
      do iblk = 1, nblk
        nmsblk = min(ca_occ(1,iblk),ca_occ(2,iblk))+1
        ndis => mel%off_op_gmox(iblk)%ndis
        do ms = 1, nmsblk
          do isym = 1, nsym
            ndisblk_mel = ndisblk_mel +
     &           ndis(isym,ms)
          end do
        end do
      end do

      return
      end
