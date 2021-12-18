*----------------------------------------------------------------------*
      subroutine set_disblkdim_mel(len_blk,ld_blk,mel)
*----------------------------------------------------------------------*
*     put len_op_gmox(<4>)%d_gam_ms(<1>,<2>,<3>)
*                        -> len_blk(<1>,<2>,<3>,<4>)
*          ld_op_gmox(<4>)%d_gam_ms(<1>,<2>,<3>)
*                         -> ld _blk(<1>,<2>,<3>,<4>)
*     i.e. for each distribution,IRREP,MS,op_block
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'


      type(me_list), intent(in) ::
     &     mel
      integer, intent(out) ::
     &     len_blk(*), ld_blk(*)

      integer ::
     &     nblk, nsym, nmsblk, iblk, ms, isym, idis, idx

      integer, pointer ::
     &     ca_occ(:,:), ndis(:,:),
     &     len_d_gam_ms(:,:,:),
     &     ld_d_gam_ms(:,:,:)

      nblk = mel%op%n_occ_cls
      nsym = mel%nsym
      ca_occ => mel%op%ica_occ

      idx = 0
      do iblk = 1, nblk
        nmsblk = min(ca_occ(1,iblk),ca_occ(2,iblk))+1
        ndis => mel%off_op_gmox(iblk)%ndis
        len_d_gam_ms => mel%len_op_gmox(iblk)%d_gam_ms
        ld_d_gam_ms  => mel%ld_op_gmox(iblk)%d_gam_ms
        do ms = 1, nmsblk
          do isym = 1, nsym
            do idis = 1, ndis(isym,ms)
              idx = idx+1
              len_blk(idx) = len_d_gam_ms(idis,isym,ms)
              ld_blk(idx)  = ld_d_gam_ms(idis,isym,ms)
            end do
          end do
        end do
      end do

      return
      end
