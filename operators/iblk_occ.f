*----------------------------------------------------------------------*
      integer function iblk_occ(iocc,dag,oper,version)
*----------------------------------------------------------------------*
*     return block of operator oper that corresponds to the
*     hpv-occupation iocc
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      logical, intent(in) ::
     &     dag
      integer, intent(in) ::
     &     iocc(ngastp,2,*), version
      type(operator) ::
     &     oper

      logical ::
     &     dago
      integer ::
     &     idx, ijoin, njoined, iblkoff, nblk
      integer, pointer ::
     &     op_occ(:,:,:)
      logical, external ::
     &     iocc_equal_n

      iblk_occ = -1
      dago = oper%dagger
      njoined = oper%njoined
      nblk = oper%n_occ_cls

      op_occ => oper%ihpvca_occ(:,:,:)

      idx_loop: do idx = 1, nblk

        iblkoff = (idx-1)*njoined
        if (iocc_equal_n(iocc,dag,
     &                   op_occ(1,1,iblkoff+1),dago,njoined)
     &     .and.(version.eq.oper%blk_version(idx).or.version.eq.0)) then
          iblk_occ = idx
          exit idx_loop
        end if

      end do idx_loop

      return

      end
