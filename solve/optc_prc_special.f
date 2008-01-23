*----------------------------------------------------------------------*
      subroutine optc_prc_special(me_grd,me_special,nspecial,
     &                            nincore,xbuf1,xbuf2,xbuf3,lenbuf)
*----------------------------------------------------------------------*
*     experimental routine for testing special preconditioners
*     if nincore>1: xbuf1 contains the gradient vector on entry
*                   and should contain the preconditioned gradient
*                   on exit
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_me_list.h'


      integer, intent(in) ::
     &     nincore, lenbuf, nspecial
      type(me_list) ::
     &     me_grd
      type(me_list_array) ::
     &     me_special(nspecial)
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf), xbuf3(lenbuf)

      integer ::
     &     ngrd, nblk_grd, nbmat, nblk, idx, iblk, len
      type(me_list), pointer ::
     &     me_bmat

      integer, pointer ::
     &     len_blk(:), ld_blk(:)

      integer, external ::
     &     ndisblk_mel

c dbg
      print *,'in optc_prc_special'
c dbg
      if (nincore.ne.3)
     &     call quit(1,'optc_prc_special',
     &     'currently: special preconditioning for nincore==3, only')

      if (nspecial.lt.1)
     &     call quit(1,'optc_prc_special',
     &     'nspecial must be >= 1')

      ngrd     = me_grd%len_op
      nblk_grd = ndisblk_mel(me_grd)

      ! R12 code: here B^-1 was passed as me_special(1)%mel
      me_bmat => me_special(1)%mel
      nbmat = me_bmat%len_op
      ! should be open
      call vec_from_da(me_bmat%fhand,1,xbuf2,nbmat)

      ! Here I use two routines to get the subblock structure
      ! so that we can loop over these without caring too much
      ! about the details
      ! one might instead loop over MS(A), GAMMA(A) [, DISTRIBUTIONS]

      ! total number of sub blocks
      nblk = ndisblk_mel(me_bmat)

      ! here: the operator lengthes must match
      !  will be different for R12-triples etc.
      if (ngrd.ne.nbmat.and.nblk.ne.nblk_grd)
     &     call quit(1,'optc_prc_special',
     &     'block structures do not match')

      ! lengthes and leading dimensions of sub-blocks
      allocate(len_blk(nblk),ld_blk(nblk))
      call set_disblkdim_mel(len_blk,ld_blk,me_bmat)

      idx = 1
      do iblk = 1, nblk
        len = ld_blk(iblk)
        call dgemm('n','n',len,len,len,
     &           1d0,xbuf2(idx),len,
     &               xbuf1(idx),len,
     &           0d0,xbuf3(idx),len)
        idx = idx + len_blk(iblk)
      end do

      deallocate(len_blk,ld_blk)

      return
      end


