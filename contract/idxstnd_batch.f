*----------------------------------------------------------------------*
      subroutine idxstnd_batch(len_batch,ibstc,ibsta,ibndc,ibnda,
     &                         idx_batch,
     &                         maxc,maxa,maxlen_batch)
*----------------------------------------------------------------------*
*     get start and end index for C/A batch #idx_batch, length len_batch
*     where maxc, maxa are the maximum lengthes for either C, A
*     initial call with ibndc, ibnda set to 0
*     C runs as leading dimension
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     idx_batch, maxc, maxa, maxlen_batch
      integer, intent(out) ::
     &     len_batch, ibstc, ibsta, ibndc, ibnda

      integer ::
     &     idxtotalst, idxtotalnd

      idxtotalst = (idx_batch-1)*maxlen_batch+1
      ibstc = mod(idxtotalst-1,maxc) + 1
      ibsta = (idxtotalst-1)/maxc + 1
      idxtotalnd = min(idxtotalst+maxlen_batch-1,maxc*maxa)
      ibndc = mod(idxtotalnd-1,maxc) + 1
      ibnda = (idxtotalnd-1)/maxc + 1
      len_batch = idxtotalnd - idxtotalst + 1

      return
      end
