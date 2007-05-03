*------------------------------------------------------------------------*
      subroutine zeroop(ffop,op)
*------------------------------------------------------------------------*
*     set operator to zero, initialize with zero
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'

      type(operator), intent(in) ::
     &     op
      type(filinf), intent(inout) ::
     &     ffop

      logical ::
     &     closeit
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
     &     idum

      real(8), pointer ::
     &     buffer(:)

      ifree = mem_setmark('zeroop')

      if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
      else
        closeit = .false.
      end if

      nblkmax = ifree/ffop%reclen
      if (nblkmax.le.0) then
        write(luout,*) 'free memory (words):  ',ifree
        write(luout,*) 'block length (words): ',ffop%reclen
        call quit(1,'zeroop','not even 1 record fits into memory?')
      end if

      len_op = op%len_op
      nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

      nbuff = min(len_op,nblk*ffop%reclen)

      ifree = mem_alloc_real(buffer,nbuff,'buffer')

      buffer(1:nbuff) = 0d0

      idxst = 1
      do while(idxst.le.len_op)
        idxnd = min(len_op,idxst-1+nbuff)
        call put_vec(ffop,buffer,idxst,idxnd)  
        idxst = idxnd+1
      end do
      
      if (closeit)
     &     call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
