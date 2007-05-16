*------------------------------------------------------------------------*
      real(8) function xnormop(ffop,op)
*------------------------------------------------------------------------*
*     get norm of operator
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
     &     idum, iblk
      real(8) ::
     &     xnrm2

      real(8), pointer ::
     &     buffer(:)

      real(8), external ::
     &     ddot

      ifree = mem_setmark('xnormop')

      if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
      else
        closeit = .false.
      end if

      xnrm2 = 0d0
      idxst = 1
      ! assuming that only the first x blocks *in a row* are
      ! incore !!
      if (ffop%buffered) then
        do iblk = 1, op%n_occ_cls
          if (ffop%incore(iblk)) then
            xnrm2 = xnrm2 +
     &           ddot(op%len_op_occ(iblk),
     &           ffop%buffer(op%off_op_occ(iblk)+1),1,
     &           ffop%buffer(op%off_op_occ(iblk)+1),1)
          else
            ! start reading file from here
            idxst = op%off_op_occ(iblk)+1
            exit
          end if
        end do
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

      ! to be fixed: should synchronize with record boundaries
      do while(idxst.le.len_op)
        idxnd = min(len_op,idxst-1+nbuff)
        call get_vec(ffop,buffer,idxst,idxnd)  
        xnrm2 = xnrm2 + ddot(idxnd-idxst+1,buffer,1,buffer,1)
        idxst = idxnd+1
      end do
      
      xnormop = sqrt(xnrm2)

      if (closeit)
     &     call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
