*------------------------------------------------------------------------*
      real(8) function xnormop(mel)
*------------------------------------------------------------------------*
*     get norm of ME-list
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      type(me_list), intent(in) ::
     &     mel

      logical ::
     &     closeit
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
     &     idum, iblk, idisc_off, jblk
      real(8) ::
     &     xnrm2, xnrm_cur

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      real(8), pointer ::
     &     buffer(:)

      real(8), external ::
     &     ddot

      op => mel%op
      ffop => mel%fhand

      ifree = mem_setmark('xnormop')

      if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
      else
        closeit = .false.
      end if

      jblk = 0
      xnrm2 = 0d0
      idxst = 1
      ! assuming that only the first x blocks *in a row* are
      ! incore !!
      if (ffop%buffered) then
        do iblk = 1, op%n_occ_cls
          if (ffop%incore(iblk).gt.0) then
c            xnrm2 = xnrm2 +
            xnrm_cur =    ddot(mel%len_op_occ(iblk),
     &           ffop%buffer(mel%off_op_occ(iblk)+1),1,
     &           ffop%buffer(mel%off_op_occ(iblk)+1),1)
            xnrm2 = xnrm2 + xnrm_cur
c dbg
c            print *,'norm of block ','iblk',': ',sqrt(xnrm_cur)
c            jblk = iblk + 1
c dbgend
          else
            ! start reading file from here
            idxst = mel%off_op_occ(iblk)+1
c dbg
c            print *,'starting to read file from ',idxst
c            jblk = iblk - 1
c dbgend
            exit
          end if
        end do
c dbg
c      else
c        print *,'ME list not buffered!'
c dbgend
      end if

      nblkmax = ifree/ffop%reclen
      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffop%reclen
        call quit(1,'xnormop','not even 1 record fits into memory?')
      end if

      do while(jblk.le.op%n_occ_cls-1)
       jblk = jblk + 1
       idxst = mel%off_op_occ(jblk)+1

      len_op = mel%len_op_occ(jblk) !mel%len_op
      nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

      nbuff = min(len_op,nblk*ffop%reclen)

      ifree = mem_alloc_real(buffer,nbuff,'buffer')

      ! offset on file (if more than one instance of operator ex.)
      idisc_off = ffop%length_of_record*(ffop%current_record-1)

      ! to be fixed: should synchronize with record boundaries
      idxst = idxst + idisc_off
      xnrm_cur = 0d0
      do while(idxst.le.mel%off_op_occ(jblk)+len_op+idisc_off)
        idxnd = min(mel%off_op_occ(jblk)+len_op+idisc_off,idxst-1+nbuff)
        call get_vec(ffop,buffer,idxst,idxnd)  
        xnrm_cur = xnrm_cur + ddot(idxnd-idxst+1,buffer,1,buffer,1)
        idxst = idxnd+1
      end do
      xnrm2 = xnrm2 + xnrm_cur
c dbg
c      print *,'norm of block ',jblk,': ',sqrt(xnrm_cur)
c dbgend
      end do
      
      xnormop = sqrt(xnrm2)

      if (closeit)
     &     call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
