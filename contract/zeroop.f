*------------------------------------------------------------------------*
!>     set ME-list to zero 
!!
!!     @param mel ME-List to be zeroed
!!
!!     prepared for buffered
*------------------------------------------------------------------------*
      subroutine zeroop(mel)
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      type(me_list), intent(inout) ::
     &     mel

      logical ::
     &     closeit
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
     &     idum, iblk, idisc_off

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      real(8), pointer ::
     &     buffer(:)

      ifree = mem_setmark('zeroop')

      ffop => mel%fhand
      op => mel%op

      if (.not.associated(ffop))
     &     call quit(1,'zeroop','No file assigned to list: '//
     &     trim(mel%label))

      if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
      else
        closeit = .false.
      end if

      nblkmax = ifree/ffop%reclen
      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffop%reclen
        call quit(1,'zeroop','not even 1 record fits into memory?')
      end if

      ! offset on file (if more than one instance of operator ex.)
      idisc_off = ffop%length_of_record*(ffop%current_record-1)

      if (.not.ffop%buffered) then

        len_op = mel%len_op
        nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffop%reclen)

        ifree = mem_alloc_real(buffer,nbuff,'buffer')

        buffer(1:nbuff) = 0d0

        idxst = idisc_off+1
        do while(idxst.le.idisc_off+len_op)
          idxnd = min(idisc_off+len_op,idxst-1+nbuff)
          call put_vec(ffop,buffer,idxst,idxnd)  
          idxst = idxnd+1
        end do
        
      else

        ! zero the buffer (= all blocks which are incore)
        ffop%buffer(1:ffop%nbuffer) = 0d0

        ! zero all blocks on disc
        len_op = 0 ! look for largest block
        do iblk = 1, op%n_occ_cls
          if (ffop%incore(iblk).le.0) 
     &         len_op = max(len_op,mel%len_op_occ(iblk))
        end do

        if (len_op.gt.0) then
          nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

          nbuff = min(len_op,nblk*ffop%reclen)
          ifree = mem_alloc_real(buffer,nbuff,'buffer')

          buffer(1:nbuff) = 0d0

          do iblk = 1, op%n_occ_cls
            if (ffop%incore(iblk).le.0) then
              len_op = mel%len_op_occ(iblk)
              idxst = idisc_off+mel%off_op_occ(iblk)+1
              len_op = idxst-1+len_op
              do while(idxst.le.len_op)
                idxnd = min(len_op,idxst-1+nbuff)
                call put_vec(ffop,buffer,idxst,idxnd)  
                idxst = idxnd+1
              end do
            end if
          end do
        end if

      end if

      if (closeit)
     &     call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
