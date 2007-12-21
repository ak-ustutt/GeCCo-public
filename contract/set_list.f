*------------------------------------------------------------------------*
      subroutine set_list(mel,idxset,valset,nset)
*------------------------------------------------------------------------*
*     initialize ME-list with zero and set the elements given
*     in idxset(1:nset) to the values valset(1:nset)
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      integer ::
     &     ntest = 100

      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     nset, idxset(nset)
      real(8), intent(in) ::
     &     valset(nset)
      
      logical ::
     &     closeit, zero_buff
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
     &     idum, iblk, idisc_off, iset

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      real(8), pointer ::
     &     buffer(:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_list')
        write(luout,*) 'nset = nset'
        do iset = 1, nset
          write(luout,'(3x,i10,g20.12)') idxset(iset), valset(iset)
        end do
      end if

      ifree = mem_setmark('set_list')

      ffop => mel%fhand
      op => mel%op

      if (.not.associated(ffop))
     &     call quit(1,'set_list','No file assigned to list: '//
     &     trim(mel%label))

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
        call quit(1,'set_list','not even 1 record fits into memory?')
      end if

      ! offset on file (if more than one instance of operator ex.)
      idisc_off = ffop%length_of_record*(ffop%current_record-1)

      if (.not.ffop%buffered) then

        len_op = mel%len_op
        nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffop%reclen)

        ifree = mem_alloc_real(buffer,nbuff,'buffer')

        buffer(1:nbuff) = 0d0

        iset = 1
        idxst = idisc_off+1
        zero_buff = .false.
        do while(idxst.le.idisc_off+len_op)
          if (zero_buff) buffer(1:nbuff) = 0d0
          zero_buff = .false.
          idxnd = min(idisc_off+len_op,idxst-1+nbuff)
          do while (idxset(iset).le.idxnd-idisc_off.and.iset.le.nset)
            buffer(idxset(iset)-idxst+idisc_off+1) = valset(iset)
            iset = iset+1
            zero_buff = .true.
          end do
          call put_vec(ffop,buffer,idxst,idxnd)  
          idxst = idxnd+1
        end do
        
      else

        call quit(1,'set_list','no incore part yet')
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
