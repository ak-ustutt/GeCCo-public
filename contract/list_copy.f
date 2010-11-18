*------------------------------------------------------------------------*
      subroutine list_copy(me_src,me_tgt)
*------------------------------------------------------------------------*
*     cp ME list from me_src to me_tgt
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      type(me_list), intent(inout) ::
     &     me_src, me_tgt

      logical ::
     &     close_src, close_tgt
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff,
     &     idxst_src, idxnd_src, idxst_tgt, idxnd_tgt,
     &     idum, iblk, idisc_off_src, idisc_off_tgt

      type(operator), pointer ::
     &     op_src, op_tgt
      type(filinf), pointer ::
     &     ffop_src, ffop_tgt
      real(8), pointer ::
     &     buffer(:)

cmh      call quit(1,'list_copy','not yet debugged (*might* work)')

      ifree = mem_setmark('list_copy')

      ffop_src => me_src%fhand
      ffop_tgt => me_tgt%fhand
      op_src => me_src%op
      op_tgt => me_tgt%op

      if (.not.associated(ffop_src))
     &     call quit(1,'list_copy','No file assigned to list: '//
     &     trim(me_src%label))
      if (.not.associated(ffop_tgt))
     &     call quit(1,'list_copy','No file assigned to list: '//
     &     trim(me_tgt%label))

      if (ffop_src%unit.le.0) then
        call file_open(ffop_src)
        close_src = .true.
      else
        close_src = .false.
      end if
      if (ffop_tgt%unit.le.0) then
        call file_open(ffop_tgt)
        close_tgt = .true.
      else
        close_tgt = .false.
      end if

      if (ffop_src%reclen.ne.ffop_tgt%reclen)
     &     call quit(1,'list_copy',
     &     'not prepared for different reclen''s: '//
     &     trim(ffop_src%name)//' '//trim(ffop_tgt%name))
      nblkmax = ifree/ffop_src%reclen
      if (nblkmax.le.0) then
        write(luout,*) 'free memory (words):  ',ifree
        write(luout,*) 'block length (words): ',ffop_src%reclen
        call quit(1,'list_copy','not even 1 record fits into memory?')
      end if

      len_op = me_src%len_op
      if (len_op.ne.me_tgt%len_op)
     &     call quit(1,'list_copy',
     &     'incompatible lists: '//
     &     trim(me_src%label)//' '//trim(me_tgt%label))

      ! offset on file (if more than one instance of list exists)
      idisc_off_src =
     &     ffop_src%length_of_record*(ffop_src%current_record-1)
      idisc_off_tgt =
     &     ffop_tgt%length_of_record*(ffop_tgt%current_record-1)

      if (.not.ffop_src%buffered.and.
     &    .not.ffop_tgt%buffered) then

        nblk = min((len_op-1)/ffop_src%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffop_src%reclen)

        ifree = mem_alloc_real(buffer,nbuff,'buffer')

        idxst_src = idisc_off_src+1
        idxst_tgt = idisc_off_tgt+1
        do while(idxst_src.le.idisc_off_src+len_op)
          idxnd_src = min(idisc_off_src+len_op,idxst_src-1+nbuff)
          idxnd_tgt = min(idisc_off_tgt+len_op,idxst_tgt-1+nbuff)
          call get_vec(ffop_src,buffer,idxst_src,idxnd_src)  
          call put_vec(ffop_tgt,buffer,idxst_tgt,idxnd_tgt)  
          idxst_src = idxnd_src+1
          idxst_tgt = idxnd_tgt+1
        end do
        
      else     
        
        call quit(1,'list_copy','adapt for buffering')
        
      end if

      call touch_file_rec(ffop_tgt)

      if (close_src)
     &     call file_close_keep(ffop_src)
      if (close_tgt)
     &     call file_close_keep(ffop_tgt)

      ifree = mem_flushmark()

      return
      end
