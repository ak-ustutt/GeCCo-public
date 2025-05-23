
!>    initialize ME-list with zero and set certain elements.
!!
!!    initalizes the ME-list to zero and sets the elements given
!!    in idxset(1:nset) to the values valset(1:nset)
!!    @param[inout] mel ME-list to be worked on
!!    @param[in] idxset list of indexes if first value negative:
!!                  cycle though valset for the whole operator
!!    @param[in] valset list of values that are written to the corresponding indices
!!    @param[in] nset length of idxset and valset


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
     &     ntest = 00

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
     &     idum, iblk, idisc_off, iset, idx

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      real(8), pointer ::
     &     buffer(:)

c      call get_argument_value('general','print',ival=iprlvl)
c      ntest=max(ntest,iprlvl)

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'set_list')
        write(lulog,*) 'nset = nset'
        do iset = 1, nset
          write(lulog,'(3x,i10,g20.12)') idxset(iset), valset(iset)
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
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffop%reclen
        call quit(1,'set_list','not even 1 record fits into memory?')
      end if

      ! offset on file (if more than one instance of operator ex.)
      idisc_off = ffop%length_of_record*(ffop%current_record-1)

c should work for simple cases
      if (.not.ffop%buffered.or.mel%len_op.eq.1) then

        len_op = mel%len_op
        nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffop%reclen)

        ifree = mem_alloc_real(buffer,nbuff,'buffer')

        buffer(1:nbuff) = 0d0

        ! all elements set to the same value?
        idx = 0
        if (idxset(1).lt.0) idx = 1

        iset = 1
        idxst = idisc_off+1
        zero_buff = .false.
        do while(idxst.le.idisc_off+len_op)
          if (zero_buff) buffer(1:nbuff) = 0d0
          zero_buff = .false.
          idxnd = min(idisc_off+len_op,idxst-1+nbuff)
          do while (max(idxset(iset),idx).le.idxnd-idisc_off
     &              .and.iset.le.nset)
            buffer(max(idxset(iset),idx)-idxst+idisc_off+1)
     &            = valset(iset)
            iset = iset+1
            if (idx.gt.0) then
              idx = idx+1
              if (iset.ge.nset) iset = 1
            end if
            zero_buff = .true.
          end do
          call put_vec(ffop,buffer,idxst,idxnd)  
          idxst = idxnd+1
        end do
        
      else

        call quit(1,'set_list','no incore part yet')
c        ! zero the buffer (= all blocks which are incore)
c        ffop%buffer(1:ffop%nbuffer) = 0d0
c
c        ! zero all blocks on disc
c        len_op = 0 ! look for largest block
c        do iblk = 1, op%n_occ_cls
c          if (ffop%incore(iblk).le.0) 
c     &         len_op = max(len_op,mel%len_op_occ(iblk))
c        end do

c        if (len_op.gt.0) then
c          nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

c          nbuff = min(len_op,nblk*ffop%reclen)
c          ifree = mem_alloc_real(buffer,nbuff,'buffer')

c          buffer(1:nbuff) = 0d0
c
c          do iblk = 1, op%n_occ_cls
c            if (ffop%incore(iblk).le.0) then
c              len_op = mel%len_op_occ(iblk)
c              idxst = idisc_off+mel%off_op_occ(iblk)+1
c              len_op = idxst-1+len_op
c              do while(idxst.le.len_op)
c                idxnd = min(len_op,idxst-1+nbuff)
c                call put_vec(ffop,buffer,idxst,idxnd)  
c                idxst = idxnd+1
c              end do
c            end if
c          end do
c        end if

      end if

      call touch_file_rec(ffop)

      if (closeit)
     &     call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
