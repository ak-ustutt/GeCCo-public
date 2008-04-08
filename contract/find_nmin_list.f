*------------------------------------------------------------------------*
      subroutine find_nmin_list(xlist,idxlist,nlist,mel)
*------------------------------------------------------------------------*
*     search ME-list and return the nlist lowest values + indices
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     nlist
      integer, intent(out) ::
     &     idxlist(nlist)
      real(8), intent(out) ::
     &     xlist(nlist)

      logical ::
     &     closeit, init
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
     &     idum, iblk, idisc_off, kdx

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      real(8), pointer ::
     &     buffer(:)

      ifree = mem_setmark('find_nmin_list')

      ffop => mel%fhand
      op => mel%op

      if (.not.associated(ffop))
     &     call quit(1,'find_nmin_list','No file assigned to list: '//
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
        call quit(1,'find_nmin_list',
     &       'not even 1 record fits into memory?')
      end if

      ! offset on file (if more than one instance of operator ex.)
      idisc_off = ffop%length_of_record*(ffop%current_record-1)

      if (.not.ffop%buffered) then

        len_op = mel%len_op
        nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffop%reclen)

        ifree = mem_alloc_real(buffer,nbuff,'buffer')

        init = .true.
        idxst = idisc_off+1
        do while(idxst.le.idisc_off+len_op)
          idxnd = min(idisc_off+len_op,idxst-1+nbuff)
          call get_vec(ffop,buffer,idxst,idxnd)  
          call find_nmin(xlist,idxlist,nlist,buffer,
     &         idxnd-idxst+1,idxst-idisc_off-1,init)
          if (ntest.ge.100) then
            write(luout,*) 'after batch: ',idxst-idisc_off,
     &                              ' to ',idxnd-idisc_off
            do kdx = 1, nlist
              write(luout,*) kdx,xlist(kdx),idxlist(kdx)
            end do
          end if
          idxst = idxnd+1
          init = .false.
        end do
        
      else

        call quit(1,'find_nmin_list','incore part not tested yet')

        ! search the buffer (= all blocks which are incore)
          call find_nmin(xlist,idxlist,nlist,ffop%buffer,
     &         ffop%nbuffer,0,.true.)

        ! search remaining blocks on disc
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
                call get_vec(ffop,buffer,idxst,idxnd)  
                call find_nmin(xlist,idxlist,nlist,buffer,
     &               idxnd-idxst+1,idxst-idisc_off-1,.false.)
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
