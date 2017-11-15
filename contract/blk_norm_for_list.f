*------------------------------------------------------------------------*
      subroutine blk_norm_for_list(xnorm,nblk_in,melist_a,melist_b)
*------------------------------------------------------------------------*
*     get norm^2 (=weight) of blocks of ME-list melist_a
*     if the second list melist_b is not equal to melist_a, 
*     the subroutine will assume that you want <A|B>) 
*     (needed f.x. if a proper norm requires a metric: pass metric x list
*     in one of the two lists)
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      character(len=*),parameter::
     &     i_am = 'blk_norm_for_list'
      type(me_list), intent(in) ::
     &     melist_a, melist_b
      integer, intent(in) ::
     &     nblk_in
      real(8), intent(out) ::
     &     xnorm(nblk_in)

      logical ::
     &     closeit,close_b,twolists
      integer ::
     &     ifree, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
     &     idum, iblk, idisc_off, jblk, idisc_off_b, idxst_b, idxnd_b
      real(8) ::
     &     xnrm2, xnrm_cur

      type(operator), pointer ::
     &     opa, opb
      type(filinf), pointer ::
     &     ffopa, ffopb
      real(8), pointer ::
     &     buffera(:), bufferb(:)

      real(8), external ::
     &     ddot

      opa => melist_a%op
      ffopa => melist_a%fhand
      opb => melist_b%op
      ffopb => melist_b%fhand

      ifree = mem_setmark('xnormop')

      if (ffopa%unit.le.0) then
        call file_open(ffopa)
        closeit = .true.
      else
        closeit = .false.
      end if

      twolists = (melist_a%label.ne.melist_b%label)

      ! some sanity checks:
      if (nblk_in.lt.opa%n_occ_cls) then
        write(luout,*) 'nblk too small!'
        write(luout,*) 'provided: ',nblk_in
        write(luout,*) 'needed:   ',melist_a%op%n_occ_cls
        call quit(1,'blk_norm_for_list','nblk_in too small!')
      end if

      if (twolists.and.opa%n_occ_cls.ne.opb%n_occ_cls) then
        call quit(1,'blk_norm_for_list','incompatible lists')
      end if

      close_b = .false.
      if (twolists) then
        if (ffopb%unit.le.0) then
          call file_open(ffopb)
          close_b = .true.
        end if
      end if

      jblk = 0
      xnrm2 = 0d0
      idxst = 1
      if (ffopa%buffered.and.ffopb%buffered) then
        do iblk = 1, opa%n_occ_cls
          if (ffopa%incore(iblk).gt.0.and.ffopb%incore(iblk).gt.0) then
c            xnrm2 = xnrm2 +
            if (melist_a%len_op_occ(iblk).ne.melist_b%len_op_occ(iblk))
     &          call quit(1,'blk_norm_for_list','incompatible lists!')
            xnrm_cur =    ddot(melist_a%len_op_occ(iblk),
     &           ffopa%buffer(melist_a%off_op_occ(iblk)+1),1,
     &           ffopb%buffer(melist_b%off_op_occ(iblk)+1),1)
            xnorm(jblk) = xnrm_cur !sqrt(xnrm_cur)
            xnrm2 = xnrm2 + xnrm_cur
          else if((ffopa%incore(iblk).gt.0.and.ffopb%incore(iblk).eq.0)
     &        .or.(ffopa%incore(iblk).eq.0.and.ffopb%incore(iblk).gt.0))
     &    then
            call quit(1,'blk_norm_for_list',
     &                'not prepared for different buffering of A and B') 
          else
            ! start reading file from here
            idxst = melist_a%off_op_occ(iblk)+1
            jblk = iblk - 1
            exit
          end if
        end do
      end if

      nblkmax = ifree/ffopa%reclen
      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffopa%reclen
        call quit(1,'blk_norm_for_list',
     &              'not even 1 record fits into memory?')
      end if

      do while(jblk.le.opa%n_occ_cls-1)
        jblk = jblk + 1
        idxst = melist_a%off_op_occ(jblk)+1

        if (twolists) idxst_b = melist_b%off_op_occ(jblk)+1

        len_op = melist_a%len_op_occ(jblk) !mel%len_op
        nblk = min((len_op-1)/ffopa%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffopa%reclen)

        ifree = mem_alloc_real(buffera,nbuff,'buffer')

        if (twolists.and.len_op.ne.melist_b%len_op_occ(jblk))
     &    call quit(1,'blk_norm_for_list','incompatible lists!')

        if (twolists) ifree = mem_alloc_real(bufferb,nbuff,'buffer2')

        ! offset on file (if more than one instance of operator ex.)
        idisc_off = ffopa%length_of_record*(ffopa%current_record-1)
        if (twolists)
     &   idisc_off_b = ffopb%length_of_record*(ffopb%current_record-1)

        
        ! to be fixed: should synchronize with record boundaries
        if (.not.twolists) then
          idxst = idxst + idisc_off
          xnrm_cur = 0d0
          do while(idxst.le.melist_a%off_op_occ(jblk)+len_op+idisc_off)
            idxnd = min(melist_a%off_op_occ(jblk)+len_op+idisc_off,
     &                  idxst-1+nbuff)
            call get_vec(ffopa,buffera,idxst,idxnd)  
            xnrm_cur = xnrm_cur+ddot(idxnd-idxst+1,buffera,1,buffera,1)
            idxst = idxnd+1
          end do
        else
          idxst = idxst + idisc_off
          idxst_b = idxst_b + idisc_off_b
          xnrm_cur = 0d0
          do while(idxst.le.melist_a%off_op_occ(jblk)+len_op+idisc_off)
            idxnd =   min(melist_a%off_op_occ(jblk)+len_op+idisc_off,
     &                    idxst-1+nbuff)
            idxnd_b = min(melist_b%off_op_occ(jblk)+len_op+idisc_off,
     &                    idxst_b-1+nbuff)
            call get_vec(ffopa,buffera,idxst,idxnd)
            call get_vec(ffopb,bufferb,idxst_b,idxnd_b)
            xnrm_cur = xnrm_cur+ddot(idxnd-idxst+1,buffera,1,bufferb,1)
            idxst = idxnd+1
            idxst_b = idxnd_b+1
          end do
        end if
        xnorm(jblk) = xnrm_cur !sqrt(xnrm_cur)
        xnrm2 = xnrm2 + xnrm_cur
      end do
      
      !xnormop = sqrt(xnrm2)

      if (closeit)
     &     call file_close_keep(ffopa)
      if (close_b)
     &     call file_close_keep(ffopb)

      ifree = mem_flushmark()

      return
      end

