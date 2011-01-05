      subroutine scale_copy_op(label_res,label_inp,fac,nfac,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00
      
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, intent(in) ::
     &     nfac
      real(8), intent(in) ::
     &     fac(nfac)
      character(*), intent(in) ::
     &     label_res, label_inp
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info

      type(me_list), pointer ::
     &     me_res, me_inp
      type(filinf), pointer ::
     &     ffop_src, ffop_tgt

      integer ::
     &     idx_res, idx_inp, idx, idxnd_src, idxnd_tgt, len_op, nbuff,
     &     ipri, idxst_tgt, idxst_src, idisc_off_src, idisc_off_tgt,
     &     ifac, nblkmax, ifree, nblk
      logical ::
     &     open_close_res, open_close_inp,
     &     same
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      real(8), pointer ::
     &     buffer(:)

      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('scale_copy_op')

      if(ntest.ge.100)then
        write(luout,*) '===================='
        write(luout,*) ' scale & copy       '
        write(luout,*) '===================='
        write(luout,*) 'Result: ',trim(label_res)
        write(luout,*) 'Input:  ',trim(label_inp)
        write(luout,*) 'The factors (applied periodically): '
        do idx = 1, nfac
          write(luout,'(3x,f12.6)') fac(idx)
        end do
      endif

      idx_res = idx_mel_list(label_res,op_info)
      idx_inp = idx_mel_list(label_inp,op_info)

      if (idx_res.lt.0) then
        write(luout,*) '"',trim(label_res),'"'
        write(luout,*) idx_res
        call quit(1,'scale_copy_op','label not on list (1)')
      end if
      if (idx_inp.lt.0) then
        write(luout,*) '"',trim(label_inp),'"'
        write(luout,*) idx_inp
        call quit(1,'scale_copy_op','label not on list (2)')
      end if

      same = idx_res.eq.idx_inp
      if (same) call quit(1,'scale_copy_op',
     &       'copy means to have a second list...')

      ! Point to the relevant operators and their associated files.
      me_res => op_info%mel_arr(idx_res)%mel
      ffop_tgt => me_res%fhand
      if (.not.associated(ffop_tgt))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = ffop_tgt%unit.le.0

      if(open_close_res)then
        call file_open(ffop_tgt)
      endif

      me_inp => op_info%mel_arr(idx_inp)%mel
      ffop_src => me_inp%fhand
      if (.not.same) then
        if (.not.associated(ffop_src))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_inp%label))
        open_close_inp = ffop_src%unit.le.0
        if (open_close_inp) then
          call file_open(ffop_src)
        endif
      else
        open_close_inp = .false.
      end if

      ! record length hopefully the same
      if (ffop_tgt%reclen.ne.ffop_src%reclen)
     &   call quit(1,'scale_copy_op',
     &             'not prepared for different reclen''s')
      nblkmax = ifree/ffop_src%reclen
      if (nblkmax.le.0) then
        write(luout,*) 'free memory (words):  ',ifree
        write(luout,*) 'block length (words): ',ffop_src%reclen
        call quit(1,'scale_copy_op',
     &            'not even 1 record fits into memory?')
      end if

      len_op = me_inp%len_op
      ! result list may be longer but not shorter
      if (me_res%len_op.lt.len_op) call quit(1,'scale_copy_op',
     &         'Result list shorter than input')

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

        ifac = 1
        idxst_src = idisc_off_src+1
        idxst_tgt = idisc_off_tgt+1
        do while(idxst_src.le.idisc_off_src+len_op)
          idxnd_src = min(idisc_off_src+len_op,idxst_src-1+nbuff)
          idxnd_tgt = min(idisc_off_tgt+len_op,idxst_tgt-1+nbuff)
          call get_vec(ffop_src,buffer,idxst_src,idxnd_src)

          ! apply scaling factors (periodically)
          do idx = 1, idxnd_src-idxst_src+1
            buffer(idx) = fac(ifac)*buffer(idx)
c dbg
            if (buffer(idx).lt.-1d-14) then
              write(luout,*) 'setting to threshold el.# ',idx
              buffer(idx) = 0.1d0
            end if
c dbgend
            ifac = ifac + 1
            if (ifac.gt.nfac) ifac = 1
          end do

          call put_vec(ffop_tgt,buffer,idxst_tgt,idxnd_tgt)
          idxst_src = idxnd_src+1
          idxst_tgt = idxnd_tgt+1
        end do

      else

        call quit(1,'scale_copy_op','adapt for buffering')

      end if

      call touch_file_rec(ffop_tgt)

      ! needed: close loop over active records

      if (open_close_res)
     &     call file_close_keep(ffop_tgt)

      if (open_close_inp)
     &     call file_close_keep(ffop_src)

      if (ntest.ge.10) then
        write(luout,*) 'dump of scaled list:'
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
        call wrt_mel_file(luout,ipri,me_res,
     &       1,me_res%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      ifree = mem_flushmark()

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for scaling',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
