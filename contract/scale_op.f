      subroutine scale_op(label_res,idx_blk,fac,label_inp,nblk,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*     scale blocks of operator by factor fac,
*     if nblk==-1, all blocks are scaled with the same factor
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100
      
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, intent(in) ::
     &     nblk, idx_blk(*)
      real(8), intent(in) ::
     &     fac(*)
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

      integer ::
     &     idx_res, idx_inp, idx, idxnd,
     &     njoined, iblk, ipri
      logical ::
     &     open_close_res, open_close_inp, same
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0, xnorm2
      integer, pointer  ::
     &     occ(:,:,:)

      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      if(ntest.ge.100)then
        write(luout,*) '===================='
        write(luout,*) ' scale operators   '
        write(luout,*) '===================='
        write(luout,*) 'Result: ',trim(label_res)
        write(luout,*) 'The factors and summands: '
        do idx = 1, nblk
          write(luout,'(3x,f12.6,x,i8)') fac(idx),idx_blk(idx)
        end do
        if (nblk.lt.0) then
          write(luout,'(3x,f12.6,x,a)') fac(idx),'applied to all blocks'
        end if
      endif

      idx_res = idx_mel_list(label_res,op_info)
      idx_inp = idx_mel_list(label_inp,op_info)

      if (idx_res.lt.0) then
        write(luout,*) '"',trim(label_res),'"'
        write(luout,*) idx_res
        call quit(1,'scale_op','label not on list')
      end if
      if (idx_inp.lt.0) then
        write(luout,*) '"',trim(label_inp),'"'
        write(luout,*) idx_inp
        call quit(1,'scale_op','label not on list')
      end if

      same = idx_res.eq.idx_inp

      ! Point to the relevant operators and their associated files.
      me_res => op_info%mel_arr(idx_res)%mel
      if (.not.associated(me_res%fhand))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = me_res%fhand%unit.le.0

      if(open_close_res)then
        call file_open(me_res%fhand)
      endif

      me_inp => op_info%mel_arr(idx_inp)%mel
      if (.not.same) then
        if (.not.associated(me_inp%fhand))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_inp%label))
        open_close_inp = me_inp%fhand%unit.le.0
        if(open_close_inp)then
          call file_open(me_res%fhand)
        endif
      else
        open_close_inp = .false.
      end if

      idxnd = nblk
      if (nblk.lt.0) idxnd = me_inp%op%n_occ_cls
      njoined = me_res%op%njoined

      do idx = 1, idxnd
        if (nblk.gt.0) then
          iblk = idx_blk(idx)
        else
          iblk = idx
        end if

        call scale_opblk(xnorm2,fac(idx),me_inp,me_res,
     &       iblk,iblk,orb_info)

      end do

      if (open_close_res)
     &     call file_close_keep(me_res%fhand)

      if (open_close_inp)
     &     call file_close_keep(me_inp%fhand)

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

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for scaling',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
