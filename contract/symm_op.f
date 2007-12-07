      subroutine symm_op(idx_opin,idx_opout,op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*
*     Routine to set up the symmetrisation of an operator matrix.
*     GWR November 2007
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_operator_array.h'
      include 'def_file_list.h'
      include 'def_file_array.h'
      include 'def_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     ntest = 1000

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      integer, intent(in) ::
     &     idx_opin, idx_opout

      type(filinf), pointer ::
     &     ffin, ffout, ffv
      type(operator), pointer ::
     &     op_in, op_out, opv
      
      logical ::
     &     bufin, bufout
      integer ::
     &     len_str, idum, ifree, lblk, nblkmax,
     &     nblk, nbuff, ioffin, ioffout, idxst, idxnd, njoined,
     &     join_off, idoffin, idoffout, idxmsa,
     &     msmax, msa, msc, igama, igamc, idx, jdx, ngam,
     &     len_gam_ms, ioff, len_blk_in, ioff_blk_in, len_blk_out,
     &     ioff_blk_out,idxv, nocc_cls, iocc_cls
      integer ::
     &     opin_temp(ngastp,2), opout_temp(ngastp,2)
      
      real(8) ::
     &     fac, cpu, sys, wall, cpu0, sys0, wall0
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     iocc_equal, irestr_equal

      call atim_csw(cpu0,sys0,wall0)

      ! Locate the file containing the necessary elements of the input
      ! operator.
      op_in => op_info%op_arr(idx_opin)%op
      ffin => op_info%opfil_arr(idx_opin)%fhand
      if(.not.associated(ffin))
     &     call quit(1,'symm_op','no file handle for'//
     &     trim(op_in%name))

      ! Assign a file to the output operator.
      op_out => op_info%op_arr(idx_opout)%op
      call assign_file_to_op(idx_opout,.true.,ffout,
     &                       1,1,1,
     &                       0,op_info)
      ffout => op_info%opfil_arr(idx_opout)%fhand

      ! Assumptions for the B-matrix in MP2-R12:
      fac = 0.5d0
      nocc_cls = op_in%n_occ_cls
      if(nocc_cls.ne.op_out%n_occ_cls)
     &     call quit(1,'symm_op','no. of occupation classes unequal')

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' operator symmetrisation: '
        write(luout,*) '=========================='
        write(luout,*) ' fac = ',fac
        write(luout,*) ' ffin: ',trim(ffin%name),
     &                   ' rec: ',ffin%current_record
        write(luout,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(luout,*) ' opinp: ',op_in%name(1:len_trim(op_in%name))
        write(luout,*) ' opinv: ',op_out%name(1:len_trim(op_out%name))
      end if

      njoined = op_in%njoined
      if(njoined.ne.op_out%njoined)
     &     call quit(1,'symm_op','in and out incompatible: njoined')

      ! Check to see whether the two operators have the same shapes in
      ! each block. Assumes that the equivalent blocks are ordered the 
      ! same in each operator.
      opin_temp(1:ngastp,1:2)=0
      opout_temp(1:ngastp,1:2)=0

      do iocc_cls = 1, nocc_cls
        join_off = (iocc_cls-1)*njoined

        do idx=1,njoined

          opin_temp(1:ngastp,1:2) =
     &         op_in%ihpvca_occ(1:ngastp,1:2,join_off+idx)
          opout_temp(1:ngastp,1:2) = 
     &         op_out%ihpvca_occ(1:ngastp,1:2,join_off+idx)

          if (.not.iocc_equal(opin_temp,op_in%dagger,
     &         opout_temp,op_out%dagger)) then
            call quit(1,'symm_op','in and out incompatible: occs.')
          endif  
        
        enddo
      enddo

      ! Call the routine which actually does the symmetrisation.
      call symmetrise(fac,ffin,op_in,ffout,op_out,nocc_cls,
     &     op_info,orb_info)


      if(ntest.ge.1000)then
        write(luout,*)'Symmetrised operator: ',trim(op_out%name)
        call wrt_op_file(luout,5,ffout,op_out,1,
     &       op_out%n_occ_cls,str_info,orb_info)
      endif

      call atim_csw(cpu,sys,wall)

      call prtim(luout,'time for symmetrisation',
     &     cpu-cpu0,sys-sys0,wall-wall0)

c dbg
c      ! Evaluate the inverse of B and multiply it by V+. Used to trick the
c      ! preconditioner of MP2-R12 1A.
c      call file_open(ffout)
c      call invert(ffout,op_out,ffout,op_out,1,
c     &     op_info,orb_info)

c      if(ntest.ge.1000)then
c        write(luout,*)'Inverted matrix: '
c        call wrt_op_file(luout,5,ffout,op_out,1,
c     &       op_out%n_occ_cls,str_info,orb_info)
c      endif
      
c      idxv = idx_oplist2(op_vbar_inter,op_info)
c      opv => op_info%op_arr(idxv)%op
c      ffv => op_info%opfil_arr(idxv)%fhand
c      if(.not.associated(ffv))
c     &     call quit(1,'symm_op','no file handle for'//
c     &     trim(opv%name))

c      call op_mult(-1d0,ffout,op_out,1,ffv,opv,1,
c     &     orb_info)

      if(ffout%unit.gt.0)then
        call file_close_keep(ffout)
      endif

c      if(ntest.ge.1000)then
c        write(luout,*)'Modified V+'
c        call wrt_op_file(luout,5,ffv,opv,1,
c     &       opv%n_occ_cls,str_info,orb_info)
c      endif
      
c      stop
c dbg

      return
      end
