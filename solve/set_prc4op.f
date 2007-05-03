*----------------------------------------------------------------------*
      subroutine set_prc4op(idxprc,idxffprc,
     &                      idxop,idxham,idxffop,idxffham,
     &                      op_info,
     &                      str_info,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up a diagonal preconditioner
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, intent(in) ::
     &     idxprc, idxffprc, idxop, idxffop, idxham, idxffham
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ifree, nops
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      type(operator), pointer ::
     &     opop, opham, opprc
      type(filinf), pointer ::
     &     ffopham, ffopprc

      real(8), pointer ::
     &     x1dia(:)

      call atim(cpu0,sys0,wall0)

      nops = op_info%nops

      ifree = mem_setmark('prc4op')
      ! this assumption is probably not too bad:
      ifree = mem_alloc_real(x1dia,2*orb_info%ntoob,'x1dia')

      if (idxham.lt.0.or.idxham.gt.nops .or.
     &    idxprc.lt.0.or.idxprc.gt.nops .or.
     &    idxop.lt.0 .or.idxop.gt.nops) then
        write(luout,*) idxham, idxop, idxprc
        write(luout,*) 'allowed: 1, ... ,',nops
        call quit(1,'set_prc4op','buggy operator indices')
      end if

      opham => op_info%op_arr(idxham)%op
      opprc => op_info%op_arr(idxprc)%op
      opop  => op_info%op_arr(idxop)%op
      ffopham => op_info%opfil_arr(idxffham)%fhand
      ffopprc => op_info%opfil_arr(idxffprc)%fhand

      if (iprlvl.ge.1)
     &     write(luout,*) 'set up diagonal for ',trim(opop%name),
     &     ' from rank 1 part of ',trim(opham%name)

      ! extract the fock-matrix diagonal
      call onedia_from_op(x1dia,ffopham,opham,
     &     orb_info)

      ! set up preconditioner
      call dia4op(ffopprc,x1dia,opop,str_info,orb_info)      

      if (ntest.ge.1000) then
        call wrt_op_file(luout,5,ffopprc,opprc,1,
     &       opprc%n_occ_cls,str_info,orb_info)
      end if

      ifree = mem_flushmark()

      call atim(cpu,sys,wall)

      call prtim(luout,'time for diagonal ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
