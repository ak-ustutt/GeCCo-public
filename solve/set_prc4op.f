*----------------------------------------------------------------------*
      subroutine set_prc4op(idxprc,idxffprc,
     &                      idxop,idxham,idxffop,idxffham,
     &                      ffops,ops,nops,
     &                      str_info,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up a diagonal preconditioner
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, intent(in) ::
     &     idxprc, idxffprc, idxop, idxffop, idxham, idxffham, nops
      type(operator), intent(in) ::
     &     ops(nops)
      type(filinf), intent(in) ::
     &     ffops(nops)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ifree
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      real(8), pointer ::
     &     x1dia(:)

      call atim(cpu0,sys0,wall0)

      ifree = mem_setmark('prc4op')
      ! this assumption is probably not too bad:
      ifree = mem_alloc_real(x1dia,2*orb_info%ntoob,'x1dia')

      if (iprlvl.ge.1)
     &     write(luout,*) 'set up diagonal for ',trim(ops(idxop)%name),
     &     ' from rank 1 part of ',trim(ops(idxham)%name)


      if (idxham.lt.0.or.idxham.gt.nops .or.
     &    idxprc.lt.0.or.idxprc.gt.nops .or.
     &    idxop.lt.0 .or.idxop.gt.nops) then
        write(luout,*) idxham, idxop, idxprc
        write(luout,*) 'allowed: 1, ... ,',nops
        call quit(1,'set_prc4op','buggy operator indices')
      end if

      ! extract the fock-matrix diagonal
      call onedia_from_op(x1dia,ffops(idxffham),ops(idxham),
     &     orb_info)

      ! set up preconditioner
      call dia4op(ffops(idxffprc),x1dia,ops(idxop),str_info,orb_info)      

      if (ntest.ge.1000) then
        call wrt_op_file(luout,5,ffops(idxffprc),ops(idxprc),1,
     &       ops(idxprc)%n_occ_cls,str_info,orb_info)
      end if

      ifree = mem_flushmark()

      call atim(cpu,sys,wall)

      call prtim(luout,'time for diagonal ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
