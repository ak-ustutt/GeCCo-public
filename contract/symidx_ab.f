*----------------------------------------------------------------------*
      subroutine symidx_ab(idxlist_out,
     &     idxlist_in,nlist,mel,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     given a list of indices (within bounds of ME-list mel)
*     return the corresponding indices if AB-symmetrization is
*     carried out incl. sign due to phase changes
*
*     nov 2008, andreas
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'mdef_operator_info.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(operator_info), intent(in) ::
     &     op_info
      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     nlist, idxlist_in(nlist)
      integer, intent(out) ::
     &     idxlist_out(nlist)

      integer ::
     &     nocc_cls, njoined, ioff_blk,
     &     ifree, nblk, nbuff, iblk,
     &     idx, ngam, ioff, lenblk
      real(8) ::
     &     fac_rel, fac_pre

      type(operator), pointer ::
     &     op

      real(8), external ::
     &     ddot
      logical, external ::
     &     list_in_bounds

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'symidx_ab')
        write(lulog,*) 'MEL:  ',trim(mel%label)
        write(lulog,*) 'idxlist(IN):'
        write(lulog,'(1x,5i6,x,5i6)') idxlist_in(1:nlist)
      end if

      idxlist_out(1:nlist) = 0

      op  => mel%op

      nocc_cls = op%n_occ_cls
      njoined  = op%njoined
 
      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Loop over occupation classes.
      iocc_loop: do iblk = 1, nocc_cls 

        if(op%formal_blk(iblk)) cycle

        ioff_blk = (iblk-1)*njoined

        ioff = mel%off_op_occ(iblk)
        lenblk = mel%len_op_occ(iblk)

        if (.not.list_in_bounds(idxlist_in,nlist,ioff+1,ioff+lenblk))
     &       cycle

        if (ntest.ge.100)
     &       write(lulog,*) 'iblk, ioff, lenblk: ',iblk, ioff, lenblk

        call symidx_ab_blk(idxlist_out,
     &                  idxlist_in, nlist,
     &                  mel,iblk,
     &                  str_info,strmap_info,orb_info)

      end do iocc_loop

      if (ntest.ge.100) then
        write(lulog,*) 'idxlist(OUT):'
        write(lulog,'(1x,5i6,x,5i6)') idxlist_out(1:nlist)
      end if

      return
      end
