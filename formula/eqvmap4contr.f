*----------------------------------------------------------------------*
      subroutine eqvfac4contr(ieqvfac,neqv_tup,idx_tup,ntup,
     &     topomap,eqv_map,neqv,idx_eqv,contr)
*----------------------------------------------------------------------*
*     using the equivalence maps, find out which tupels of vertices
*     are really equivalent
*
*     topomap,eqv_map,neqv,idx_eqv are provided by topomap4contr
*
*     ntup is the number of tupels of equivalent vertices
*     neqv_tup(ntup) is the number of equivalent vertices per tupel
*     idx_tup(ntup) is the index of the first vertex in the tupel
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     ieqvfac, ntup, neqv_tup(*), idx_tup(*)
      type(contraction), intent(in) ::
     &     contr
c      logical, intent(in) ::
c     &     fix_vtx(contr%nvtx)
      integer, intent(in) ::
     &     topomap(contr%nvtx,contr%nvtx),eqv_map(contr%nvtx),
     &     neqv(contr%nvtx),idx_eqv(contr%nvtx,contr%nvtx)

      integer ::
     &     nvtx, ivtx, kvtx, lvtx, idx, itopo, itup

      integer, external ::
     &     ifac

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'eqvmap4contr at work')

      nvtx = contr%nvtx
      ntup = 0
      ivtx_loop: do ivtx = 1, nvtx
        if (neqv(ivtx).lt.0) cycle
        ntup = ntup + 1
        idx_tup(ntup) = ivtx
        neqv_tup(ntup) = 1
        kvtx = ivtx
        do idx = 2, neqv(ivtx)
          lvtx = idx_eqv(idx,ivtx)
          itopo = topo_cmp(topomap(1:nvtx,kvtx),topomap(1:nvtx,lvtx),
     &         eqv_map,contr%arc,nvtx)
          kvtx = lvtx
c dbg
c          print *,'comparing'
c          print *,'1: ',topomap(1:nvtx,ivtx)
c          print *,'2: ',topomap(1:nvtx,kvtx)
c          print *,'result -> ',itopo
c dbg
          if (itopo.lt.0) then
            ntup = ntup+1
            idx_tup(ntup) = lvtx
            neqv_tup(ntup) = 1
          else if (itopo.eq.0) then
            neqv_tup(ntup) = neqv_tup(ntup) + 1
          else if (itopo.gt.0) then
            ! put a minus sign to indicate that this
            ! contraction is not in standard order
            neqv_tup(ntup) = -abs(neqv_tup(ntup))
            ntup = ntup+1
            idx_tup(ntup) = lvtx
            neqv_tup(ntup) = 1
          end if
        end do
      end do ivtx_loop

      if (ntest.ge.100) then
        write(luout,*) 'ntup = ',ntup
        write(luout,*) 'neqv_tup: ',neqv_tup(1:ntup)
        write(luout,*) 'idx_tup:  ',idx_tup(1:ntup)
      end if

      ieqvfac = 1
      do itup = 1, ntup
        ieqvfac = ieqvfac*ifac(abs(neqv_tup(itup)))
        if (neqv_tup(itup).lt.0) ieqvfac = -abs(ieqvfac)
      end do

      if (ntest.ge.100) then
        write(luout,*) '-> fac = ',ieqvfac
      end if

      return

      contains

      logical function topo_cmp(top1,top2,eqv,arc,nel)

      implicit none

      integer, intent(in) ::
     &     nel, top1(nel), top2(nel), eqv(nel)
      type(cntr_arc), intent(in) ::
     &     arc(*)

      integer ::
     &     iel, jel, ieqv, npick, nconn1, nconn2
      integer ::
     &     pick1(nel), pick2(nel), eqv_cp(nel)

      integer, external ::
     &     imltlist

      ! number of vertices connected to
      nconn1 = nvtx-imltlist(0,top1,nvtx,1)
      nconn2 = nvtx-imltlist(0,top2,nvtx,1)

      if (nconn1.gt.nconn2) then
        topo_cmp = -1
        return
      else if (nconn1.lt.nconn2) then
        topo_cmp = +1
        return
      end if

      eqv_cp = eqv ! get a copy
      do iel = 1, nel
        ! already processed this one?
        if (eqv_cp(iel).le.0) cycle
        ! current eqv number is:
        ieqv = eqv_cp(iel)
        ! pick all elements from eqv places:
        npick = 1
        pick1(1) = top1(iel)
        pick2(1) = top2(iel)
        do jel = iel+1, nel
          if (eqv_cp(jel).eq.ieqv) then
            eqv_cp(jel) = -1 ! mark this element as "already processed"
            npick = npick+1
            pick1(npick) = top1(jel)
            pick2(npick) = top2(jel)
          end if
        end do
        ! compare
        if (npick.eq.1) then
          if (pick1(1).lt.pick2(1)) then
            topo_cmp = -1
          else if (pick1(1).eq.pick2(1)) then
            topo_cmp = 0
          else
            topo_cmp = +1
          end if
        else
          ! sort the two lists
          call isort(pick1,npick,+1)
          call isort(pick2,npick,+1)
          topo_cmp = 0
          do jel = 1, npick
            if (pick1(jel).lt.pick2(jel)) then
              topo_cmp = -1
              exit
            else if (pick1(jel).gt.pick2(jel)) then
              topo_cmp = +1
              exit
            end if
          end do
        end if
        ! if different, we can leave the loop
        if (topo_cmp.ne.0) exit
      end do

      return
      end function
      
      end
