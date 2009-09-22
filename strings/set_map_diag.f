*----------------------------------------------------------------------*
      subroutine set_map_diag(map,len,ms,igam,mel,iblk,
     &                        str_info,orb_info)
*----------------------------------------------------------------------*
*     set map for strings of diagonal elements of diagonal op. block
*
*     matthias 2009
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'ifc_baserout.h'

      integer, intent(in) ::
     &     ms, igam, iblk, len
      integer, intent(out) ::
     &     map(len)
      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     idxorb(2), idxspn(2), idxdss(2), igrph, ihpv,
     &     nsym, ifree, istr
      integer, pointer ::
     &     igas_restr(:,:,:,:,:),
     &     mostnd(:,:,:), idx_gas(:), ngas_hpv(:), igamorb(:)

      integer, pointer ::
     &     buffer(:)
      type(operator), pointer ::
     &     op
      logical ::
     &     first
      logical, external ::
     &     next_string

      if (ntest.ge.100) then
        write(luout,*) '=================='
        write(luout,*) 'this is set_map_diag'
        write(luout,*) '=================='
      end if

      ifree = mem_setmark('set_map_diag')

      op => mel%op
      igas_restr => str_info%igas_restr
      nsym = orb_info%nsym
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      idx_gas => orb_info%idx_gas
      ngas_hpv => orb_info%ngas_hpv


      if (.not.all((op%ihpvca_occ(1:ngastp,1,iblk)-
     &             op%ihpvca_occ(1:ngastp,2,iblk)).eq.0)
     &    .or.sum(op%ihpvca_occ(1:ngastp,1,iblk)).ne.2)
     &   call quit(0,'set_map_diag','only diagonal operator blocks!')
      do ihpv = 1, ngastp
        if (op%ihpvca_occ(ihpv,1,iblk).eq.0) then
          cycle
        else if (op%ihpvca_occ(ihpv,1,iblk).ne.2) then
          call quit(0,'set_map_diag','only diagonal operator blocks!')
        else
          igrph = mel%idx_graph(ihpv,1,iblk)
          exit
        end if
      end do

      first = .true.
      istr = 0

      str_loop: do
        if (.not.next_string(idxorb,idxspn,idxdss,
     &                 2,ms,igam,first,
     &                 igas_restr(1,1,1,1,igrph),
     &                 mostnd(1,1,idx_gas(ihpv)),igamorb,
     &                 nsym,ngas_hpv(ihpv))
     &                 ) exit str_loop

        first = .false.

        istr = istr+1
        map(istr) = (idxorb(1)-1)*orb_info%ntoob+idxorb(2)
        if (ntest.ge.100)
     &       write(luout,'(x,a,i4,x,2i4,x,i4)') 'istr, idxorb, map: ',
     &             istr,idxorb(1:2),map(istr)

      end do str_loop
      if (istr.ne.len) call quit(0,'set_map_diag',
     &                 'inconsistent map length')

      ifree = mem_flushmark()

      return
      end
