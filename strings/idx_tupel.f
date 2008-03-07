*----------------------------------------------------------------------*
      subroutine idx_tupel(idxhpv,
     &     iocc,idx_graph,idspc,idorb,idspn,idgam,set_spc_gam,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     given a C or A string, obtain the string numbers for 
*     H P V (seperately in array idxhpv, so outside routines
*     decides on actual storage order)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, intent(out) ::
     &     idxhpv(ngastp)
      logical, intent(in) ::
     &     set_spc_gam
      integer, intent(in) ::
     &     iocc(ngastp), idx_graph(ngastp),
     &     idorb(*), idspn(*)
      integer, intent(inout) ::
     &     idspc(*), idgam(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     nelt, nel, igraph, ipos, ihpv, ihpvdx, igtp

      type(graph), pointer ::
     &     curgraph

      logical, external ::
     &     allow_sbsp_dis
      integer, external ::
     &     idx4sg


      idxhpv(1:ngastp) = 0

      nelt = 0
      do ihpv = 1, ngastp
        nelt = nelt + iocc(ihpv)
      end do

      if (nelt.eq.0) return

      if (set_spc_gam) then
        do ipos = 1, nelt
          ! get IRREP of orbital
          idgam(ipos) = orb_info%igamorb(idorb(ipos))
          ! get subspace of orbital
          idspc(ipos) = orb_info%igasorb(idorb(ipos))
          ! correct counting: only inside H/P/V
          igtp = orb_info%ihpvgas(idspc(ipos),1)
          if (orb_info%nspin.eq.2.and.
     &        orb_info%ihpvgas(idspc(ipos),1).ne.
     &        orb_info%ihpvgas(idspc(ipos),2) .and.
     &        idspn(ipos).eq.-1)
     &         igtp = orb_info%ihpvgas(idspc(ipos),2)
          idspc(ipos) = idspc(ipos)-orb_info%idx_gas(igtp)+1
        end do
      end if

      ipos = 1
      ihpv_loop: do ihpvdx = 1, ngastp
        if (ipos.gt.nelt) exit ihpv_loop
        ! not quite efficient fix to get actual H/P/V
        ihpv = orb_info%ihpvgas(orb_info%igasorb(idorb(ipos)),1)
        if (orb_info%nspin.eq.2.and.
     &      orb_info%ihpvgas(idspc(ipos),1).ne.
     &      orb_info%ihpvgas(idspc(ipos),2) .and.
     &      idspn(ipos).eq.-1)
     &         ihpv = orb_info%ihpvgas(idspc(ipos),2)
        nel = iocc(ihpv)
        if (nel.eq.0) then
          call quit(1,'idx_tupel','strange event')
        end if

        ! point to graph needed for current string
        igraph = idx_graph(ihpv) !,ica,iblk_op)
        curgraph => str_info%g(igraph)

        ! check for restrictions
        if (.not.allow_sbsp_dis(idspc(ipos),nel,
     &                          orb_info%ngas_hpv(ihpv),
     &                        str_info%igas_restr(1,1,1,1,igraph))) then
                         !        ADAPT FOR OPEN-SHELL ^^^
          return
        end if

        ! get string index
        idxhpv(ihpv) = idx4sg(nel,idspc(ipos),idorb(ipos),
     &                               idspn(ipos),idgam(ipos),
     &             curgraph%y4sg,curgraph%yinf,
     &             curgraph%yssg,curgraph%wssg,
     &             orb_info%mostnd(1,1,orb_info%idx_gas(ihpv)),
     &             str_info%ispc_occ(igraph),orb_info%nsym,
     &             orb_info%ngas_hpv(ihpv))

        ! string number is actual index - 1

        ipos = ipos+nel
      end do ihpv_loop

      return
      end

