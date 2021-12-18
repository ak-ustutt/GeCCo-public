*----------------------------------------------------------------------*
      subroutine set_integral_graph(iy_int,nints,nelmax,mnmxspc_in,
     &                              orb_info)
*----------------------------------------------------------------------*
*     set graph for addressing of raw integral lists
*     the routine generates a weight array for addressing the integral
*     via its (ordered) index n-tuple
*
*     andreas, march 2008
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     nelmax,
     &     mnmxspc_in(2,orb_info%ngas)
      integer, intent(out) ::
     &     iy_int(nelmax,orb_info%nsym,orb_info%ntoob+orb_info%caborb),
     &     nints(orb_info%nsym)

      integer ::
     &     norb, ngam, gamorb_cur, gasorb_cur,
     &     igam, jgam, iel, iorb, igas, nelmax_sp, nelmin_sp

      integer, pointer ::
     &     gamorb(:), gasorb(:), iw_int(:,:,:)

      integer ::
     &     mnmxspc(2,orb_info%ngas)


      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'set_integral_graph')
        write(lulog,*) 'nelmax: ',nelmax
      end if
      
      if (mnmxspc_in(1,1).eq.-1) then
        ! default: no restriction
        do igas = 1, orb_info%ngas
          mnmxspc(1:2,igas) = (/0,nelmax/)
        end do
      else
        mnmxspc = mnmxspc_in
      end if

      if (ntest.ge.100) then
        if (mnmxspc_in(1,1).eq.-1)
     &       write(lulog,*) 'using standard restrictions'
        do igas = 1, orb_info%ngas
          write(lulog,'(2x,i3,3x,2i3)') igas,mnmxspc(1:2,igas)
        end do
      end if

      norb = orb_info%ntoob+orb_info%caborb
      ngam = orb_info%nsym
      gamorb => orb_info%igamorb
      gasorb => orb_info%igasorb
      ! a) set graph for full four-fold ordered quadruples
      allocate(iw_int(0:nelmax,1:ngam,1:norb))
      iw_int(0:nelmax,1:ngam,1:norb) = 0

      !  i) weight array
      ! set head of graph
      iw_int(0,1,1) = 1
      ! first row
      gasorb_cur = gasorb(1)
      gamorb_cur = gamorb(1)
      nelmax_sp = mnmxspc(2,gasorb_cur)
      do iel = 1, nelmax_sp
        do igam = 1, ngam
          jgam = multd2h(igam,gamorb_cur)
          iw_int(iel,igam,1) =
     &         iw_int(iel,igam,1) + iw_int(iel-1,jgam,1)
        end do
      end do
            
      do iorb = 2, norb
        gamorb_cur = gamorb(iorb)        
        gasorb_cur = gasorb(iorb)
        nelmin_sp = 0
        if (gasorb_cur.gt.1) nelmin_sp = mnmxspc(1,gasorb_cur-1)
        nelmax_sp = mnmxspc(2,gasorb_cur)
        iw_int(nelmin_sp,1:ngam,iorb) = iw_int(nelmin_sp,1:ngam,iorb-1)
        do iel = nelmin_sp+1, nelmax_sp
          do igam = 1, ngam
            jgam = multd2h(igam,gamorb_cur)
            iw_int(iel,igam,iorb) =
     &           iw_int(iel,igam,iorb) + iw_int(iel-1,jgam,iorb) 
     &                                 + iw_int(iel  ,igam,iorb-1)
          end do
        end do        
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'generated vertex weights:'
        do iorb = 1, norb
          write(lulog,'(2x,i2,x,i2,5x,10i6)')
     &         iorb, 1, iw_int(0:nelmax,1,iorb)
          do igam = 2, ngam
            write(lulog,'(5x,i2,5x,10i6)')
     &           igam, iw_int(0:nelmax,igam,iorb)
          end do
        end do
      end if

      nints(1:ngam) = iw_int(nelmax,1:ngam,norb)

      !  ii) transform to arc weight array 
      iy_int = 0
      do iorb = 2, norb
        do igam = 1, ngam
          do iel = 1, nelmax
            iy_int(iel,igam,iorb) = iw_int(iel  ,igam,iorb-1)
          end do
        end do
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'generated arc weights:'
        do iorb = 1, norb
          write(lulog,'(2x,i2,x,i2,5x,10i6)')
     &         iorb, 1, iy_int(1:nelmax,1,iorb)
          do igam = 2, ngam
            write(lulog,'(5x,i2,5x,10i6)')
     &           igam, iy_int(1:nelmax,igam,iorb)
          end do
        end do
      end if

      deallocate(iw_int)

      return
      end
