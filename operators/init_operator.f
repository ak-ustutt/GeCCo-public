*----------------------------------------------------------------------*
      subroutine init_operator(op,orb_info)
*----------------------------------------------------------------------*
*     allocate operator sub-arrays
*      on entry, at least op%n_occ_cls must be set
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      type(operator), intent(inout) ::
     &     op
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     nblk, nblkt, iblk, nsym, ndis, nexc, ncount, ifree

      ! allocate in operator section:
      call mem_pushmark()
      ifree = mem_gotomark(operator_def)

      if (op%n_occ_cls.le.0.or.op%n_occ_cls.ge.1000) then
        write(luout,*) 'n_occ_cls = ',op%n_occ_cls
        call quit(1,'init_operator',
     &              'suspicious number of blocks (bug?)')
      end if

      nblk = op%n_occ_cls
      op%ngas = orb_info%ngas  ! remember dimension
      op%nspin = orb_info%nspin

        ! some arrays run over 1..njoined as second index
      nblkt = nblk * op%njoined
      ifree = mem_register(2*ngastp*nblkt
     &                      +2*nblk
     &                      +8*orb_info%ngas*nblkt*orb_info%nspin
     &                      +nblk,
     &       trim(op%name))
      allocate(op%ihpvca_occ(ngastp,2,nblkt),
     &           op%ica_occ(2,nblk),
     &           op%igasca_restr(2,orb_info%ngas,2,2,
     &                             orb_info%nspin,nblkt),
     &           op%formal_blk(nblk))

      call mem_popmark()

      return
      end
