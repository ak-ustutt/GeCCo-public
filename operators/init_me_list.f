*----------------------------------------------------------------------*
      subroutine init_me_list(ipass,mel,orb_info)
*----------------------------------------------------------------------*
*     allocate me_list sub-arrays
*      on entry, mel%op must be set
*      ipass==1:   set arrays needed in first pass of set_op_dim
*      ipass==2:   set arrays needed in second pass of set_op_dim
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      integer, intent(in) ::
     &     ipass
      type(me_list), intent(inout) ::
     &     mel
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     nblk, nblkt, iblk, nsym, ndis, nexc, ncount, ifree

      ! allocate in me_list section:
      call mem_pushmark()
      ifree = mem_gotomark(me_list_def)

      if (.not.associated(mel%op))
     &     call quit(1,'init_me_list',
     &     'no operator referenced by ME-list "'//trim(mel%label)//'"')

      if (mel%op%n_occ_cls.le.0.or.mel%op%n_occ_cls.ge.1000) then
        write(luout,*) 'n_occ_cls = ',mel%op%n_occ_cls
        call quit(1,'init_me_list',
     &              'suspicious number of blocks (bug?)')
      end if

      mel%frequency = 0d0      ! set assigned frequency to zero

      nblk = mel%op%n_occ_cls
      select case(ipass)
      case(1)
        ! some arrays run over 1..njoined as second index
        nblkt = nblk * mel%op%njoined
        nsym = orb_info%nsym
        mel%nsym = nsym    !  remember dimension
        allocate(mel%len_op_occ(nblk),
     &           mel%idx_graph(ngastp,2,nblkt),
     &           mel%off_op_occ(nblk),
     &           mel%len_op_occ(nblk),
     &           mel%off_op_gmo(nblk),
     &           mel%len_op_gmo(nblk),
     &           mel%off_op_gmox(nblk),
     &           mel%len_op_gmox(nblk),
     &           mel%ld_op_gmox(nblk))
        ncount = nblk+2*ngastp*nblkt
        do iblk = 1, nblk
c          if (mel%op%formal_blk(iblk)) cycle
          nexc = min(mel%op%ica_occ(1,iblk),
     &               mel%op%ica_occ(2,iblk))
          allocate(mel%len_op_gmo(iblk)%gam_ms(nsym,nexc+1),
     &             mel%off_op_gmo(iblk)%gam_ms(nsym,nexc+1))
          ncount = ncount+2*nsym*(nexc+1)
        end do        
        ifree = mem_register(6*nblk+ncount,
     &       trim(mel%label)//'-1')
      case(2)
        nsym = orb_info%nsym
        if (nsym.ne.mel%nsym)
     &       call quit(1,'init_me_list',
     &       'setting orb_info%nsym does not conform with '//
     &       'nsym associated with current list')
        ncount = 0
        do iblk = 1, nblk
c          if (mel%op%formal_blk(iblk)) cycle
          nexc = min(mel%op%ica_occ(1,iblk),
     &               mel%op%ica_occ(2,iblk))
          ndis = mel%off_op_gmox(iblk)%maxd
          allocate(mel%len_op_gmox(iblk)%
     &                d_gam_ms(ndis,nsym,nexc+1),
     &             mel%ld_op_gmox(iblk)%
     &                d_gam_ms(ndis,nsym,nexc+1),
     &             mel%off_op_gmox(iblk)%
     &                d_gam_ms(ndis,nsym,nexc+1),
     &             mel%off_op_gmox(iblk)%
     &                did(ndis,nsym,nexc+1),
     &             mel%off_op_gmox(iblk)%ndis(nsym,nexc+1))
          ncount = ncount+ndis*nsym*(nexc+1)*3+nsym*(nexc+1)
        end do
        ifree = mem_register(ncount,trim(mel%label)//'-2')
      case default
        write(luout,*) 'ipass = ',ipass
        call quit(1,'init_me_list',
     &              'illegal ipass-flag')
      end select

      call mem_popmark()

      return
      end
