*----------------------------------------------------------------------*
      subroutine init_operator(ipass,op,orb_info)
*----------------------------------------------------------------------*
*     allocate operator sub-arrays
*      on entry, at least op%n_occ_cls must be set
*      ipass==0:   set primary arrays
*      ipass==1:   set arrays needed in first pass of set_op_dim
*      ipass==2:   set arrays needed in second pass of set_op_dim
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'

      integer, intent(in) ::
     &     ipass
      type(operator), intent(inout) ::
     &     op
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     nblk, iblk, nsym, ndis, nexc

      if (op%n_occ_cls.le.0.or.op%n_occ_cls.ge.1000) then
        write(luout,*) 'n_occ_cls = ',op%n_occ_cls
        call quit(1,'init_operator',
     &              'suspicious number of blocks (bug?)')
      end if

      nblk = op%n_occ_cls
      select case(ipass)
      case(0)
        allocate(op%ihpvca_occ(ngastp,2,nblk),
     &           op%ica_occ(2,nblk),
     &           op%igasca_restr(2,orb_info%ngas,2,2,nblk),
     &           op%len_op_occ(nblk),
     &           op%idx_graph(ngastp,2,nblk))        
        
        op%off_op_occ => null()
        op%len_op_occ => null()
        op%off_op_gmo => null()
        op%len_op_gmo => null()
        op%off_op_gmox => null()
        op%len_op_gmox => null()

      case(1)
        nsym = orb_info%nsym
        allocate(op%off_op_occ(nblk),
     &           op%len_op_occ(nblk),
     &           op%off_op_gmo(nblk),
     &           op%len_op_gmo(nblk),
     &           op%off_op_gmox(nblk),
     &           op%len_op_gmox(nblk))
        do iblk = 1, nblk
          nexc = min(op%ica_occ(1,iblk),
     &               op%ica_occ(2,iblk))
          allocate(op%len_op_gmo(iblk)%gam_ms(nsym,nexc+1),
     &             op%off_op_gmo(iblk)%gam_ms(nsym,nexc+1))
        end do        
      case(2)
        do iblk = 1, nblk
          nexc = min(op%ica_occ(1,iblk),
     &               op%ica_occ(2,iblk))
          ndis = op%off_op_gmox(iblk)%maxd
          allocate(op%len_op_gmox(iblk)%
     &                d_gam_ms(ndis,nsym,nexc+1),
     &             op%off_op_gmox(iblk)%
     &                d_gam_ms(ndis,nsym,nexc+1),
     &             op%off_op_gmox(iblk)%
     &                did(ndis,nsym,nexc+1),
     &             op%off_op_gmox(iblk)%ndis(nsym,nexc+1))
        end do
      case default
        write(luout,*) 'ipass = ',ipass
        call quit(1,'init_operator',
     &              'illegal ipass-flag')
      end select

      return
      end
