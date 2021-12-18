      subroutine make_ext2typ_reo(orb_info)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'

      character(len=16) :: i_am = 'make_ext2typ_reo'
      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     ngam, ntoob, nact, igam, ii, jj, idx_gam(8)
      integer, pointer ::
     &     ext_gamorb(:), int_gamorb(:), ireost(:)


      if (ntest.ge.100) call write_title(lulog,wst_dbg_subr,i_am)

      ngam = orb_info%nsym
      ntoob = orb_info%ntoob
      nact  = ntoob-orb_info%ncore_ext
      ext_gamorb => orb_info%ext_gamorb
      int_gamorb => orb_info%igamorb
      ireost => orb_info%ireost

      if (ntest.ge.100) then
        write(lulog,'(x,"external ordering:")')
        write(lulog,'(x,5i4,x,5i4)') ext_gamorb(1:ntoob)
        write(lulog,'(x,"internal ordering:")')
        write(lulog,'(x,5i4,x,5i4)') int_gamorb(1:ntoob)
        if (orb_info%ncore_ext.gt.0) then
          write(lulog,'(x,"reordering fc -> ae")')
          write(lulog,'(x,5i4,x,5i4)') orb_info%ext_fcreo(1:nact)
        end if
      end if

      ! run over the ext_gamorb array
      idx_gam(1:ngam) = 0
      do ii = 1, ntoob
        igam = ext_gamorb(ii) ! current orbital
        if (ntest.ge.100) write(lulog,*) ' > ',ii,igam
        jj = idx_gam(igam)+1
        search: do  ! find matching lowest orbital of same symmetry
          if (int_gamorb(jj).eq.igam.or.jj.gt.ntoob) exit search
          jj = jj+1
        end do search
        if (ntest.ge.100) write(lulog,*) '   > ',jj 
        if (jj.gt.ntoob) call quit(1,i_am,'something is strange')
        ireost(ii) = jj  ! set entry of reordering array
        idx_gam(igam) = jj
      end do

      end
