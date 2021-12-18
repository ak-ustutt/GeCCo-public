*----------------------------------------------------------------------*
      subroutine put_orbinfo(orb_info, fforbinf)
*----------------------------------------------------------------------*
*     Create a formated file with the content of orb_inf
*     The file handle is returned via fforbinf
*
*     The purpose is create a file that will be read by a Python target
*     script.
*
*     yuri, oct 2014
*----------------------------------------------------------------------*
      implicit none

      include 'def_orbinf.h'
      include 'opdim.h'
      include 'def_filinf.h'

      type(orbinf), intent(in) ::
     &     orb_info
      type(filinf), intent(out) ::
     &     fforbinf

      character(20), parameter ::
     &     name_orbinf = "orb_info.gecco"
      integer ::
     &     idum, luorb, i, j, k

      call file_init(fforbinf,trim(name_orbinf),ftyp_sq_frm,idum)
      call file_open(fforbinf)
      luorb = fforbinf%unit

      write(luorb,fmt='(a,x,i0)') "ngastp", ngastp
      write(luorb,fmt='(a,x,i0)') "two", 2
      write(luorb,fmt='(a,x,i0)') "ntoob_caborb",
     &     orb_info%ntoob+orb_info%caborb


      write(luorb,fmt='(a,x,i0)') "nsym",  orb_info%nsym
      write(luorb,fmt='(a,x,i0)') "ngas",  orb_info%ngas
      write(luorb,fmt='(a,x,i0)') "nspin", orb_info%nspin


      write(luorb,fmt='(a,x,i0)') "ntoob",  orb_info%ntoob
      write(luorb,fmt='(a,x,i0)') "caborb", orb_info%caborb
      write(luorb,fmt='(a,x,i0)') "nbast",  orb_info%nbast
      write(luorb,fmt='(a,x,i0)') "nxbast", orb_info%nxbast
      write(luorb,fmt='(a,x,i0)') "nactel", orb_info%nactel
      write(luorb,fmt='(a,x,i0)') "nactorb",orb_info%nactorb
      write(luorb,fmt='(a,x,i0)') "lsym",   orb_info%lsym
      write(luorb,fmt='(a,x,i0)') "imult",  orb_info%imult
      write(luorb,fmt='(a,x,i0)') "ims",    orb_info%ims


      write(luorb,fmt='(a)', advance='no') "igassh nsym ngas"
      do j = 1,orb_info%ngas
       do i = 1,orb_info%nsym
        write (luorb,fmt='(x,i0)', advance='no') orb_info%igassh(i,j)
       end do
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "nbas nsym"
      do i = 1,orb_info%nsym
       write (luorb,fmt='(x,i0)', advance='no') orb_info%nbas(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "ntoobs nsym"
      do i = 1,orb_info%nsym
       write (luorb,fmt='(x,i0)', advance='no') orb_info%ntoobs(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "ireots ntoob_caborb"
      do i = 1,orb_info%ntoob+orb_info%caborb
       write (luorb,fmt='(x,i0)', advance='no') orb_info%ireots(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "ireost ntoob_caborb"
      do i = 1,orb_info%ntoob+orb_info%caborb
       write (luorb,fmt='(x,i0)', advance='no') orb_info%ireost(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "igamorb ntoob_caborb"
      do i = 1,orb_info%ntoob+orb_info%caborb
       write (luorb,fmt='(x,i0)', advance='no') orb_info%igamorb(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "igasorb ntoob_caborb"
      do i = 1,orb_info%ntoob+orb_info%caborb
       write (luorb,fmt='(x,i0)', advance='no') orb_info%igasorb(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "mostnd two nsym ngas"
      do k = 1,orb_info%ngas
       do j = 1,orb_info%nsym
        do i = 1,2
         write (luorb,fmt='(x,i0)', advance='no') orb_info%mostnd(i,j,k)
        end do
       end do
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "iad_gas ngas"
      do i = 1,orb_info%ngas
       write (luorb,fmt='(x,i0)', advance='no') orb_info%iad_gas(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "gas_reo ngas"
      do i = 1,orb_info%ngas
       write (luorb,fmt='(x,i0)', advance='no') orb_info%gas_reo(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "ihpvgas ngas nspin"
      do j = 1,orb_info%nspin
       do i = 1,orb_info%ngas
        write (luorb,fmt='(x,i0)', advance='no') orb_info%ihpvgas(i,j)
       end do
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "ngas_hpv ngastp"
      do i = 1,ngastp
       write (luorb,fmt='(x,i0)', advance='no') orb_info%ngas_hpv(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "nactt_hpv ngastp"
      do i = 1,ngastp
       write (luorb,fmt='(x,i0)', advance='no') orb_info%nactt_hpv(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "norb_hpv ngastp nspin"
      do j = 1,orb_info%nspin
       do i = 1,ngastp
        write (luorb,fmt='(x,i0)', advance='no') orb_info%norb_hpv(i,j)
       end do
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "idx_gas ngastp"
      do i = 1,ngastp
       write (luorb,fmt='(x,i0)', advance='no') orb_info%idx_gas(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "ioff_gas ngastp"
      do i = 1,ngastp
       write (luorb,fmt='(x,i0)', advance='no') orb_info%ioff_gas(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "cab_orb nsym"
      do i = 1,orb_info%nsym
       write (luorb,fmt='(x,i0)', advance='no') orb_info%cab_orb(i)
      end do
      write (luorb,*)

      write(luorb,fmt='(a)', advance='no') "nxbas nsym"
      do i = 1,orb_info%nsym
       write (luorb,fmt='(x,i0)', advance='no') orb_info%nxbas(i)
      end do
      write (luorb,*)



      write(luorb,fmt='(a,x,i0)') "n_bound_orbs", orb_info%n_bound_orbs
      write(luorb,fmt='(a,x,i0)') "n_freeze_rcmd",orb_info%n_freeze_rcmd


      write(luorb,fmt='(a)', advance='no')"isym_bound_orbs n_bound_orbs"
      do i = 1,orb_info%n_bound_orbs
       write (luorb,fmt='(x,i0)', advance='no')
     &      orb_info%isym_bound_orbs(i)
      end do
      write (luorb,*)


      call file_close_keep(fforbinf)
      end
