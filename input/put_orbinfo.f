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
     &     idum, luorb

      call file_init(fforbinf,trim(name_orbinf),ftyp_sq_frm,idum)
      call file_open(fforbinf)
      luorb = fforbinf%unit

      write(luorb,*) "ngastp ", ngastp
      write(luorb,*) "two ", 2
      write(luorb,*) "ntoob_caborb ",orb_info%ntoob+orb_info%caborb

      write(luorb,*) "nsym ",  orb_info%nsym
      write(luorb,*) "ngas ",  orb_info%ngas
      write(luorb,*) "nspin ", orb_info%nspin

      write(luorb,*) "ntoob ",  orb_info%ntoob
      write(luorb,*) "caborb ", orb_info%caborb
      write(luorb,*) "nbast ",  orb_info%nbast
      write(luorb,*) "nxbast ", orb_info%nxbast
      write(luorb,*) "nactel ", orb_info%nactel
      write(luorb,*) "nactorb ",orb_info%nactorb
      write(luorb,*) "lsym ",   orb_info%lsym
      write(luorb,*) "imult ",  orb_info%imult
      write(luorb,*) "ims ",    orb_info%ims

      write(luorb,*) "igassh nsym ngas ",
     &       orb_info%igassh(1:orb_info%nsym,1:orb_info%ngas)
      write(luorb,*) "nbas nsym ",
     &       orb_info%nbas(1:orb_info%nsym)
      write(luorb,*) "ntoobs nsym ",
     &       orb_info%ntoobs(1:orb_info%nsym)
      write(luorb,*) "ireots ntoob_caborb ",
     &       orb_info%ireots(1:orb_info%ntoob+orb_info%caborb)
      write(luorb,*) "ireost ntoob_caborb ",
     &       orb_info%ireost(1:orb_info%ntoob+orb_info%caborb)
      write(luorb,*) "igamorb ntoob_caborb ",
     &       orb_info%igamorb(1:orb_info%ntoob+orb_info%caborb)
      write(luorb,*) "igasorb ntoob_caborb ",
     &       orb_info%igasorb(1:orb_info%ntoob+orb_info%caborb)
      write(luorb,*) "mostnd two nsym ngas ",
     &       orb_info%mostnd(1:2,1:orb_info%nsym,1:orb_info%ngas)
      write(luorb,*) "iad_gas ngas ",
     &       orb_info%iad_gas(1:orb_info%ngas)
      write(luorb,*) "gas_reo ngas ",
     &       orb_info%gas_reo(1:orb_info%ngas)
      write(luorb,*) "ihpvgas ngas nspin ",
     &       orb_info%ihpvgas(1:orb_info%ngas,1:orb_info%nspin)
      write(luorb,*) "ngas_hpv ngastp ",
     &       orb_info%ngas_hpv (1:ngastp)
      write(luorb,*) "nactt_hpv ngastp ",
     &       orb_info%nactt_hpv(1:ngastp)
      write(luorb,*) "norb_hpv ngastp nspin ",
     &       orb_info%norb_hpv(1:ngastp,1:orb_info%nspin)
      write(luorb,*) "idx_gas ngastp ",
     &       orb_info%idx_gas (1:ngastp)
      write(luorb,*) "ioff_gas ngastp ",
     &       orb_info%ioff_gas(1:ngastp)
      write(luorb,*) "cab_orb nsym ",
     &       orb_info%cab_orb(1:orb_info%nsym)
      write(luorb,*) "nxbas nsym ",
     &       orb_info%nxbas(1:orb_info%nsym)

      write(luorb,*) "n_bound_orbs ", orb_info%n_bound_orbs
      write(luorb,*) "n_freeze_rcmd ",orb_info%n_freeze_rcmd

      write(luorb,*) "isym_bound_orbs n_bound_orbs",
     &     orb_info%isym_bound_orbs(1:orb_info%n_bound_orbs)

      call file_close_keep(fforbinf)
      end
