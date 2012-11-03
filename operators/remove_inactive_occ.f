      subroutine remove_inactive_occ(occ,ndef,freeze,njoined,orb_info)
*
*     occ   input: defined occ.s    
*           output: dto. with unneccessary ones rem.
*     ndef  input: number of occ.s
*           output: dto. minus removed ones
*
*     freeze: number of electrons allowed in frozen shell
*             (for F12 and gradient stuff)
*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'

      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     njoined,
     &     freeze(2,njoined)
      integer, intent(inout) ::
     &     ndef,
     &     occ(ngastp,2,njoined,ndef)

      integer ::
     &     ihpvx, ijoin, ica, idef, jdef, ndel
      logical ::
     &     active(2,ngastp,njoined)
      logical, pointer ::
     &     delete(:)

      if (orb_info%nspin.ne.1) 
     &    call quit(1,'remove_inactive_occ','not adapted for nspin<>1')
  
      do ijoin = 1, njoined
        do ihpvx = 1, ngastp
          if (orb_info%norb_hpv(ihpvx,1).gt.0) then ! ADAPT FOR ISPIN!
            if (orb_info%nactt_hpv(ihpvx).eq.0) then
              active(1,ihpvx,ijoin) = freeze(1,ijoin).gt.0
              active(2,ihpvx,ijoin) = freeze(2,ijoin).gt.0
            else
              active(1:2,ihpvx,ijoin)=.true.
            end if
          else
            active(1:2,ihpvx,ijoin)=.false.
          end if
        end do
      end do

      allocate(delete(ndef))
      delete(1:ndef) = .false.
      ndel = 0
      do idef = 1, ndef
        do ijoin = 1, njoined
          do ihpvx = 1, ngastp
            do ica = 1, 2
              if (occ(ihpvx,ica,ijoin,idef).eq.0) cycle
              delete(idef)=delete(idef).or..not.active(ica,ihpvx,ijoin)
            end do        
          end do
        end do
        if (delete(idef)) ndel = ndel+1
      end do

      ! remove occ.s if needed:
      if (ndel.gt.0) then
        jdef = 0
        do idef = 1, ndef
          if (delete(idef)) cycle
 
          jdef = jdef+1
          if (jdef.lt.idef) then
            occ(1:ngastp,1:2,1:njoined,jdef) =
     &          occ(1:ngastp,1:2,1:njoined,idef)
          end if

        end do
      end if

      ndef = ndef-ndel

      deallocate(delete)
      return
      end

