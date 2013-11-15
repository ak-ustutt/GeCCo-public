*----------------------------------------------------------------------*
      subroutine set_spprjmap_kernel(strmap,
     &     iocc,irestr,idxms,igam,
     &     g_y4sg,g_yinf,
     &     g_yssg,g_wssg,
     &     g_off_dgm,g_ndis,
     &     mostnd_cur,nsym,ngas_cur,igamorb)
*----------------------------------------------------------------------*
*     core routine for setting up the (str)->(str_other_spins) mapping
*
*     get for -+ -> ++ +- --
*     get for -++ -> +++ +-+ ++- --+ -+- +-- --- 
*
*     andreas, april 2013
*      
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000

      integer, parameter ::
     &     idspn_prm2(6) = (/+1,+1,
     &                       +1,-1,
     &                       -1,-1/)
      integer, parameter ::
     &     idspn_prm3(21) = (/+1,+1,+1,
     &                        +1,-1,+1,
     &                        +1,+1,-1,
     &                        -1,-1,+1,
     &                        -1,+1,-1,
     &                        +1,-1,-1,
     &                        -1,-1,-1/)

      integer, intent(out) ::
     &     strmap(*)
      integer, intent(in) ::
     &     iocc, irestr(2,ngas_cur,2), idxms, igam,
     &     nsym, ngas_cur,
     &     mostnd_cur(2,nsym),igamorb(*),
     &     g_y4sg(*),g_yinf(*),
     &     g_yssg(*),g_wssg(*),
     &     g_off_dgm(*),g_ndis

      logical ::
     &     first
      integer ::
     &     idorb(iocc),idspn(iocc),idgam(iocc),idss(iocc),
     &     idspn_flipped(iocc),
     &     istr, nstr, idx, ms, idxmap, nmaps, imaps

      integer, external ::
     &     idx4sg, std_spsign
      logical, external ::
     &     next_string

      if (iocc.eq.1) 
     &     call quit(1,'set_spprjmap_kernel','occ=1: use trivial map!')
      if (iocc.gt.3)
     &     call quit(1,'set_spprjmap_kernel','sorry, only occ=2 or 3')

      nmaps = 2**iocc-1
      
      ms = iocc-(idxms-1)*2

      if (iocc.eq.2.and.ms.ne.0)
     & call quit(1,'set_spprjmap_kernel','sorry, only ms=0 for occ=2')
      if (iocc.eq.3.and.ms.ne.1)
     & call quit(1,'set_spprjmap_kernel','sorry, only ms=1/2 for occ=3')

      idxmap = 0

      first = .true.
      ! loop over original strings
      do while(next_string(idorb,idspn,idss,
     &     iocc,ms,igam,first,
     &     irestr,
     &     mostnd_cur,igamorb,
     &     nsym,ngas_cur))

        first = .false.

        ! must be one of: -+, 22  or -++, 22+ 
        if (idspn(1).ne.-1.and.idspn(1).ne.2) then
          strmap(idxmap+1:idxmap+nmaps) = 0
          idxmap = idxmap+nmaps
          cycle
        end if
        
        do idx = 1, iocc
          idgam(idx) = igamorb(idorb(idx))
        end do

        do imaps = 1, nmaps
          ! flip spin indices
          if (iocc.eq.2) then
            idspn_flipped(1:2) = 
     &               idspn_prm2(((imaps-1)*2)+1:((imaps-1)*2)+2)
          else
            idspn_flipped(1:3) = 
     &               idspn_prm3(((imaps-1)*3)+1:((imaps-1)*3)+3)
          end if

          ! check for paired indices:
          ! they can never become a ++ or -- (Pauli says)
          if (idspn(1).eq.2) then
            if (idspn_flipped(1).eq.idspn_flipped(2)) then
              idxmap = idxmap+1
              strmap(idxmap) = 0
              cycle
            end if
          end if
          if (iocc.eq.3.and.idspn(2).eq.2) then
            if (idspn_flipped(2).eq.idspn_flipped(3)) then
              idxmap = idxmap+1
              strmap(idxmap) = 0
              cycle
            end if
          end if
          ! reset all paired indices to 22 notation
          do idx = 1, iocc
            if (idspn(idx).eq.2) idspn_flipped(idx)=2
          end do

          idxmap = idxmap + 1
          ! get sign change
          strmap(idxmap) = std_spsign(idspn        ,iocc)*
     &                   std_spsign(idspn_flipped,iocc)
c dbg
c        print *,'idorb:   ',idorb(1:iocc)
c        print *,'idspn:   ',idspn(1:iocc)
c        print *,'idss:    ',idss (1:iocc)
c        print *,'idspn_fl:',idspn_flipped(1:iocc)
c        print *,'sign: ',strmap(idxmap)
c dbg
          ! get new address
          strmap(idxmap) = strmap(idxmap)*
     &        (idx4sg(iocc,idss,idorb,idspn_flipped,idgam,
     &                g_y4sg,g_yinf,
     &                g_yssg,g_wssg,
     &                g_off_dgm,g_ndis,
     &                mostnd_cur,
     &                iocc,nsym,ngas_cur)+1)

        end do
c dbg
c        print *,'idx = ',abs(strmap(idxmap))
c dbg

      end do

      if (ntest.ge.1000) then
        write(luout,*) 'length = ',idxmap
        write(luout,*) 'set_spprjmap_kernel: the result of my work is'
        nstr = idxmap/nmaps
        do istr = 1, nstr
          idx = (istr-1)*nmaps+1
          write(luout,'(x,i5,"->",20i5)') istr,strmap(idx:idx-1+nmaps)
        end do
      end if

      return
      end
