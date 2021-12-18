*----------------------------------------------------------------------*
      subroutine set_spprjmap_kernel(strmap,
     &     iocc,idspn_prm,irestr,idxms,igam,
     &     g_y4sg,g_yinf,
     &     g_yssg,g_wssg,
     &     g_off_dgm,g_ndis,
     &     mostnd_cur,nsym,ngas_cur,igamorb)
*----------------------------------------------------------------------*
*     core routine for setting up the (str)->(str_other_spins) mapping
*
*     get for -+ -> ++ +- --
*     get for +-+ -> +++ ++- +-- -++ -+- --+ ---
*     (not for -++, because orbital pairing +22 would then be forbidden)
*
*     andreas, april 2013
*      
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000

      integer, intent(out) ::
     &     strmap(*)
      integer, intent(in) ::
     &     iocc, irestr(2,ngas_cur,2), idxms, igam,
     &     nsym, ngas_cur, idspn_prm(0:iocc,0:2**iocc-1),
     &     mostnd_cur(2,nsym),igamorb(*),
     &     g_y4sg(*),g_yinf(*),
     &     g_yssg(*),g_wssg(*),
     &     g_off_dgm(*),g_ndis

      logical ::
     &     first
      integer ::
     &     idorb(iocc),idspn(iocc),idgam(iocc),idss(iocc),
     &     idspn_flipped(iocc),
     &     istr, nstr, idx, ms, idxmap, nmaps, imaps, isgn

      integer, external ::
     &     idx4sg, std_spsign
      logical, external ::
     &     next_string

      if (iocc.eq.1) 
     &     call quit(1,'set_spprjmap_kernel','occ=1: use trivial map!')

      nmaps = 2**iocc-1
      
      ms = iocc-(idxms-1)*2

      if (mod(iocc,2).eq.0.and.ms.ne.0)
     & call quit(1,'set_spprjmap_kernel',
     &           'sorry, only ms=0 for even occ')
      if (mod(iocc,2).eq.1.and.ms.ne.1)
     & call quit(1,'set_spprjmap_kernel',
     &           'sorry, only ms=1/2 for odd occ')

      idxmap = 0

      first = .true.
      ! loop over original strings
      str_loop: do while(next_string(idorb,idspn,idss,
     &     iocc,ms,igam,first,
     &     irestr,
     &     mostnd_cur,igamorb,
     &     nsym,ngas_cur))

        first = .false.

        ! must be alternatingly - and + (always + at the right end)
        isgn = +1
        do idx = iocc, 1, -1
          if (idspn(idx).ne.isgn.and.idspn(idx).ne.2) then
            strmap(idxmap+1:idxmap+nmaps) = 0
            idxmap = idxmap+nmaps
            cycle str_loop
          end if
          isgn = -isgn
        end do
        
        do idx = 1, iocc
          idgam(idx) = igamorb(idorb(idx))
        end do

        maps_loop: do imaps = 1, nmaps
          ! flip spin indices
          idspn_flipped(1:iocc) = idspn_prm(1:iocc,imaps)

          ! reset all paired indices to 22 notation and
          ! check for paired indices:
          ! they can never become a ++ or -- (Pauli says)
          do idx = 1, iocc
            if (idspn(idx).eq.2.and.idspn_flipped(idx).ne.2) then
              if (idspn(idx+1).ne.2) call quit(1,'set_spprjmap_kernel',
     &                                         'broken spin pair?')
              if (idspn_flipped(idx)*idspn_flipped(idx+1).ne.-1) then
                idxmap = idxmap+1
                strmap(idxmap) = 0
                cycle maps_loop
              end if
              idspn_flipped(idx:idx+1)=2
            end if
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

        end do maps_loop
c dbg
c        print *,'idx = ',abs(strmap(idxmap))
c dbg

      end do str_loop

      if (ntest.ge.1000) then
        write(lulog,*) 'length = ',idxmap
        write(lulog,*) 'set_spprjmap_kernel: the result of my work is'
        nstr = idxmap/nmaps
        do istr = 1, nstr
          idx = (istr-1)*nmaps+1
          write(lulog,'(x,i5,"->",20i5)') istr,strmap(idx:idx-1+nmaps)
        end do
      end if

      return
      end
