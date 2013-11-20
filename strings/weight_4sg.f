*----------------------------------------------------------------------*
      subroutine weight_4sg(iy4sg,iw4sg,ientry,
     &     norb,igamorb,nel,maxms,ngam)
*----------------------------------------------------------------------*
*     set up arc-weights on iy4sg (and vertex-weights on iw4sg) for
*     a four-slope graph
*
*     iy4sg(1:3,n,ms,igam,k) is the weight of the arc for occupying
*                           orbital k with alpha(+1),beta(-1),2 el.s
*                           when up to now n orbitals were filled
*                           giving M_s=ms and an IRREP igam
*     iw4sg(n,ms,igam,k)    is the weight of the vertex reached at orb.
*                           level k after filling n orbitals 
*                           giving M_s=ms and an IRREP igam
*
*     ientry(ms,igam)       is the entry point; should be set to
*                           ientry(0,1) = 1, all other 0, if new
*                           graph starts. Else, the exit points of
*                           the preceeding graph are passed
*
*     andreas, june 2006
*
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel,                 ! number of electrons to be considered
     &     maxms,               ! ms levels to be considered
     &     ientry(-maxms:maxms,ngam),  ! weigths of entry point
     &     norb,    ! number of orbital levels
     &     ngam,    ! number of IRREPS (1, 2, 4, or 8)
     &     igamorb(norb)     ! IRREP of orbitals

      integer, intent(out) ::
     &     iy4sg(3,0:nel,-maxms:maxms,ngam,norb),
     &     iw4sg(0:nel,-maxms:maxms,ngam,norb)


*----------------------------------------------------------------------*

      integer ::
     &     nelmin, nelmax, msmin, msmax,
     &     igammin, igammax, igammin_, igammax_,
     &     iorb, igam, ims, iel, ispc, ia,
     &     ielnext, igamnext, imsnext, ianum,
     &     ielprev, igamprev, imsprev, 
     &     iw, iwcur, iw4sgprev

*----------------------------------------------------------------------*

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'weight_4sg')
        write(lulog,'(x,a,4i4)')
     &       'norb, nel, maxms, ngam: ',norb, nel, maxms, ngam
        write(lulog,*) 'igamorb:'
        write(lulog,'(x,2(x,5i5))') igamorb(1:norb)
      end if

      ! set everything to zero
      iy4sg(1:3,0:nel,-maxms:maxms,1:ngam,1:norb)=0
      iw4sg(0:nel,-maxms:maxms,1:ngam,1:norb)=0
      
      nelmin = 0
      nelmax = 0
      ispc = 1
      ! loop over orbital levels
      do iorb = 1, norb
        ! set up vertex weights using weights from previous level
        ! --> special for first level
        if (iorb.eq.1) then
c          iw4sg(0,0,1,1)=1
c          igam=igamorb(1)
c          iw4sg(1,+1,igam,1)=1
c          iw4sg(1,-1,igam,1)=1
c          if (nel.ge.2) iw4sg(2,0,1,1)=1
          nelmax = min(2,nel)
          ! to be improved:
          msmin = -maxms
          msmax = +maxms
          igammin=1
          igammax=ngam
        end if
c        else
          igammin_ = igammin
          igammax_ = igammax
          do igam = igammin_, igammax_
            do ims = msmin, msmax
              do iel = nelmin, nelmax
                if (iorb.gt.1) then
                  iw4sgprev=iw4sg(iel,ims,igam,iorb-1)
                else if (iel.eq.0) then
                  iw4sgprev=ientry(ims,igam)
                else
                  iw4sgprev=0
                end if
                ! add weights of previous level to vertices that
                ! can be reached by adding 0 to 2 electrons
                ia_loop: do ia = -1, 2
                  ielnext = iel+iabs(ia)
                  if (ielnext.gt.nel) cycle ia_loop
                  if (iabs(ia).eq.1) then
                    imsnext = ims+ia
                    igamnext = multd2h(igam,igamorb(iorb))
                  else
                    imsnext = ims
                    igamnext = igam
                  end if
                  if (iabs(imsnext).gt.maxms) then
                    cycle ia_loop
                  end if
                  igammin = min(igammin,igamnext)
                  igammax = max(igammax,igamnext)

                  iw4sg(ielnext,imsnext,igamnext,iorb) = 
     &                 iw4sg(ielnext,imsnext,igamnext,iorb) +
     &                 iw4sgprev
c     &                 iw4sg(iel,ims,igam,iorb-1)

                end do ia_loop

              end do  ! iel
            end do  ! ims
          end do  ! igam

          ! simple version
          msmin = max(-maxms,msmin-1)
          msmax = min(maxms,msmax+1)
          nelmax = min(nel,nelmax+2)

c        end if

*----------------------------------------------------------------------*

        ! set up arc weights using vertex weights from this 
        ! and previous level
c        if (iorb.eq.1) then
c
c          ! nothing to do as first set of weights is zero
c
c        else

          ! weight of 0 arc is 0
          ! weight of +1 arc is weight reached by arc 0
          ! weight of -1 arc is sum of weights reached by arc 0 and +1
          ! weight of +2 arc is sum of weights reached by arc 0, +1, and -1

          do igam = 1, ngam
            do ims = -maxms, maxms
              do iel = 0, nel
                iwcur = iw4sg(iel,ims,igam,iorb)
                if (iwcur.eq.0) cycle

                if (iorb.gt.1) then
                  iw = iw4sg(iel,ims,igam,iorb-1)
                else if (iel.eq.0) then
                  iw = ientry(ims,igam)
                else
                  iw = 0
                end if
                do ianum = 1, 3
                  if (ianum.eq.1) ia=+1
                  if (ianum.eq.2) ia=-1
                  if (ianum.eq.3) ia=+2
                  ielprev = iel-iabs(ia)
                  if (ielprev.lt.0) cycle
                  if (iabs(ia).eq.1) then
                    imsprev = ims-ia
                    igamprev = multd2h(igam,igamorb(iorb))
                  else
                    imsprev = ims
                    igamprev = igam
                  end if
                  if (iabs(imsprev).gt.maxms) cycle

                  iy4sg(ianum,ielprev,imsprev,igamprev,iorb) = iw 
                  
                  if (iorb.gt.1) then
                    iw = iw + iw4sg(ielprev,imsprev,igamprev,iorb-1)
                  else if (ielprev.eq.0) then
c                    if (ielprev.ne.0) stop 'error, ielprev.ne.0!'
                    iw = iw + ientry(imsprev,igamprev)
                  end if

                end do ! ianum

              end do ! iel
            end do ! ims
          end do ! igam
                
c        end if

      end do ! iorb

      if (ntest.ge.100) then
        write(lulog,*) 'Vertex weight array:'
        do iorb = 1, norb
          do ims = -maxms, maxms
            do igam = 1, ngam
              write(lulog,'(x,i2,"-",i1,",",i2,12i5,/,6x,12i5)') 
     &             iorb,igam,
     &             ims,iw4sg(0:nel,ims,igam,iorb)
            end do
          end do
          write(lulog,'(/)')
        end do

        write(lulog,*) 'Arc weight array:'
        do iorb = 1, norb
          do ims = -maxms, maxms
            do igam = 1, ngam
              write(lulog,
     &              '(x,i2,"-",i1,",",i2,4("("i4,",",i4,",",i4")"),'//
     &              '/,6x,4("("i4,",",i4,",",i4")"))')
     &             iorb,igam,
     &             ims,iy4sg(1:3,0:nel,ims,igam,iorb)
            end do
          end do
          write(lulog,'(/)')
        end do

      end if

      return
      
      end
