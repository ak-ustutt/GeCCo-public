*----------------------------------------------------------------------*
      logical function nextstr(nel,itms,idss,idorb,idspn,
     &     mostnd,ngam,nspc)
*----------------------------------------------------------------------*
*     generate lexically next string for given number of elements
*     (nel), total Ms, and subspace distribution idss(1:nel)
*     the total IRREP is not considered, call function igamstr to
*     set the IRREP distribution and get the total IRREP, and skip
*     if the total irrep is not correct
*----------------------------------------------------------------------*

      implicit none
       
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel, itms, ngam, nspc,
     &     idss(nel), mostnd(2,ngam,nspc)
      integer, intent(inout) ::
     &     idorb(nel), idspn(nel)

      logical ::
     &     lnext, lsame
      integer ::
     &     nalph, nbeta, ipos, ipairc, nelsb, ms_pos, ms_sub,
     &     idorbc, idspnc, idorbn, idspnn, iorbmax,
     &     iosssb(nspc)

      logical, external ::
     &     lexlstr

      if (ntest.gt.0) then
        write(lulog,*) '---------------'
        write(lulog,*) 'This is nextstr'
        write(lulog,*) '---------------'
      end if
      if (ntest.ge.50) then
        write(lulog,*) 'ON ENTRY:'
        write(lulog,*) ' nel = ',nel
        write(lulog,*) ' itms = ',itms
        write(lulog,*) ' idss = ',idss(1:nel)
        write(lulog,*) ' idorb = ',idorb(1:nel)
        write(lulog,*) ' idspn = ',idspn(1:nel)
      end if

      ipos = 1
      ms_pos = 0
      ipairc = 0 ! counter for pairs
      lnext = .true.
      lsame = .false. ! flag whether same ipos is visited again
      string_loop: do

        idorbc = idorb(ipos)
        idspnc = idspn(ipos)

        if (.not.lsame) then      
          ! get Ms of components 1 through ipos:
          if (idspnc.eq.2.and.ipairc.eq.0) then
            ! no change of Ms for doubly occupied orbital
            ! we can skip to next index directly
            ipos = ipos+1
            ipairc=1
            cycle string_loop
          else if (idspnc.eq.2.and.ipairc.eq.1) then
            ipairc=0
            ! we treat this orbital as beta case
            idspnc=-1
          else
            ms_pos = ms_pos + idspnc
          end if
        else
          lsame = .false.
        end if

        !  -->  number of allowed alphas and allowed betas
        nalph = (ipos+ms_pos)/2
        nbeta = (ipos-ms_pos)/2

        iorbmax=mostnd(2,ngam,idss(ipos))
        if (ipos.lt.nel) idorbn = idorb(ipos+1)
c dbg fix by mh
        if (ipos.ge.nel) idorbn = 100
c dbg end fix
        if (ipos.lt.nel) idspnn = idspn(ipos+1)
        if (ntest.ge.100) then
          write(lulog,*) 'ipos,nalph,nbeta: ',ipos,nalph,nbeta
          write(lulog,*) 'iorbmax: ',iorbmax
        end if

        ! possibility a) alpha->beta
        if (idspnc.eq.+1.and.nbeta.gt.0) then
          if (ntest.ge.100) then
            write(lulog,*) 'case 1'
          end if
          idspn(ipos)=-1
          ! Ms of substring
          ms_sub = ms_pos+1
c          exit string_loop
        ! else we try to increase orbital counter
        ! difference to next orbital > 1 or last ipos
        else if ( (ipos.eq.nel.or.idorbc.lt.idorbn-1)
     &            .and.idorbc.lt.iorbmax ) then
          ! increase orbital component ...
          idorb(ipos) = idorb(ipos)+1
          if (idspnc.eq.-1.and.nalph.gt.0) then
            ! ... and go to alpha
            if (ntest.ge.100) then
              write(lulog,*) 'case 2a'
            end if
            idspn(ipos)=+1
            ms_sub = ms_pos-1
          else
            if (ntest.ge.100) then
              write(lulog,*) 'case 2b'
            end if
            ms_sub = ms_pos-idspnc
          end if
c          exit string_loop
          ! difference is only one and next position has
          ! beta spin (IMPORTANT) --> make a pair
        else if (idorbc.lt.idorbn.and.idorbc.lt.iorbmax.and.
     &         idspnn.eq.-1.and.nalph.gt.0) then
          if (ntest.ge.100) then
            write(lulog,*) 'case 3'
          end if
          idorb(ipos) = idorb(ipos)+1
          idspn(ipos)=2
          idspn(ipos+1)=2
          ms_sub = ms_pos-1
c          exit string_loop
        else if (ipos.eq.nel) then
          ! no further increment possible
          if (ntest.ge.100) then
            write(lulog,*) 'case x'
          end if
          lnext = .false.
          exit string_loop
        else
          ! go to next ipos
          if (ntest.ge.100) then
            write(lulog,*) 'case n'
          end if
          ipos = ipos+1
          cycle string_loop
        end if

        ! this part of loop is only entered, if a index
        ! at ipos was changed ...
        if (ipos.eq.1) exit string_loop
        ! ... call lexlstr for leftmost remainder
        nelsb = ipos-1
        iosssb(1:nspc) = 0
        do ipos = 1, nelsb
          iosssb(idss(ipos)) = iosssb(idss(ipos))+1
        end do
        lnext = lexlstr(nelsb,ms_sub,iosssb,idorb,idspn,
     &       mostnd,ngam,nspc)
        ! test whether nelsb and nelsb+1 are paired now
        if (lnext.and.idspn(nelsb+1).eq.-1.and.
     &       (idorb(nelsb).eq.idorb(nelsb+1))) then
          idspn(nelsb:nelsb+1) = 2
        end if

        ! if everything OK we can go home ...
        if (lnext) exit string_loop
        ! else we have to increase counter starting 
        ! at THE SAME ipos as before
        lsame = .true.
                
      end do string_loop
      
      nextstr = lnext

      if (ntest.ge.50) then
        write(lulog,*) 'ON EXIT from nextstr:'
        if (lnext) then
          write(lulog,*) ' idorb = ',idorb(1:nel)
          write(lulog,*) ' idspn = ',idspn(1:nel)
        else
          write(lulog,*) ' no next string !'
        end if
      end if

      return
      end 
