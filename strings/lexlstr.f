*----------------------------------------------------------------------*
      logical function lexlstr(nel,itms,ioss,idorb,idspn,
     &     mostnd,ngam,nspc)
*----------------------------------------------------------------------*
*     generate lexically lowest string for given number of elements
*     (nel), total Ms, and subspace occupation ioss(1:nspc)
*     the total IRREP is not considered, call function igamstr to
*     set the IRREP distribution and get the total IRREP, see also
*     nextstr
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel, itms, ngam, nspc,
     &     ioss(nspc),mostnd(2,ngam,nspc)
      integer, intent(out) ::
     &     idorb(nel), idspn(nel)

      logical ::
     &     succ
      integer ::
     &     npair, nlone, nalpht, nbetat, ispc,
     &     ipos, ipos0, inum, ipair, ilone,
     &     maxnum, lenstr, npair4ss, nalph4ss, nbeta4ss,
     &     npair_res(nspc)


      if (ntest.gt.0) then
        write(luout,*) '---------------'
        write(luout,*) 'This is lexlstr'
        write(luout,*) '---------------'
      end if
      if (ntest.ge.50) then
        write(luout,*) 'ON ENTRY:'
        write(luout,*) ' nel = ',nel
        write(luout,*) ' itms = ',itms
        write(luout,*) ' ioss = ',ioss(1:nspc)
      end if

      succ = .true.

      ! get maximum number of paired electrons
      npair = (nel-abs(itms))/2
      ! rest is purely alpha or beta, currently
      nalpht = npair
      nbetat = nalpht
      if (itms.gt.0) nalpht = nalpht+itms
      if (itms.lt.0) nbetat = nbetat-itms

      ! scan size of subspaces --> maybe there is a min. number
      ! of pairs necessary to form an allowed string at all
      npair = 0
      do ispc = 1, nspc
        maxnum = mostnd(2,ngam,ispc)-mostnd(1,1,ispc)+1
        lenstr = ioss(ispc)
        ! store number of pairs that are reserved by lower spaces
        npair_res(ispc) = npair
        if (maxnum.lt.lenstr) then
          npair4ss = lenstr-maxnum
          if (npair4ss.gt.maxnum) then
            ! impossible to generate any string at all!
            succ=.false.
            exit
          end if
          npair = npair+npair4ss
        end if
      end do

      if (succ) then
        ! loop over subspaces
        ipos0 = nel+1
        ispc_loop: do ispc = nspc, 1, -1
          if (ioss(ispc).eq.0) cycle ispc_loop          
          inum = mostnd(1,1,ispc)
          ipos0 = ipos0-ioss(ispc)
          ipos = ipos0
          ! paired electrons ...
          npair = min(nalpht,nbetat)-npair_res(ispc)  ! possible pairs
          if (npair.lt.0) then
            succ = .false.
            exit ispc_loop
          end if
          npair4ss = min(npair,ioss(ispc)/2) ! number of pairs in subspace
c        npair = npair-npair4ss
          nalpht = nalpht-npair4ss ! decrease counters
          nbetat = nbetat-npair4ss
          do ipair = 1, npair4ss
            idorb(ipos:ipos+1) = inum
            idspn(ipos:ipos+1) = 2
            ipos = ipos + 2
            inum = inum + 1
          end do
          ! ... and their unpaired cousins:
          nlone = ioss(ispc)-2*npair4ss ! number of lone electrons
          nalph4ss = min(nalpht,nlone) ! as many alphas as possible
          nbeta4ss = nlone-nalph4ss ! rest as betas
          if (min(nalpht-nalph4ss,nbetat-nbeta4ss).lt.npair_res(ispc))
     &         then
            npair = npair_res(ispc)-min(nalpht-nalph4ss,nbetat-nbeta4ss)
            nalph4ss = nalph4ss-npair
            nbeta4ss = nalph4ss+npair
          end if
          nalpht = nalpht-nalph4ss ! decrease counters
          nbetat = nbetat-nbeta4ss
          do ilone = 1, nbeta4ss
            idorb(ipos) = inum
            idspn(ipos) = -1
            inum = inum+1
            ipos = ipos+1
          end do
          do ilone = 1, nalph4ss
            idorb(ipos) = inum
            idspn(ipos) = +1
            inum = inum+1
            ipos = ipos+1
          end do
          ! was last orbital number larger than allowed maximum?
          if (inum-1.gt.mostnd(2,ngam,ispc)) then
            succ = .false.
            exit ispc_loop
          end if
        end do ispc_loop

      end if

      lexlstr = succ

      if (ntest.ge.50) then
        write(luout,*) 'ON EXIT from lexlstr:'
        if (succ) then
          write(luout,*) ' idorb = ',idorb(1:nel)
          write(luout,*) ' idspn = ',idspn(1:nel)
        else
          write(luout,*) ' no lexically lowest string possible!'
        end if
      end if
      
      return
      end
