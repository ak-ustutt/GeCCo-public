      logical function next_part_triple(init,restrict,ipart,
     &     inum,nsum,inummin,inummax)
*-----------------------------------------------------------------------*
*     Create next partitioning of a triplet of integer numbers, inum(1:3)
*     into nsum summands (inummin.le.value.le.inummax for each element)
*     returned in ipart(1:3,nsum).
*     
*     If init=.true. then the first possible partitioning is created.
*-----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest=00
      logical, intent(in) ::
     &     init,restrict
      integer, intent(in) ::
     &     inum(3),nsum,inummin,inummax
      integer, intent(inout) ::
     &     ipart(3,nsum)

      logical ::
     &     succ,init_loc
      integer ::
     &     idx,inumidx,inummaxidx,imod1,imod2,ip,
     &     ipscr(nsum),ipscr2(3,nsum),isum1,isum2,isum3

      logical, external ::
     &     next_part_number

      ! Define the moduli used to encode the members of inum into a 
      ! single number.
      imod1=inum(1)+1
      inumidx=inum(1)+imod1*inum(2)
      imod2=max(inumidx+1,imod1+1)

      inumidx=inumidx+inum(3)*imod2
      inummaxidx=inummax*(imod1+imod2+1)

      if(.not.init)then
        do idx=1,nsum
          ipscr(idx)=ipart(1,idx)+ipart(2,idx)*imod1+ipart(3,idx)*imod2
        enddo
      endif  

      init_loc=init
      succ=.false.
      pn_loop : do while(next_part_number(init_loc,restrict,ipscr,
     &     inumidx,nsum,inummin,inummaxidx))
        
        init_loc=.false.
        succ=.true.

        isum1=0
        isum2=0
        isum3=0
        do idx=1,nsum
          ipscr2(3,idx)=ipscr(idx)/imod2
          ip=ipscr(idx)-ipscr2(3,idx)*imod2
          ipscr2(2,idx)=ip/imod1
          ipscr2(1,idx)=mod(ip,imod1)
          isum1=isum1+ipscr2(1,idx)
          isum2=isum2+ipscr2(2,idx)
          isum3=isum3+ipscr2(3,idx)

          succ=succ.and.ipscr2(1,idx).le.inummax.and.
     &         ipscr2(2,idx).le.inummax.and.ipscr2(3,idx).le.inummax
        enddo  
        if(isum1.ne.inum(1).or.isum2.ne.inum(2).or.isum3.ne.inum(3))then
          if(ipscr(1).eq.(inumidx-nsum+1))then
            succ=.false.
            exit pn_loop
          else  
            cycle pn_loop
          endif
        endif
        if(succ)then
          ipart(1:3,1:nsum)=ipscr2(1:3,1:nsum)
          exit pn_loop
        endif  
      enddo  pn_loop

      if (ntest.ge.100) then
        if (succ) then
          write(lulog,*) '  created partitioning: '
          write(lulog,'(2x,">",10i4)') ipart(1,1:nsum)
          write(lulog,'(2x,">",10i4)') ipart(2,1:nsum)
          write(lulog,'(2x,">",10i4)') ipart(3,1:nsum)
        else
          write(lulog,*) ' no (further) partitioning possible'
        end if
      end if

c      succ=succ.and.(ipscr(1).ne.(inumidx-nsum+1))
      next_part_triple = succ

      return
      end
