*----------------------------------------------------------------------*
      logical function next_part_pair(init,restrict,ipart,
     &     inum,nsum,inummin,inummax)
*----------------------------------------------------------------------*
*     create next partitioning of a pair of integer numbers inum(1:2) 
*     into nsum summands (inummin.le.value.le.inummax for each element) 
*      -- returned on ipart(1:2,nsum)
*     the numbers are ordered in increasing sequence, with the convention
*     inum=(/j,i+1/) is for all i,j,k greater than inum=(/k,i/)
*     if init==.true. the first possible partitioning is created
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     init, restrict
      integer, intent(in) ::
     &     inum(2), nsum, inummin, inummax
      integer, intent(inout) ::
     &     ipart(2,nsum)

      logical ::
     &     succ, init_loc
      integer ::
     &     ipscr(nsum), ipscr2(2,nsum), idx, isum, ibase,
     &     inumidx, inummaxidx

      logical, external ::
     &     next_part_number

      if (ntest.ge.100) then
        write(lulog,*) '=========================='
        write(lulog,*) ' this is next_part_pair'
        write(lulog,*) '=========================='
        write(lulog,*) '  inum = /',inum(1),char(92)
        write(lulog,*) '        ',char(92),inum(2),'/'
        write(lulog,*) '  nsum,inummin,inummax: ',
     &                nsum,inummin,inummax
        if (init) then
          write(lulog,*) 'initiating ...'
        else
          write(lulog,*) 'partitioning on entry: '
          write(lulog,'(2x,">",10i4)') ipart(1,1:nsum)
          write(lulog,'(2x,">",10i4)') ipart(2,1:nsum)
        end if
      end if

      ! base for element ordering function:
      ibase = inummax*nsum+1

      ! create auxiliary array
      if (.not.init) then
        do idx = 1, nsum
          ipscr(idx) = ipart(1,idx)+ipart(2,idx)*ibase
        end do
      end if

      inumidx = inum(1)+inum(2)*ibase
      inummaxidx = inummax*(ibase+1)

      init_loc = init
      succ = .false.
      do while(next_part_number(init_loc,restrict,ipscr,
     &     inumidx,nsum,inummin,inummaxidx))

        init_loc = .false.
        succ = .true.
        do idx = 1, nsum
          ipscr2(1,idx) = mod(ipscr(idx),ibase)
          ipscr2(2,idx) = ipscr(idx)/ibase
          succ = succ.and.
     &         ipscr2(1,idx).le.inummax.and.ipscr2(2,idx).le.inummax
        end do
        if (succ) then
          ipart(1:2,1:nsum) = ipscr2(1:2,1:nsum)
          exit
        end if

      end do

      if (ntest.ge.100) then
        if (succ) then
          write(lulog,*) '  created partitioning: '
          write(lulog,'(2x,">",10i4)') ipart(1,1:nsum)
          write(lulog,'(2x,">",10i4)') ipart(2,1:nsum)
        else
          write(lulog,*) ' no (further) partitioning possible'
        end if
      end if

      next_part_pair = succ

      return
      end
