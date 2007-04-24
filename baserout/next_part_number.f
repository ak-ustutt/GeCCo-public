
*----------------------------------------------------------------------*
      logical function next_part_number(init,restrict,ipart,
     &     inum,nsum,inummin,inummax)
*----------------------------------------------------------------------*
*     create next partitioning of an integer number inum into nsum
*     summands (inummin.le.value.le.inummax) -- returned on ipart(nsum)
*     the numbers are ordered in increasing sequence
*     if init==.true. the first possible partitioning is created
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     init, restrict
      integer, intent(in) ::
     &     inum, nsum, inummin, inummax
      integer, intent(inout) ::
     &     ipart(nsum)

      logical ::
     &     succ
      integer ::
     &     ipscr(nsum), idx, isum 
      
      logical, external ::
     &     nxtptn_rc
      integer, external ::
     &     ifndmax

      if (ntest.ge.100) then
        write(6,*) '=========================='
        write(6,*) ' this is next_part_number'
        write(6,*) '=========================='
        write(6,*) '  inum,nsum,inummin,inummax: ',
     &                inum,nsum,inummin,inummax
        if (init) then
          write(6,*) 'initiating ...'
        else
          write(6,*) 'partitioning on entry: '
          write(6,'(2x,">",10i4)') ipart(1:nsum)
        end if
      end if

      ! start out: nsum-1 times the minimum, rest in last entry
      if (init) then
        ipscr(1:nsum-1)=inummin
        ipscr(nsum)=inum-(nsum-1)*inummin
        succ = .true.
      else
        ipscr(1:nsum)=ipart(1:nsum)
        succ = .false.
      end if

      if (ipscr(nsum).ge.inummin) then
        ip_loop: do
          ! check upper bounds:
          if (succ.and.ifndmax(ipscr,1,nsum,1).le.inummax) then
            ! return current partitioning ...
            ipart(1:nsum) = ipscr(1:nsum)
            exit ip_loop
          end if

          succ = nxtptn_rc(restrict,ipscr,inum,nsum,inummin)
          if (.not.succ) exit ip_loop

        end do ip_loop
      else
        ! no partitionings possible
        succ = .false.
      end if

      if (ntest.ge.100) then
        if (succ) then
          write(6,*) '  created partitioning: '
          write(6,'(2x,">",10i4)') ipart(1:nsum)
        else
          write(6,*) ' no (further) partitioning possible'
        end if
      end if

      next_part_number = succ

      return
      end

*----------------------------------------------------------------------*
      logical recursive function nxtptn_rc(restrict,
     &     ipart,inum,nsum,inummin)
*----------------------------------------------------------------------*
*     recursive slave function for generating next partitioning of inum
*----------------------------------------------------------------------*

      implicit none

      logical, intent(in) ::
     &     restrict
      integer, intent(in) ::
     &     inum, nsum, inummin
      integer, intent(inout) ::
     &     ipart(nsum)

      logical ::
     &     succ
      integer ::
     &     isumm1

      if (nsum.eq.2) then

        ! n.eq.2 is easy:
        succ = ((.not.restrict.and.ipart(2)-1.ge.inummin).or.
     &       (ipart(1)+1).le.(ipart(2)-1))
        if (succ) then
          ipart(1) = (ipart(1)+1)
          ipart(2) = (ipart(2)-1)
        end if

      else if (nsum.gt.2) then

        ! sum of first n-1 elements:
        isumm1 = sum(ipart(1:nsum-1))
        ! get next partitioning for these elements
        succ = nxtptn_rc(restrict,ipart,isumm1,nsum-1,inummin)
        ! if unsuccessful ...
        if (.not.succ) then
          ! ... increase sum of first n-1 elements: 
          ipart(1:nsum-2) = inummin
          ipart(nsum-1) = isumm1+1-(nsum-2)*inummin
          ! ... and decrease last element
          ipart(nsum) = ipart(nsum)-1
          succ = ipart(nsum).ge.inummin
          ! check, whether all ordering conditions can be satisfied:
          if (restrict) then
            do while (ipart(nsum-1).gt.ipart(nsum))
              ! another recursive call 
              succ = nxtptn_rc(restrict,ipart,isumm1+1,nsum-1,inummin)
              if (.not.succ) exit
            end do
          end if
        end if

      else
        stop 'illegal value of nsum'
      end if

      nxtptn_rc = succ

      return
      end
