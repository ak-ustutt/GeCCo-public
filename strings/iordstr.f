*----------------------------------------------------------------------*
      integer function iordstr(idorb12,idspn12,
     &     idorb1,idspn1,n1,
     &     idorb2,idspn2,n2)
*----------------------------------------------------------------------*
*     given are 2 strings -> reorder to rev.lex.order and give
*     sign factor (-1)^<number of permutations>
*     input strings must be ordered !
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     n1, n2,
     &     idorb1(n1), idspn1(n1), idorb2(n2), idspn2(n2)
      integer, intent(inout) ::
     &     idorb12(n1+n2), idspn12(n1+n2)
      
      integer ::
     &     iperm, ipos1, ipos2, ipos12, isp1sp2

      if (ntest.gt.0) then
        write(lulog,*) '--------------'
        write(lulog,*) 'this is ordstr'
        write(lulog,*) '--------------'
      end if
      if (ntest.ge.100) then
        write(lulog,*) 'input strings: '
        write(lulog,'(x,a,15i4)') ' 1 o: ',idorb1(1:n1)
        write(lulog,'(x,a,15i4)') '   s: ',idspn1(1:n1)
        write(lulog,'(x,a,15i4)') ' 2 o: ',idorb2(1:n2)
        write(lulog,'(x,a,15i4)') '   s: ',idspn2(1:n2)
      end if

      ! handle some special cases:
      iordstr = 1
      if (n1+n2.eq.0) return
      if (n2.eq.0) then
        idorb12(1:n1) = idorb1(1:n1)
        idspn12(1:n1) = idspn1(1:n1)
        return
      end if
      if (n1.eq.0) then
        idorb12(1:n2) = idorb2(1:n2)
        idspn12(1:n2) = idspn2(1:n2)
        return
      end if

      iordstr = 0

      ipos1 = 1
      ipos2 = 1
      iperm = 0
      ! better: this loop only until ipos1 or 2 is exhausted,
      ! the rest in post-loop
      do ipos12 = 1, n1+n2
        ! pairing requires this to be == -1:
        if (ipos1.le.n1.and.ipos2.le.n2) then
          isp1sp2 = idspn1(ipos1)*idspn2(ipos2)
          if (idorb1(ipos1).eq.idorb2(ipos2).and.isp1sp2.ne.-1) return
        end if

        if ( ipos2.gt.n2 .or.
     &      (ipos1.le.n1 .and.
     &       ( ( idorb1(ipos1).lt.idorb2(ipos2) .or.
     &          (idorb1(ipos1).eq.idorb2(ipos2) .and.
     &           idspn1(ipos1).eq.+1)          ) ) ) ) then
          idorb12(ipos12) = idorb1(ipos1)
          idspn12(ipos12) = idspn1(ipos1)
          ipos1 = ipos1+1
        else
          idorb12(ipos12) = idorb2(ipos2)
          idspn12(ipos12) = idspn2(ipos2)
          ipos2 = ipos2+1
          iperm = iperm+n1-ipos1+1
        end if

      end do
        
      ! post-processing for spin-coupling:
      do ipos1 = 1, n1+n2-1
        if (idorb12(ipos1).eq.idorb12(ipos1+1)) then
          idspn12(ipos1:ipos1+1) = 2
        end if
      end do

      iordstr = (-1)**iperm

      if (ntest.ge.100) then
        write(lulog,*) 'output value: ',iordstr
        if (iordstr.ne.0) then
          write(lulog,*) 'output string: '
          write(lulog,'(x,a,15i4)') ' > o: ',idorb12(1:n1+n2)
          write(lulog,'(x,a,15i4)') '   s: ',idspn12(1:n1+n2)
        end if
      end if

      return
      end
