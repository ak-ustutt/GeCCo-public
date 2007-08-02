*----------------------------------------------------------------------*
      logical recursive function next_fact(iperm,imult,nel,nuniq,
     &     iconn)
*----------------------------------------------------------------------*
*     generate next contraction sequence (factorization)
*     comment on algorithm: "nicht schoen aber selten ...."
*     to get the basic idea, cf. recursive algorithm in permute.f
*     here, the possible permutations are cut by considering the
*     number of arcs eliminated from a graph when one contraction
*     of vertices is done (i.e. fusion of contraction indices)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel, nuniq,
     &     iconn(3,nuniq), imult(nel)
      integer, intent(inout) ::
     &     iperm(nel)

      logical ::
     &     succ
      integer ::
     &     inum1, inum2, imlt1, imlt2, ivtxnew

      integer ::
     &     icoscr(3,nuniq), iscr(nel), iscr2(nel),
     &     iperm_next(nel), imult_next(nel), iarceq(2,nuniq),
     &     nfused,  idx, jdx, kdx, idxn,
     &     ipos, iposn, ifuse, imin, idxmin, ibase,
     &     nuniq_next, ioff, nel_next, idxls
      

      integer, external ::
     &     maxvtx, idxlist, imltlist

      if (ntest.ge.100) then
        write(luout,*) '------------------------------'
        write(luout,*) ' next_fact speaks recursively'
        write(luout,*) '------------------------------'
        write(luout,*) ' current nel = ',nel
        write(luout,*) ' current nuniq = ',nuniq
        write(luout,'(x,a,20i3)') '  iperm = ',iperm(1:nel)
        write(luout,'(x,a,20i3)') '  imult = ',imult(1:nel)
      end if

      ! two unique elements is always easy
      if (nuniq.eq.2) then
        if (ntest.ge.100)
     &       write(luout,*) ' processing n==2 case ...'
        ! read elements and their multiplicity
        inum1 = iperm(1)
        imlt1 = imult(1)
        inum2 = iperm(imlt1+1)
        imlt2 = imult(imlt1+1)
        ! if still in ordered sequence ...
        if (inum1.lt.inum2) then
          ! ... obtain next permutation
          iscr(1:imlt2) = iperm(imlt1+1:imlt1+imlt2)
          iscr(imlt2+1:imlt1+imlt2) = iperm(1:imlt1)
          iperm(1:nel) = iscr(1:nel)
          succ = .true.
        else
          succ = .false.
        end if
      ! more than two unique elements
      else if (nuniq.gt.2) then
        if (ntest.ge.100)
     &       write(luout,*) ' processing n>2 case ...'

        ! read first element and multiplicity
        inum1 = iperm(1)
        imlt1 = imult(1)

        if (ntest.ge.100)
     &       write(luout,*) 'getting reduced graph for ',inum1
        ! obtain graph reduced by the present element
        ! we need a new vertex number for the result
        ivtxnew = maxvtx(iconn,nuniq)+1
c        nfused(1:nuniq) = 0
        
        do idx = 1, nuniq
          iarceq(1,idx) = iconn(1,idx)
          iarceq(2,idx) = iconn(1,idx)
        end do
        call reduce_graph(icoscr,nuniq_next,iarceq,
     &       inum1,ivtxnew,iconn,nuniq,nuniq)
        
        if (ntest.ge.100)
     &       write(luout,*) 'nuniq_next = ', nuniq_next

        if (nuniq_next.ge.2) then
          ! set up next iperm imult
          if (ntest.ge.100) write(luout,*)
     &         'setting up iperm and imult for next level'

          ioff = imlt1
          nel_next = nel-imlt1
          ! on iscr, we set up some information
          ! on equivalent elements using iarceq
          do idxn = 1, nel_next
            idxls = idxlist(iperm(ioff+idxn),iarceq,nuniq,2)
            iscr(idxn) = iarceq(2,idxls)
          end do

c          idx = imlt1+1
          idx = 1 ! index in iscr/iperm 
          kdx = 1 ! index in iperm_next
          do while(kdx.le.nel_next)
            inum1 = iscr(idx)
            ! make sure that element was not already processed
            succ = .true.
            do jdx = 1, idx-1
              if (iscr(jdx).eq.inum1) succ = .false.
            end do
            if (.not.succ) then
              idx = idx+1
              cycle
            end if

            iperm_next(kdx) = inum1
            ! look after equivalent elements
            imlt1 = 1
            do jdx = idx, nel_next
              inum2 = iscr(jdx)
              if (inum1.eq.inum2.and.iperm(ioff+jdx).ne.inum1) then
                iperm_next(kdx+imlt1) = iperm(ioff+jdx) 
                imlt1 = imlt1+1 
              end if
            end do
            imult_next(kdx) = imlt1
            if (imlt1.gt.1) imult_next(kdx+1:kdx+imlt1-1) = 0
            kdx = kdx+imlt1
            idx = idx+1
          end do

          if (ntest.ge.100) then
            write(luout,*) 'call to next_fact with:'
            write(luout,*) 'iperm_next = ',iperm_next(1:nel_next)
            write(luout,*) 'imult_next = ',imult_next(1:nel_next)
          end if

          succ = next_fact(iperm_next,imult_next,nel-ioff,
     &         nuniq_next,icoscr)
        else
          succ = .false.
        end if

        if (ntest.ge.100)
     &       write(luout,*) 'succes status from rightmost permut.: ',
     &       succ

        if (succ) then
          iperm(ioff+1:nel) = iperm_next(1:nel-ioff)
        else
          ! find next free arc number on iperm
          if (ntest.ge.100)
     &         write(luout,*) 'trying to increment current 1st position'
          imin = huge(imin)
          idxmin = -1
          do idx = 2, nel
            if (iperm(idx).gt.iperm(1).and.iperm(idx).lt.imin
     &           .and.imult(idx).gt.0) then
              imin = iperm(idx)
              idxmin = idx
            end if
          end do
          succ = idxmin.gt.0
          if (succ) then
            if (ntest.ge.100)
     &           write(luout,*) 'imin = ',imin, ' ... now rearranging'
            if (ntest.ge.100) then
              write(luout,*) 'original iperm, imult:'
              write(luout,*) iperm(1:nel)
              write(luout,*) imult(1:nel)
            end if
            inum1 = iperm(1)
            imlt1 = imult(1)
            inum2 = iperm(idxmin)
            imlt2 = imult(idxmin)
            iperm_next(1:imlt2)=iperm(idxmin:idxmin+imlt2-1)
            imult_next(1:imlt2)=imult(idxmin:idxmin+imlt2-1)
            iperm_next(imlt2+1:imlt2+idxmin-1)=iperm(1:idxmin-1)
            imult_next(imlt2+1:imlt2+idxmin-1)=imult(1:idxmin-1)
c            iperm_next(idxmin+imlt2:nel)=iperm(idxmin+imlt2:nel)
            imult_next(idxmin+imlt2:nel)=imult(idxmin+imlt2:nel)
            iperm(1:idxmin+imlt2-1) = iperm_next(1:idxmin+imlt2-1)
c            imult(1:idxmin+imlt2-1) = imult_next(1:idxmin+imlt2-1)
            if (ntest.ge.100) then
              write(luout,*) 'iperm after swap'
              write(luout,*) iperm(1:nel)
              write(luout,*) 'imult_next:'
              write(luout,*) imult_next(1:nel)
            end if
            ! sort remaining list
            idx = imlt2+2
            ioff = imlt2
            do while(idx.le.nel)
              inum1 = iperm(idx)
              imlt1 = imult_next(idx)
              if (imlt1.eq.0) then
                idx = idx+1
                cycle
              end if
              iscr(1:imlt1) = iperm(idx:idx+imlt1-1)
              iscr2(1:imlt1) = imult_next(idx:idx+imlt1-1)
              jdx = idx-1
              do while (jdx.gt.ioff.and.
     &             (imult_next(jdx).eq.0.or.iperm(jdx).gt.inum1))
                imlt2 = imult_next(jdx)
                if (imlt2.ne.0.and.iperm(jdx).gt.inum1) then
                  ibase=jdx+imlt1-imlt2+1
                  call icp_save(iperm(jdx),iperm(ibase),imlt2)
                  call icp_save(imult_next(jdx),imult_next(ibase),imlt2)
                end if
                jdx = jdx-1
              end do
              imlt2 = imult_next(jdx)
              iperm(jdx+imlt2:jdx+imlt2-1+imlt1) = iscr(1:imlt1)
              imult_next(jdx+imlt2:jdx+imlt2-1+imlt1) = iscr2(1:imlt1)
              idx = idx+imlt1
            end do
            if (ntest.ge.100) then
              write(luout,*) 'iperm after sort'
              write(luout,*) iperm(1:nel)
              write(luout,*) 'imult_next (for checking): '
              write(luout,*) imult_next(1:nel)              
            end if
          end if

        end if

      else
        write(luout,*) 'programming error (obviously)...'
        stop 'programming error'
      end if
      
      next_fact = succ

      if (ntest.ge.100) then
        write(luout,*) 'result: '
        if (succ) then
          write(luout,*) 'Generated the following: nel = ',nel
          write(luout,'(x,a,20i3)') '  iperm = ',iperm(1:nel)
c          write(luout,'(x,a,20i3)') '  imult = ',imult(1:nel)
        else
          write(luout,*) 'No success! returning to higher level'
        end if
      end if

      return
      end
