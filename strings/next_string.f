*----------------------------------------------------------------------*
      logical function next_string(idorb,idspn,idss,
     &     nidx,ms,igam,first,
     &     igas_restr,
     &     mostnd_cur,igamorb,
     &     nsym,ngas_cur)
*----------------------------------------------------------------------*
*     generate next string with for current Ms, IRREP
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nidx, ms, igam,
     &     igas_restr(2,ngas_cur,2),
     &     nsym,ngas_cur,
     &     mostnd_cur(2,nsym),igamorb(*)
      integer, intent(inout) ::
     &     idorb(nidx), idspn(nidx), idss(nidx)
      logical, intent(in) ::
     &     first

      logical ::
     &     succ
      integer ::
     &     idx, idxoff, idxnd, idxst
      integer ::
     &     ioss(ngas_cur)

      logical, external ::
     &     lexlstr, nextstr, next_ssd
      integer, external ::
     &     igamstr2, ielsum


      if (ntest.gt.0) then
        write(luout,*) '-----------------------'
        write(luout,*) ' info from next_string'
        write(luout,*) '-----------------------'
        write(luout,*) ' nidx = ',nidx
      end if

c dbg
c      print *,'idorb',idorb(1:nidx)
c dbg

      if (first) then

        next_string = .true.
        if (nidx.eq.0) return

        ! index range within tupel
        idxst = 1
        idxnd = nidx
        ! first subspace distribution
        idss(1:nidx) = 1
        dss_loop: do
          ! lexically lowest string
          ! reform distribution to occupation
          ioss(1:ngas_cur) = 0
          do idx = idxst, idxnd
            ioss(idss(idx)) = ioss(idss(idx))+1
          end do

          succ = lexlstr(nidx,ms,
     &           ioss,idorb,idspn,
     &           mostnd_cur,nsym,ngas_cur)

          if (succ) then
            str_loop: do
              ! check symmetry
c dbg
c            print *,'nidx,idorb: ',nidx,idorb(1:nidx)
c            print *,'igamorb',igamorb(1:8)
c dbg
              succ = igamstr2(nidx,idorb,igamorb).eq.igam
              ! exit if successful
              if (succ) exit dss_loop
              ! else: get next possible string
              succ = nextstr(nidx,ms,
     &               idss,idorb,idspn,
     &               mostnd_cur,nsym,ngas_cur) 

              if (.not.succ) exit str_loop
            end do str_loop
          end if
          ! if no appropriate string found: 
          !  get next subspace distr.
          succ = next_ssd(idss,nidx,nidx,
     &             ngas_cur,igas_restr)
          ! nothing: we have to give up
          if (.not.succ) exit dss_loop
        end do dss_loop

        next_string = succ

        if (ntest.ge.100) then
          if (succ) then
            write(luout,*) 'first string: ',idorb(1:nidx)
            write(luout,*) '              ',idspn(1:nidx)
            write(luout,*) '              ',idss(1:nidx)
          else
            write(luout,*) 'no string exists'
          end if
        end if

        return
      else

        next_string = .false.
        if (nidx.eq.0) return

        if (ntest.ge.100) then
          write(luout,*) ' input string:',idorb(1:nidx)
          write(luout,*) '              ',idspn(1:nidx)
          write(luout,*) '              ',idss(1:nidx)
        end if

        dss_loop2: do
          str_loop2: do
            ! try next string ...
            succ = nextstr(nidx,ms,
     &           idss,idorb,idspn,
     &           mostnd_cur,nsym,ngas_cur) 

            ! ... no further string ? ....
            if (.not.succ) exit str_loop2

            ! ... else check symmetry ...
            if (succ)
     &           succ = igamstr2(nidx,idorb,igamorb).eq.igam

            ! if symmetry was alright, we are done
            if (succ) exit dss_loop2
          end do str_loop2
                
          ! here we land, if we need to try the next
          ! subspace distribution
          succ = next_ssd(idss,nidx,nidx,
     &             ngas_cur,igas_restr)
          ! nothing: we have to give up
          if (.not.succ) exit dss_loop2

          ! else we have to get the lexically lowest string
          ! for the new subspace distribution

          ! reform distribution to occupation
          ioss(1:ngas_cur) = 0
          do idx = idxst, idxnd
            ioss(idss(idx)) = ioss(idss(idx))+1
          end do

          succ = lexlstr(nidx,ms,
     &             ioss,idorb,idspn,
     &             mostnd_cur,nsym,ngas_cur)

c dbg
c            print *,'nidx,idorb: ',nidx,idorb(1:nidx)
c            print *,'igamorb',igamorb(1:8)
c dbg

          ! ... check symmetry ...
          if (succ)
     &           succ = igamstr2(nidx,idorb,igamorb).eq.igam

          ! ... we are done?
          if (succ) exit dss_loop2

          ! ... else go to top of loop and increment string

        end do dss_loop2
        
        next_string = succ

        if (ntest.ge.100) then
          if (succ) then
            write(luout,*) ' next string: ',idorb(1:nidx)
            write(luout,*) '              ',idspn(1:nidx)
            write(luout,*) '              ',idss(1:nidx)
          else
            write(luout,*) ' no further string exists'
          end if
        end if

        return

      end if

      end
*----------------------------------------------------------------------*

