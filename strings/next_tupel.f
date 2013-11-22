*----------------------------------------------------------------------*
      logical function next_tupel(idorb,idspn,idss,
     &     nidx,iocc,igrph,msdst,igamdst,first,
     &     igas_restr,
     &     mostnd,igamorb,
     &     nsym,ngas,ngas_hpv,idx_gas,
     &     ihpvseq,lexlscr)
*----------------------------------------------------------------------*
*     as next_tupel_ca, but only for either c or a alone
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nidx,
     &     iocc(ngastp), igrph(ngastp),
     &     msdst(ngastp), igamdst(ngastp),
     &     igas_restr(2,ngas,2,*),
     &     nsym,ngas,ngas_hpv(*),idx_gas(*),
     &     mostnd(2,nsym,ngas),igamorb(*),
     &     ihpvseq(*)
      integer, intent(inout) ::
     &     idorb(nidx), idspn(nidx), idss(nidx),
     &     lexlscr(nidx,3)
      logical, intent(in) ::
     &     first

      logical ::
     &     succ
      integer ::
     &     ihpvdx, ihpv, idx, idxoff(ngastp), idxnd, idxst,
     &     nel, ms, igam, igr
      integer ::
     &     ioss(ngas)

      logical, external ::
     &     lexlstr, nextstr, next_ssd
      integer, external ::
     &     igamstr2, ielsum


      if (ntest.gt.0) then
        write(lulog,*) '----------------------'
        write(lulog,*) ' info from next_tupel'
        write(lulog,*) '----------------------'
        write(lulog,*) ' nidx = ',nidx
      end if

      ! set offset array
      idx = 0
      do ihpv = 1, ngastp
        idxoff(ihpv)=idx
        idx = idx+iocc(ihpv)
      end do

      if (first) then

        next_tupel = .true.
        if (nidx.eq.0) return

        ihpv_loop: do ihpvdx = 1, ngastp
          ! retrieve actual index to be incremented next
          ihpv = ihpvseq(ihpvdx)
          nel = iocc(ihpv)
          if (nel.le.0) cycle
          ! index range within tupel
          idxst = idxoff(ihpv)+1
          idxnd = idxoff(ihpv)+nel
          igr = igrph(ihpv)     ! number of graph
          ms  = msdst(ihpv)     ! ms of string
          igam = igamdst(ihpv)  ! IRREP of string
          ! first subspace distribution
          idss(idxst:idxnd) = 1
          dss_loop: do
            ! lexically lowest string
            ! reform distribution to occupation
            ioss(1:ngas_hpv(ihpv)) = 0
            do idx = idxst, idxnd
              ioss(idss(idx)) = ioss(idss(idx))+1
            end do

            succ = lexlstr(nel,ms,
     &           ioss,idorb(idxst:idxnd),
     &           idspn(idxst:idxnd),
     &           mostnd(1,1,idx_gas(ihpv)),nsym,ngas_hpv(ihpv))

            if (succ) then
              str_loop: do
                ! check symmetry
                succ = igamstr2(nel,idorb(idxst),igamorb).eq.igam
                ! exit if successful
                if (succ) exit dss_loop
                ! else: get next possible string
                succ = nextstr(nel,ms,
     &               idss(idxst:idxnd),idorb(idxst:idxnd),
     &               idspn(idxst:idxnd),
     &               mostnd(1,1,idx_gas(ihpv)),nsym,ngas_hpv(ihpv)) 

                if (.not.succ) exit str_loop
              end do str_loop
            end if
            ! if no appropriate string found: 
            !  get next subspace distr.
            succ = next_ssd(idss(idxst:idxnd),nel,nel,
     &               ngas_hpv(ihpv),igas_restr(1,1,1,igr))
            ! nothing: we have to give up
            if (.not.succ) exit dss_loop
          end do dss_loop

          ! no success for current ihpv => give up
          if (.not.succ) exit ihpv_loop

        end do ihpv_loop

        ! remember for next passes
        lexlscr(1:nidx,1) = idorb(1:nidx)
        lexlscr(1:nidx,2) = idspn(1:nidx)
        lexlscr(1:nidx,3) = idss(1:nidx)
        
        next_tupel = succ

        if (ntest.ge.100) then
          if (succ) then
            write(lulog,*) 'first string: ',idorb(1:nidx)
            write(lulog,*) '              ',idspn(1:nidx)
            write(lulog,*) '              ',idss(1:nidx)
          else
            write(lulog,*) 'no string exists'
          end if
        end if

        return
      else

        next_tupel = .false.
        if (nidx.eq.0) return

        if (ntest.ge.100) then
          write(lulog,*) ' input string:',idorb(1:nidx)
          write(lulog,*) '              ',idspn(1:nidx)
          write(lulog,*) '              ',idss(1:nidx)
        end if

        ihpv_loop2: do ihpvdx = 1, ngastp
          ! retrieve actual index to be incremented next
          ihpv = ihpvseq(ihpvdx)
          nel = iocc(ihpv)
          if (nel.le.0) cycle
          ! index range within tupel
          idxst = idxoff(ihpv)+1
          idxnd = idxoff(ihpv)+nel
          igr = igrph(ihpv)     ! number of graph
          ms  = msdst(ihpv)     ! ms of string
          igam = igamdst(ihpv)  ! IRREP of string

          dss_loop2: do
            str_loop2: do
              ! try next string ...
              succ = nextstr(nel,ms,
     &           idss(idxst:idxnd),idorb(idxst:idxnd),
     &           idspn(idxst:idxnd),
     &           mostnd(1,1,idx_gas(ihpv)),nsym,ngas_hpv(ihpv)) 

              ! ... no further string ? ....
              if (.not.succ) exit str_loop2

              ! ... else check symmetry ...
              if (succ)
     &             succ = igamstr2(nel,idorb(idxst),igamorb).eq.igam

              ! if symmetry was alright, we are done
              if (succ) exit dss_loop2
            end do str_loop2
                
            ! here we land, if we need to try the next
            ! subspace distribution
            succ = next_ssd(idss(idxst:idxnd),nel,nel,
     &             ngas_hpv(ihpv),igas_restr(1,1,1,igr))
            ! nothing: we have to give up
            if (.not.succ) exit dss_loop2

            ! else we have to get the lexically lowest string
            ! for the new subspace distribution

            ! reform distribution to occupation
            ioss(1:ngas_hpv(ihpv)) = 0
            do idx = idxst, idxnd
              ioss(idss(idx)) = ioss(idss(idx))+1
            end do

            succ = lexlstr(nel,ms,
     &             ioss,idorb(idxst:idxnd),
     &             idspn(idxst:idxnd),
     &             mostnd(1,1,idx_gas(ihpv)),nsym,ngas_hpv(ihpv))

            ! ... check symmetry ...
            if (succ)
     &           succ = igamstr2(nel,idorb(idxst),igamorb).eq.igam

            ! ... we are done?
            if (succ) exit dss_loop2

            ! ... else go to top of loop and increment string

          end do dss_loop2

          ! success for current ica,ihpv => exit, we are done
          if (succ) exit ihpv_loop2

          ! ... else, reset current string to lexically
          ! lowest (as stored in lexlscr) and increment next
          ! hpv block
          idorb(idxst:idxnd) = lexlscr(idxst:idxnd,1)
          idspn(idxst:idxnd) = lexlscr(idxst:idxnd,2)
          idss(idxst:idxnd)  = lexlscr(idxst:idxnd,3)

        end do ihpv_loop2
        
        next_tupel = succ

        if (ntest.ge.100) then
          if (succ) then
            write(lulog,*) ' next string: ',idorb(1:nidx)
            write(lulog,*) '              ',idspn(1:nidx)
            write(lulog,*) '              ',idss(1:nidx)
          else
            write(lulog,*) ' no further string exists'
          end if
        end if

        return

      end if

      end
*----------------------------------------------------------------------*

