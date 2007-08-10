*----------------------------------------------------------------------*
      subroutine set_idx_seq(idxseq,iocc,njoined,sep_ca,hpvxseq)
*----------------------------------------------------------------------*
*     for a given (multi-component) occupation set a similar array
*     which instead of occupations contains the sequence number of
*     the corresponding string, e.g.
*
*     / 1 0 0 1 \/2 0 0 0\   ->  /4 0 0 1\/5 0 0 0\
*     \ 0 2 0 0 /\0 1 0 0/       \0 2 0 0/\0 3 0 0/
*
*     if sep_ca==.false., the behaviour is as in the above example
*     if sep_ca==.true.,  the behaviour is instead (C/A separate):
*
*     / 1 0 0 1 \/2 0 0 0\   ->  /2 0 0 1\/3 0 0 0\
*     \ 0 2 0 0 /\0 1 0 0/       \0 1 0 0/\0 2 0 0/
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined),
     &     hpvxseq(ngastp)
      logical, intent(in) ::
     &     sep_ca
      integer, intent(out) ::
     &     idxseq(ngastp,2,njoined)

      integer ::
     &     idx, idx_hpvx, ica, ijoin, hpvx

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_idx_seq')
        write(luout,*) 'sep_ca = ',sep_ca
        write(luout,*) 'input occupation:'
        call wrt_occ_n(luout,iocc,njoined)
      end if

      idxseq = 0

      if (sep_ca) then

        idx = 0
        do idx_hpvx = 1, ngastp
          hpvx = hpvxseq(idx_hpvx)
          do ijoin = 1, njoined
            if (iocc(hpvx,1,ijoin).eq.0) cycle
            idx = idx+1
            idxseq(hpvx,1,ijoin) = idx
          end do
        end do

        idx = 0
        do idx_hpvx = 1, ngastp
          hpvx = hpvxseq(idx_hpvx)
          do ijoin = 1, njoined
            if (iocc(hpvx,2,ijoin).eq.0) cycle
            idx = idx+1
            idxseq(hpvx,2,ijoin) = idx
          end do
        end do

      else

        idx = 0
        do idx_hpvx = 1, ngastp
          hpvx = hpvxseq(idx_hpvx)
          do ica = 1, 2
            do ijoin = 1, njoined
              if (iocc(hpvx,ica,ijoin).eq.0) cycle
              idx = idx+1
              idxseq(hpvx,ica,ijoin) = idx
            end do
          end do
        end do
        
      end if

      if (ntest.ge.100) then
        write(luout,*) 'generated array::'
        call wrt_occ_n(luout,idxseq,njoined)
      end if

      return
      end
