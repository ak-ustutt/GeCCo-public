************************************************************************
!>alternative version to read a vector from file
!!
!!Reads a vector from file starting at irecst
!!If lenvec is not a multiple of the record length
!!only lenvec%reclen entries of the last record a read.
!!\param[in] ffvec Filhandler
!!\param[in] lenvec length of buffer
!!\param[in] irecst record, where reading should start.
!!\param[out] vec output vector
      subroutine get_vec_a(ffvec,vec,lenvec,irecst)

************************************************************************
*     alternative version of get_vec
************************************************************************
      
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      type(filinf), intent(in) ::
     &     ffvec
      integer, intent(in) ::
     &     lenvec, irecst
      real(8), intent(out) ::
     &     vec(lenvec)

      integer ::
     &     lblk, luvec, nrec, len_last, len, irec, irecnd, idxst, idxnd

      lblk = ffvec%reclen
      luvec = ffvec%unit
      nrec = lenvec/lblk
      len_last = mod(lenvec,lblk)
      if (len_last.gt.0) nrec = nrec+1
      if (len_last.eq.0) len_last = lblk
      len = lblk
      idxst = 1
      irecnd = irecst-1 + nrec
      do irec = irecst, irecnd
        if (irec.eq.irecnd) len = len_last
        idxnd = idxst-1 + len
        read(luvec,rec=irec) vec(idxst:idxnd)
        idxst = idxnd + 1
      end do

      return
      end
