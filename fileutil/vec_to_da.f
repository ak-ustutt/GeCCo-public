*----------------------------------------------------------------------*
      subroutine vec_to_da(ffda,idxvec,vec,len)
*----------------------------------------------------------------------*
*     punch a vector of length len to disc
*     ffda is the file handle
*     idxvec is the number of the vector, in case of several vectors
*     in the same file; the vectors are always aligned such that they
*     start at the beginning of a block
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      integer ::
     &     ntest = 00

      type(filinf), intent(in) ::
     &     ffda
      integer, intent(in) ::
     &     len, idxvec
      real(8), intent(in) ::
     &     vec(len)

      integer ::
     &     lblk, luda,
     &     nrecs, len_rest, irec, irecst, irecnd, idxst, idxnd, len_wr
      real(8) ::
     &     xnrm

      real(8), external ::
     &     ddot

      luda = ffda%unit
      lblk = ffda%reclen

      if (luda.lt.0)
     &     call quit(1,'vec_to_da',
     &     'file is not open: '//trim(ffda%name))

      nrecs = len/lblk
      len_rest = mod(len,lblk)
      if (len_rest.gt.0) nrecs = nrecs+1
      if (len_rest.eq.0) len_rest=lblk

      irecst = (idxvec-1)*nrecs+1
      irecnd = (idxvec-1)*nrecs+nrecs

      len_wr = lblk
      idxst = 1
      do irec = irecst, irecnd
        if (irec.eq.irecnd) len_wr = len_rest
        idxnd = idxst-1 + len_wr
        write(luda,rec=irec) vec(idxst:idxnd)
        idxst = idxnd+1
      end do

      if (ntest.eq.100) then
        write(lulog,*) 'wrote ',nrecs,' blocks'
        write(lulog,*) 'length of vector: ',len
        xnrm = sqrt(ddot(len,vec,1,vec,1))
        write(lulog,*) 'norm of vector: ',xnrm
      end if

      return

      end
