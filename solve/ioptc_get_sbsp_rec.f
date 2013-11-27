*----------------------------------------------------------------------*
      integer function ioptc_get_sbsp_rec(ivec,iord,ndim,maxdim)
*----------------------------------------------------------------------*
*
*     get record number for desired record, and update iord and ndim,
*     if necessary
*
*     ivec = -1,-2,...,-ndim : get record of last, 2nd-last, etc. vect.
*     ivec =  1, 2,..., ndim : get record of first, 2nd, etc. vect.
*
*     ivec = 0   add new vector
*
*     iord(irec) contains the number of the vector in the subspace
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ivec, maxdim
      integer, intent(inout) ::
     &     ndim, iord(maxdim)

      integer ::
     &     irec, itarget

      if (ntest.ge.10) then
        write(lulog,*) '------------------------------'
        write(lulog,*) ' info from ioptc_get_sbsp_rec'
        write(lulog,*) '------------------------------'
        write(lulog,*) ' ivec = ',ivec
        write(lulog,*) ' ndim, maxdim: ', ndim, maxdim
        write(lulog,*) ' iord = ',iord(1:maxdim)
      end if
      
      if (ivec.eq.0) then

        if (ndim.lt.maxdim) then
          ndim = ndim+1
          iord(ndim) = ndim
          ioptc_get_sbsp_rec = ndim
        else
          do irec = 1, ndim
            iord(irec) = iord(irec)-1
            if (iord(irec).eq.0) then
              iord(irec) = ndim
              ioptc_get_sbsp_rec = irec
            end if
          end do
        end if

      else
        if (ivec.lt.0.and.ivec.ge.-ndim) then
          itarget = ndim-ivec+1
        else if (ivec.gt.0.and.ivec.le.ndim) then
          itarget = ivec
        else
          write(lulog,*) 'invalid ivec passed to ioptc_get_sbsp_rec: '
          write(lulog,*) ' ndim, maxdim: ',ndim,maxdim
          write(lulog,*) ' ivec:         ',ivec
          call quit(1,'ioptc_get_sbsp_rec',
     &         'invalid ivec passed to function')
        end if

        do irec = 1, ndim
          if (iord(irec).eq.itarget) then
            ioptc_get_sbsp_rec = irec
          end if
        end do

      end if

      if (ntest.ge.10) then
        write(lulog,*) 'function returns ',ioptc_get_sbsp_rec
        if (ivec.eq.0) then
          write(lulog,*) 'Updated subspace:'
          write(lulog,*) ' ndim: ', ndim
          write(lulog,*) ' iord = ',iord(1:maxdim)
        end if

      end if

      return

      end
