*----------------------------------------------------------------------*
      subroutine fh_init(iprint)
*----------------------------------------------------------------------*
!> initializes the freeunit mechanic and the blocklength parameter
!! 
*----------------------------------------------------------------------*

      implicit none
      include "freeunits.h"
      include "stdunit.h"
      include "ioparam.h"
      
      integer, intent(in) ::
     &     iprint             ! print level

      integer :: iunit
      logical :: lopen

! set standard values for direct access
      nrecfc  = 8               ! can be compiler dependent --> improve that
      
      ! currently 32K blocks
      ! could be adapted to size of current calculation
      !lblk_da is record_length of a direct access file ( reals)
      lblk_da = 32*1024/nrecfc

      ifrunit(1:6)       = 3 ! reserved
      ifrunit(7:mxpunit) = 0
      do iunit = 7, mxpunit
        if (ifrunit(iunit).eq.0) then
          inquire(unit=iunit,opened=lopen)
          if (lopen) then
            ifrunit(iunit) = -2
          end if
        end if
      end do

!      if (iprint.gt.10) then
!        call fh_prstat(lulog)
!      end if

      end

*----------------------------------------------------------------------*
      subroutine fh_prstat(lulog)
*----------------------------------------------------------------------*
      implicit none
      include "freeunits.h"

      integer, intent(in) ::
     &     lulog         ! output unit

      integer ::
     &     iunit
      logical ::
     &     lopen
      character ::
     &     cstat*20, name*256

      write(lulog,'(2(/x,a))')
     &       'Status of unit numbers',
     &       '======================'
      do iunit = 1, mxpunit
        if (ifrunit(iunit).eq.0)  cstat = 'FREE                '
        if (ifrunit(iunit).eq.3)  cstat = 'RESERVED BY SYSTEM  '
        if (ifrunit(iunit).eq.2)  cstat = 'RESERVED BY USER    '
        if (ifrunit(iunit).eq.-2) cstat = 'OCCUPIED ON INIT    '
        if (ifrunit(iunit).eq.1)  cstat = 'RESERVED BY IGETUN  '
        if (ifrunit(iunit).eq.-1) cstat = 'OCCUPIED AFTER INIT '
        name = ''
        if (ifrunit(iunit).eq.1) then
          inquire(unit=iunit,name=name,opened=lopen)
        end if
        write(lulog,'(x,i3,2x,a,2x,l,2x,a)')
     &        iunit, cstat, lopen, trim(name)
      end do

      end
*----------------------------------------------------------------------*
      integer function igetunit(iwish)
*----------------------------------------------------------------------*
      
      implicit none
      include "freeunits.h"
      include "stdunit.h"

      integer, intent(in) ::
     &     iwish

      integer :: iunit
      logical :: lopen

      ! trying to fulfil your wishes (if you asked for that)
      if (iwish.ge.7) then
        iunit = iwish
        if (ifrunit(iwish).eq.0) then
          inquire(unit=iunit,opened=lopen)
          if (lopen) then
            ifrunit(iunit) = -1
          else
            ifrunit(iunit) = 1
            igetunit = iunit
            return
          end if
        end if
      end if

      ! looking for the lowest free unit
      do iunit = 7, mxpunit
        if (ifrunit(iunit).eq.0) then
          inquire(unit=iunit,opened=lopen)
          if (lopen) then
            ifrunit(iunit) = -1
          else
            ifrunit(iunit) = 1
            igetunit = iunit
            return
          end if
        end if
      end do

      write(lulog,*) 'Fatal: run out of unit numbers in igetun'
      call fh_prstat(6)
      call quit(1,'igetunit','run out of unit numbers')

      end
*----------------------------------------------------------------------*
      subroutine relunit(iunit,cstat)
*----------------------------------------------------------------------*

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      integer, intent(in) :: iunit
      character, intent(in) :: cstat*(*)

      logical :: lopen, lnamed
      integer :: len
      character :: filnam*256, cscr*8

      if (iunit.le.mxpunit.and.iunit.gt.0.and.ifrunit(iunit).eq.1) then
        inquire(unit=iunit,opened=lopen)
        if (lopen) then
          ! check cstat, as a forgotten or wrong cstat can lead to
          ! strange problems (e.g. my laptop stopped [!?!])
          len = len_trim(cstat)
          if (len.le.6) then 
            cscr(1:len) = cstat(1:len)
            call uppcas(cscr)            
            if (cscr(1:4).ne.'KEEP'.and.cscr(1:6).ne.'DELETE')
     &           len = -len
          else
            len = -len
          end if
          if (len.le.0) then
            len = min(40,abs(len))
            write(lulog,*) ' illegal status line: "',cstat(1:len),'"'
            call quit(1,'relunit','illegal status line')
          end if
          close(iunit,status=cscr(1:len))
        end if
        ifrunit(iunit) = 0
      else
        write(lulog,*) 'relunit called for illegal unit ',
     &       iunit, ifrunit(iunit)
        inquire(unit=iunit,opened=lopen,named=lnamed,name=filnam)
        if (lopen) then
          write(lulog,*) ' the file was open'
        else
          write(lulog,*) ' the file was not open'
        end if
        if (lnamed) then
          len = len_trim(filnam)
          write(lulog,*) ' the file was named: "',filnam(1:len),'"'
        else
          write(lulog,*) ' the file was not named'
        end if
        call quit(1,'relunit','illegal unit')
      end if
      
      end

*----------------------------------------------------------------------*
      integer function iopen_uus()
*----------------------------------------------------------------------*
* open  unnamed unformatted sequential

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      integer :: iunit
      integer, external :: igetunit

      iunit = igetunit(-1)
      open(unit=iunit,form="unformatted",access="sequential")

      iopen_uus = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nus(filename)
*----------------------------------------------------------------------*
* open  named unformatted sequential

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      character, intent(in) :: filename*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)

      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="sequential")

      iopen_nus = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nfs(filename)
*----------------------------------------------------------------------*
* open  named formatted sequential

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      character, intent(in) :: filename*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="formatted",access="sequential")

      iopen_nfs = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nustr(filename)
*----------------------------------------------------------------------*
* open  named unformatted stream

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      character, intent(in) :: filename*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="stream")

      iopen_nustr = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nuss(filename,statstr)
*----------------------------------------------------------------------*
* open  named unformatted sequential in status 'statstr'

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      character, intent(in) :: filename*(*), statstr*(*)

      integer :: iunit, len
      integer, external :: igetunit

      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="sequential",status=statstr)

      iopen_nuss = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nud(filename,lrec)
*----------------------------------------------------------------------*
* open  named unformatted direct access

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      character, intent(in) :: filename*(*)
      integer, intent(in)   :: lrec

      integer :: iunit, len
      integer, external :: igetunit, ngtrecfc


      iunit = igetunit(-1)
      len = len_trim(filename)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="direct",recl=lrec*ngtrecfc())

      iopen_nud = iunit

      end

*----------------------------------------------------------------------*
      integer function iopen_nuds(filename,lrec,statstr)
*----------------------------------------------------------------------*
* open  named unformatted direct access

      implicit none
      include "freeunits.h"
      include "stdunit.h"

      character, intent(in) :: filename*(*), statstr*(*)
      integer, intent(in)   :: lrec

      integer :: iunit, len, len2
      integer, external :: igetunit, ngtrecfc


      iunit = igetunit(-1)
      len  = len_trim(filename)
      len2 = len_trim(statstr)
      open(unit=iunit,file=filename(1:len),
     &     form="unformatted",access="direct",
     &     status=statstr(1:len2),recl=lrec*ngtrecfc())

      iopen_nuds = iunit

      end

*----------------------------------------------------------------------*
      subroutine unit_info(lu)
*----------------------------------------------------------------------*
      implicit none
      include "freeunits.h"
      include "stdunit.h"
      
      integer, intent(in) :: lu

      logical :: lnamed, lopen
      integer :: len, lrec, nrec
      character :: str*256

      inquire(unit=lu,named=lnamed,opened=lopen,name=str)

      write(lulog,'(/x,78("="))')
      write(lulog,*) ' Information on unit ',lu
      if (lopen) then
        write(lulog,*) '  file is open'
        if (lnamed) then
          len = len_trim(str)
          write(lulog,*) '  filename: ',str(1:len)
        else
          write(lulog,*) '  unnamed'
        end if
        inquire(unit=lu,form=str)
        len = len_trim(str)
        write(lulog,*) '  format: ',str(1:len)
        inquire(unit=lu,access=str)
        len = len_trim(str)
        write(lulog,*) '  access: ',str(1:len)
        if (str(1:len).eq.'DIRECT'.or.str(1:len).eq.'direct') then
          inquire(unit=lu,recl=lrec,nextrec=nrec)
          write(lulog,*) '  record-length: ',lrec
          write(lulog,*) '  next record  : ',nrec
        end if
      else
        write(lulog,*) ' file is not open'
      end if
      write(lulog,'(x,78("="))')

      end

*----------------------------------------------------------------------*
      integer function ngtrecfc()
*----------------------------------------------------------------------*
*     little wrapper to access ioparam.h
*----------------------------------------------------------------------*

      implicit none

      include 'ioparam.h'

      ngtrecfc = nrecfc

      end
