* declare interfaces to functions in baserout library
      interface
        integer function ielsum(ivec,nelmnt)
        integer, intent(in) ::
     &       nelmnt,ivec(nelmnt)
        end function ielsum

        integer(4) function i4elsum(ivec,nelmnt)
        integer(4), intent(in) ::
     &       nelmnt,ivec(nelmnt)
        end function i4elsum

        integer(8) function i8elsum(ivec,nelmnt)
        integer(8), intent(in) ::
     &       nelmnt,ivec(nelmnt)
        end function i8elsum

        integer function ielsqsum(ivec,nelmnt)
        integer, intent(in) ::
     &       nelmnt,ivec(nelmnt)
        end function ielsqsum

        integer(4) function i4elsqsum(ivec,nelmnt)
        integer(4), intent(in) ::
     &       nelmnt,ivec(nelmnt)
        end function i4elsqsum

        integer function ielprd(ivec,nel)
        integer, intent(in) ::
     &       nel,ivec(nel)
        end function ielprd

        integer function ifac(n)
        integer, intent(in) ::
     &       n
        end function ifac

        real(8) function xfac(n)
        integer, intent(in) ::
     &       n
        end function xfac

        integer function first_nonblank(str)
        implicit none 
        character, intent(in) ::
     &     str*(*)
        end function first_nonblank

        integer function zirat()
        end function zirat

        logical function cmpxarr(xval,arr,nel)
        implicit none
        integer, intent(in) ::
     &     nel
        real(8), intent(in) ::
     &     arr(nel), xval
        end function

        logical function cmpiarr(ival,iarr,nel)
        implicit none
        integer, intent(in) ::
     &     nel,
     &     iarr(nel), ival
        end function

        logical function next_part_number(init,restrict,ipart,
     &     inum,nsum,inummin,inummax)
        implicit none
        logical, intent(in) ::
     &       init, restrict
        integer, intent(in) ::
     &       inum, nsum, inummin, inummax
        integer, intent(inout) ::
     &       ipart(nsum)
        end function

        logical function next_part_pair(init,restrict,ipart,
     &     inum,nsum,inummin,inummax)
        implicit none
        logical, intent(in) ::
     &       init, restrict
        integer, intent(in) ::
     &       inum(2), nsum, inummin, inummax
        integer, intent(inout) ::
     &       ipart(2,nsum)
        end function

        logical function list_cmp(list1,list2,nel)
        implicit none
        integer, intent(in) ::
     &     nel, list1(nel), list2(nel)        
        end function

        integer function idxlist(inum,ilist,nel,inc)
        implicit none
        integer, intent(in) ::
     &     inum, nel, inc, ilist(nel)
        end function

        integer function imltlist(inum,ilist,nel,inc)
        implicit none
        integer, intent(in) ::
     &       inum, nel, inc, ilist(nel)
        end function

        logical function next_dist(idist,nel,imnmx,inc)
        implicit none
        integer, intent(in) ::
     &     inc, nel, imnmx(2,nel)
        integer, intent(inout) ::
     &     idist(nel)
        end function

        logical function next_dist2(idist,nel,imn,imx,inc)
        implicit none
        integer, intent(in) ::
     &     inc, nel, imn(nel), imx(nel)
        integer, intent(inout) ::
     &     idist(nel)
        end function

      end interface
