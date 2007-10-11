*----------------------------------------------------------------------*
      integer function sign_hpvx(mode,iocc1,dag1,iocc2,dag2)
*----------------------------------------------------------------------*
*     return the sign originating from bringing of H/P/V/X blocks
*     of the strings into correct order when merging the primitive
*     vertices with occupations iocc1 and iocc2     
*      mode==1:  1C2C 2A1A merge
*      mode==2:  1C2C 1A2A merge
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     mode,
     &     iocc1(ngastp,2), iocc2(ngastp,2)
      logical, intent(in) ::
     &     dag1,dag2

      logical ::
     &     inv
      integer ::
     &     ica1, ica2, hpvx, hpvx2, sign0

      inv = dag1.and..not.dag2 .or. .not.dag1.and.dag2

      sign0 = 0

      select case(mode)
      case(1)

        ! ----------------------------------------------------------- !
        ! sign for bringing both C and A                              !
        ! from H1,P1,V1,X1;H2,P2,V2,X2  to  H1H2,P1P2,V1V2,X1X2       !
        ! (A is stored in reverse sequence, therefore this reordering !
        ! is consistent with formal C1,C2;A2,A1 merge)                !
        ! ----------------------------------------------------------- !
        do ica1 = 1, 2
          ica2 = ica1
          if (inv) ica2 = 3-ica1
          do hpvx = 2, ngastp
            if (iocc1(hpvx,ica1).gt.0) then
              do hpvx2 = 1, hpvx-1
                sign0 = mod(sign0 
     &               +iocc1(hpvx,ica1)*iocc2(hpvx2,ica2),2)
              end do
            end if
          end do
        end do

      case(2)

        ! ----------------------------------------------------------- !
        ! sign for bringing C                                         !
        ! from H1,P1,V1,X1;H2,P2,V2,X2  to  H1H2,P1P2,V1V2,X1X2       !
        ! ----------------------------------------------------------- !
        ica2 = 1
        if (inv) ica2 = 2
        do hpvx = 2, ngastp
          if (iocc1(hpvx,1).gt.0) then
            do hpvx2 = 1, hpvx-1
              sign0 = mod(sign0 
     &             +iocc1(hpvx,1)*iocc2(hpvx2,ica2),2)
            end do
          end if
        ! ----------------------------------------------------------- !
        ! ... and for bringing A                                      !
        ! from H2,P2,V2,X2;H1,P1,V1,X1  to  H2H1,P2P1,V2V1,X2X1       !
        ! (A is stored in reverse sequence, therefore this reordering !
        ! is consistent with formal C1,C2;A1,A2 merge)                !
        ! ----------------------------------------------------------- !
          if (iocc2(hpvx,2).gt.0) then
            do hpvx2 = 1, hpvx-1
              sign0 = mod(sign0 
     &             +iocc2(hpvx,3-ica2)*iocc1(hpvx2,2),2)
            end do
          end if
        end do

      case default

        call quit(1,'sign_hpvx','call with unknown mode-parameter')

      end select

      sign_hpvx = 1
      if (sign0.eq.1) sign_hpvx = -1

      return
      end
