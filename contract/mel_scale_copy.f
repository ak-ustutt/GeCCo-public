!-----------------------------------------------------------------------
!> Subroutine to scale matrix element lists
!!
!!This suroutine has multiple uses depending on mode. (parameter set in par_scale_copy_modes.h)
!!+ MOD_SCALE:  scales all matrix elements by fac.
!!+ MOD_MULT: multiplies the me_res elementwise with me_inp (requires 2 buffer)
!!+ MOD_SQUARE: squares all matrix elements and multiplies them by fac.
!!+ MOD_PRECOND: divides label_res elementwise by label_inp and applies fac(1). (requires 2 buffer)
!!+ MOD_PRC_THRESH: sets all elements to at least fac.
!!  This subroutine also applies a sign fix to the action.
!!
!!\param[inout] me_inp input list 
!!\param[inout] me_res result list ()
!!\param[inout] buf1, buf2  buffer of length lbuf (min largest section with specific sign)
!!\param[in] lbuf length of buffers
!!\param[in] nbuf 1 if only buf1 is allocated, 2 otherwise
!!\param[in] mode one of MOD_SCALE MOD_MULT  MOD_SQUARE MOD_PRECOND MOD_PRC_THRESH
!!\param[in] fac list of factors which are applied (repeatedly)
!!\param[in] nfac length of factor list
!!\param[in] opti_info an opti info object for sections with sign change
!-----------------------------------------------------------------------
      subroutine mel_scale_copy(
     &     me_inp, me_res,
     &     buf1,buf2, lbuf, nbuf,
     &     fac, nfac,
     &     mode,
     &     opti_info)
!-----------------------------------------------------------------------
      
      implicit none
      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_optimize_info.h'
      include 'par_scale_copy_modes.h'
      
      !using includes of parent
      integer, parameter ::
     &     ntest = 00
      character(len=*), parameter ::
     &     i_am = 'mel_scale_copy'
      integer,intent(in)::
     &     lbuf, nbuf, nfac,
     &     mode
      type(me_list), intent(inout)::
     &     me_inp, me_res
      
      type(optimize_info),intent(in)::
     &     opti_info
      real(8), intent(inout) ::
     &     buf1(:), buf2(:)
      real(8),intent(in)::
     &     fac(nfac)
      
      integer::
     &     len_op,
     &     idisc_off_src, idisc_off_tgt,
     &     idxst_src,idxst_tgt,
     &     idxnd_src,idxnd_tgt,
     &     nsec, isec,
     &     ifac,
     &     idx,
     &     smapre_num, negpre_num
      real(8),pointer::
     &     signsec(:)
      integer, pointer ::
     &     lensec(:), idstsec(:)
      type(filinf),pointer::
     &     ffop_src,ffop_tgt
      

      if ((mode .eq. MOD_MULT.or. mode .eq. MOD_PRECOND )
     &     .and. nbuf .lt. 2 )
     &     call quit(2,i_am,"this mode requires two buffer")
      
      if (ntest.ge.10)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      if(ntest.ge.20) then
         write(lulog,*)me_res%label
         write(lulog,*)me_inp%label
      end if

      ffop_tgt => me_res%fhand
      ffop_src => me_inp%fhand
      smapre_num=0
      negpre_num=0 
      ! offset on file (if more than one instance of list exists)
      idisc_off_src =
     &     ffop_src%length_of_record*(ffop_src%current_record-1)
      idisc_off_tgt =
     &     ffop_tgt%length_of_record*(ffop_tgt%current_record-1)

      nsec = opti_info%nsec(1)
      lensec => opti_info%nwfpsec(1:nsec)
      idstsec => opti_info%idstsec(1:nsec)
      signsec => opti_info%signsec(1:nsec)


      do isec = 1, nsec
         len_op=lensec(isec)
         ifac=1
         idxst_src = idisc_off_src+idstsec(isec)
         idxst_tgt = idisc_off_tgt+idstsec(isec)
         
         do while(idxst_src.le.idisc_off_src+idstsec(isec)-1+len_op)
            idxnd_src = min(idisc_off_src+idstsec(isec)-1+len_op,
     &           idxst_src-1+lbuf)
            idxnd_tgt = min(idisc_off_tgt+idstsec(isec)-1+len_op,
     &           idxst_tgt-1+lbuf)
            
            call get_vec(ffop_src,buf1,idxst_src,idxnd_src)
            if (mode .eq. MOD_MULT
     &           .or. mode .eq. MOD_PRECOND
     &           .or. mode .eq. MOD_ADD_VEC)
     &           call get_vec(ffop_tgt,buf2,idxst_tgt,idxnd_tgt)

            if (mode .eq. MOD_PRECOND)then
               call diavc(buf1,buf2,
     &              signsec(isec)*fac(1),buf1,
     &              0d0,idxnd_src-idxst_src+1)
            else if (nfac .eq.1)then
!     optimization:
!     nfac=1 is most common case and can effectively be vectorized
               select case(mode)
            case(MOD_SQUARE)
               do idx = 1, idxnd_src-idxst_src+1
                  buf1(idx) = signsec(isec)*fac(1)*(buf1(idx)**2)
               end do
            case(MOD_ADD_VEC)
               do idx = 1, idxnd_src-idxst_src+1
                  buf1(idx) = buf1(idx)+signsec(isec)*fac(1)*buf2(idx)
               end do
            case(MOD_MULT)
               do idx = 1, idxnd_src-idxst_src+1
                  buf1(idx) = signsec(isec)*
     &                 fac(1)*buf1(idx)*buf2(idx)
               end do
            case(MOD_PRC_THRESH)
               do idx = 1, idxnd_src-idxst_src+1
                  if (buf1(idx).lt.fac(1))smapre_num=smapre_num+1
                  if (buf1(idx).lt.0d0) negpre_num = negpre_num + 1
                  buf1(idx) = max(fac(1),buf1(idx))
               end do
            case(MOD_SCALE)
               do idx = 1, idxnd_src-idxst_src+1
                  buf1(idx) = signsec(isec)*fac(1)*buf1(idx)
               end do
            case default
               call quit(2,i_am,"unknown mode")
               end select 
!     could also be optimized for large or small nfac (no application yet)
            else                !nfac != 1
               do idx = 1, idxnd_src-idxst_src+1
                  select case(mode)
               case(MOD_SQUARE)
                  buf1(idx) = signsec(isec)*fac(ifac)*(buf1(idx)**2)
               case(MOD_ADD_VEC)
                     buf1(idx) = buf1(idx)+signsec(isec)*fac(ifac)*buf2(idx)
               case(MOD_MULT)
                  buf1(idx) = signsec(isec)*
     &                 fac(ifac)*buf1(idx)*buf2(idx)
               case(MOD_PRC_THRESH)
                  if (buf1(idx).lt.fac(ifac))smapre_num=smapre_num+1
                  if (buf1(idx).lt.0d0) negpre_num = negpre_num + 1
                  buf1(idx) = max(fac(ifac),buf1(idx))
               case(MOD_SCALE)
                  buf1(idx) = signsec(isec)*fac(ifac)*buf1(idx)
                  
               case default
                  call quit(2,i_am,"unknown mode")
               end select 
               ifac = ifac + 1
               if (ifac.gt.nfac) ifac = 1
            end do
            end if
            call put_vec(ffop_tgt,buf1,idxst_src,idxnd_src)
            idxst_src = idxnd_src+1
            idxst_tgt = idxnd_tgt+1
            
            if (mode.eq.MOD_PRC_THRESH) then
               if (smapre_num.gt.0) then
                  write(lulog,'(1x,a,i9,a,g9.2)')
     &                 'number of small preconditioner elements: ',
     &                 smapre_num, '; set to ',fac(1:nfac)
            end if
            if (negpre_num.gt.0) then
               write(lulog,'(1x,a,i9)')
     &              'number of negative preconditioner elements: ',
     &              negpre_num
               call warn(i_am,
     &                 'negative preconditioner elements!')
            end if
         end if
      end do
      end do
      return
      end subroutine
