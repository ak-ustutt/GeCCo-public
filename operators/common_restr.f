      logical function common_restr(rst_res,rst1,dag1,rst2,dag2,
     &                        hpvxgas,ngas,nspin)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00
      character(len=12), parameter ::
     &     i_am = 'common_restr'

      integer, intent(in) ::
     &     ngas, nspin, hpvxgas(ngas,nspin),
     &     rst1(2,ngas,2,2,nspin),
     &     rst2(2,ngas,2,2,nspin)
      logical, intent(in) ::
     &     dag1, dag2
      integer, intent(out) ::
     &     rst_res(2,ngas,2,2,nspin)

      logical ::
     &     ok
      integer ::
     &     ispin, ica, ica1, ica2

      logical, external ::
     &     is_proper_restr_for_hpvx, zero_ivec

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_func,i_am)
        write(luout,*) 'on entry:'
        call wrt_rstr(luout,rst1,ngas)
        call wrt_rstr(luout,rst2,ngas)
      end if

      if (nspin.gt.1) call quit(1,i_am,'nspin.gt.1 -> adapt me fully')

      rst_res = 0
      ok = .true.
      do ispin = 1, nspin
        ! mask = 1: paths included in occupation graph
        do ica = 1, 2
          ica1 = ica
          if (dag1) ica1 = 3-ica
          ica2 = ica
          if (dag2) ica2 = 3-ica
          call common_restr_for_hpvx(rst_res(:,:,ica ,1,ispin),
     &                                  rst1(:,:,ica1,1,ispin),
     &                                  rst2(:,:,ica2,1,ispin),
     &                                  IHOLE,hpvxgas(:,ispin),ngas)
          call common_restr_for_hpvx(rst_res(:,:,ica ,1,ispin),
     &                                  rst1(:,:,ica1,1,ispin),
     &                                  rst2(:,:,ica2,1,ispin),
     &                                  IPART,hpvxgas(:,ispin),ngas)
          call common_restr_for_hpvx(rst_res(:,:,ica ,1,ispin),
     &                                  rst1(:,:,ica1,1,ispin),
     &                                  rst2(:,:,ica2,1,ispin),
     &                                  IVALE,hpvxgas(:,ispin),ngas)
          call common_restr_for_hpvx(rst_res(:,:,ica ,1,ispin),
     &                                  rst1(:,:,ica1,1,ispin),
     &                                  rst2(:,:,ica2,1,ispin),
     &                                  IEXTR,hpvxgas(:,ispin),ngas)

          if (ntest.ge.100) then
            write(luout,*) 'present ',ispin,ica
            call wrt_rstr(luout,rst_res,ngas)
          end if

          ok = ok.and.is_proper_restr_for_hpvx(rst_res(:,:,ica,1,ispin),
     &                                      IHOLE,hpvxgas(:,ispin),ngas)
          ok = ok.and.is_proper_restr_for_hpvx(rst_res(:,:,ica,1,ispin),
     &                                      IPART,hpvxgas(:,ispin),ngas)
          ok = ok.and.is_proper_restr_for_hpvx(rst_res(:,:,ica,1,ispin),
     &                                      IVALE,hpvxgas(:,ispin),ngas)
          ok = ok.and.is_proper_restr_for_hpvx(rst_res(:,:,ica,1,ispin),
     &                                      IEXTR,hpvxgas(:,ispin),ngas)

        end do

        ! mask = 2: paths excluded in occupation graph
        if (.not.zero_ivec(rst1(:,:,:,2,ispin),2*ngas*2).or.
     &      .not.zero_ivec(rst2(:,:,:,2,ispin),2*ngas*2)) then
          call quit(1,i_am,'adapt me for mask restrictions!')
        end if

      end do

      common_restr = ok

      if (ntest.ge.100) then
        write(luout,*) 'on exit: ok = ',ok
        call wrt_rstr(luout,rst_res,ngas)
      end if

      end
