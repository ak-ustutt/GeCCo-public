      logical function is_dalton64() result(is64)

      implicit none
      
      include 'ioparam.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'ifc_baserout.h'
      include 'par_dalton.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 100

      type(filinf) ::
     &     ffsir
      integer ::
     &     lusir, luerr, iprint, ngas, loop, i, j, caborb, idx, jdx,
     &     n_frozen, n_act, n_as1, n_as2, n_as3, nspin
      logical ::
     &     logaux

      integer, parameter ::
     &     mxsym = 8
      integer(4) ::
     &     istate,ispin,nactel,lsym,nsym,
     &     nisht,nasht,nocct,norbt,nbast,nxbast,nconf,nwopt,nwoph,
     &     ncdets, ncmot,nnashx,nnashy,nnorbt,n2orbt,
     &     nelmn1, nelmx1, nelmn3, nelmx3, mctype,
     &     nish(mxsym),nash(mxsym),norb(mxsym),nbas(mxsym),
     &     nfro(mxsym),nrhf(mxsym),
     &     muld2h(mxsym,mxsym), nas1(mxsym), nas2(mxsym), nas3(mxsym),
     &     auxbas(mxsym),linind(mxsym),totbas(mxsym)
      integer ::
     &     inactorb(mxsym),cas(mxsym),actel
      real(8) ::
     &     potnuc,emy,eactiv,emcscf

      real(8), pointer ::
     &     orb_en(:), orb_en_sort(:)
      integer(4), pointer ::
     &     orb_sym_i4(:)

      call file_init(ffsir,sirifc,ftyp_sq_unf,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      luerr = lulog
      rewind lusir
      call mollab('SIR IPH ',lusir,luerr)

      read (lusir) potnuc,emy,eactiv,emcscf,istate,ispin,nactel,lsym
      read (lusir) nisht,nasht,nocct,norbt,nbast,nconf,nwopt,nwoph,
     &     ncdets, ncmot,nnashx,nnashy,nnorbt,n2orbt,
     &     nsym,muld2h(1:mxsym,1:mxsym),nrhf(1:mxsym),nfro(1:mxsym),
     &     nish(1:mxsym),nash(1:mxsym),norb(1:mxsym),nbas(1:mxsym),
     &     nelmn1, nelmx1, nelmn3, nelmx3, mctype,
     &     nas1(1:mxsym), nas2(1:mxsym), nas3(1:mxsym)

      ! check for inconsistent input -> set to 64 bit then
      if (ispin.eq.0.or.lsym.eq.0.or.norb(1).eq.0) then
        write(lulog,*) 'Trying 64-bit read-in routines ...'
        is64 = .true.
      else
        write(lulog,*) 'Trying 32-bit read-in routines ...'
        is64 = .false.
      end if
      call file_close_keep(ffsir)

      return
      end
