*----------------------------------------------------------------------*
*     a set of routines to convert parameters to "actions" into strings
*     and vice versa
*
*     the following routines are contained
*
*       genop_parameters
*       hop_parameters
*       xop_parameters
*       dens_parameters
*       r12gem_parameters
*       cloneop_parameters
*       form_parameters
*       me_list_parameters
*       import_parameters
*       solve_parameters
*
*----------------------------------------------------------------------*
      subroutine genop_parameters(rw,parameters,n_par_str,
     &     dagger,min_rank,max_rank,ncadiff,iformal,explicit,
     &     hpvx_constr,ngastp,gas_constr,ngas)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      logical, intent(inout) ::
     &     explicit, dagger
      integer, intent(inout) ::
     &     min_rank,max_rank,ncadiff,iformal,ngas,ngastp,
     &     hpvx_constr(*),gas_constr(*) 
      character(*), intent(inout) ::
     &     parameters(n_par_str)

      if (n_par_str.ne.3)
     &       call quit(1,'genop_parameters','n_par_str must be 3!')
      if (rw.lt.0) then
        parameters(1)(1:len(parameters(1))) = ' '
        parameters(2)(1:len(parameters(2))) = ' '
        parameters(3)(1:len(parameters(3))) = ' '
        write(parameters(1),'(4(i5,x),2l,2(i5,x))')
     &        min_rank,max_rank,ncadiff,iformal,dagger,explicit,
     &        ngastp,ngas
        if (ngastp.gt.8.or.ngas.gt.12)
     &       call quit(1,'genop_parameters',
     &       'ngastp or ngas larger than expected')
        write(parameters(2),'(64(i2))')
     &        hpvx_constr(1:2*ngastp*2)
        write(parameters(3),'(96(i2))')
     &        gas_constr(1:2*ngas*2*2)
        
      else
        read(parameters(1),'(4(i5,x),2l,2(i5,x))')
     &       min_rank,max_rank,ncadiff,iformal,dagger,explicit,
     &       ngastp,ngas
        read(parameters(2),'(64(i2))')
     &        hpvx_constr(1:2*ngastp*2)
        read(parameters(3),'(96(i2))')
     &        gas_constr(1:2*ngas*2*2)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine op_from_occ_parameters(rw,parameters,n_par_str,
     &     occ_def,ndef,njoined,nmax)

      implicit none
      
      include 'opdim.h'

      integer, intent(in) ::
     &     rw, n_par_str, nmax
      integer, intent(inout) ::
     &     occ_def(*), ndef, njoined
      character(*), intent(inout) ::
     &     parameters(n_par_str)

      if (n_par_str.ne.2)
     &       call quit(1,'genop_parameters','n_par_str must be 2!')
      if (rw.lt.0) then
        parameters(1)(1:len(parameters(1))) = ' '
        parameters(2)(1:len(parameters(2))) = ' '
        write(parameters(1),'(2(i5,x))')
     &        njoined, ndef
c        if (2*ngastp*ndef*njoined.gt.120)
c        if (2*ngastp*ndef*njoined.gt.240)
        if (2*ngastp*ndef*njoined.gt.480)
     &      call quit(1,'op_from_occ_parameters','2*ngastp*ndef.gt.480')
c        write(parameters(2),'(120(i2))')
c        write(parameters(2),'(240(i1))')
c dbg
        print *,'njoined,ndef',njoined,ndef
        print *,'nmax',nmax
        print *,'passed: ',occ_def(1:2*ngastp*ndef*njoined)
c dbg
        write(parameters(2),'(480(i1))')
     &        occ_def(1:2*ngastp*ndef*njoined)        
      else
        read(parameters(1),'(2(i5,x))')
     &       njoined, ndef
c dbg
        print *,'njoined,ndef',njoined,ndef
        print *,'nmax',nmax
c dbg
        if (ndef*njoined.gt.nmax)
     &       call quit(1,'op_from_occ_parameters','nmax too small')
c        read(parameters(2),'(120(i2))')
c        read(parameters(2),'(240(i1))')
        read(parameters(2),'(480(i1))')
     &        occ_def(1:2*ngastp*ndef*njoined)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine hop_parameters(rw,parameters,
     &     min_rank,max_rank,iformal,explicit)

      implicit none
      
      logical, intent(inout) ::
     &     explicit
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(8(i5,x),l)')
     &        min_rank,max_rank,iformal,explicit
      else
        read(parameters,'(8(i5,x),l)')
     &       min_rank,max_rank,iformal,explicit
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine xop_parameters(rw,parameters,
     &     dagger,min_rank,max_rank,ncadiff,iformal)

      implicit none
      
      logical, intent(inout) ::
     &     dagger
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,ncadiff,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(l,x,4(i5,x))')
     &       dagger,min_rank,max_rank,ncadiff,iformal
      else
        read(parameters,'(l,x,4(i5,x))')
     &       dagger,min_rank,max_rank,ncadiff,iformal
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine dens_parameters(rw,parameters,
     &     min_rank,max_rank,iformal)

      implicit none
      
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(8(i5,x))')
     &        min_rank,max_rank,iformal
      else
        read(parameters,'(8(i5,x))')
     &        min_rank,max_rank,iformal
      end if

      return
      end


*----------------------------------------------------------------------*
      subroutine r12gem_parameters(rw,parameters,
     &     n_ap,min_rank,max_rank,ansatz)

      implicit none
      
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,ansatz,n_ap
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(i2,x,3(i5,x))')
     &       n_ap,min_rank,max_rank,ansatz
      else
        read(parameters,'(i2,x,3(i5,x))')
     &       n_ap,min_rank,max_rank,ansatz
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine r12int_parameters(rw,parameters,
     &     n_ap,min_rank,max_rank,ncadiff,iformal)

      implicit none
      
      integer, intent(inout) ::
     &     n_ap
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,ncadiff,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(i2,x,4(i5,x))')
     &       n_ap,min_rank,max_rank,ncadiff,iformal
      else
        read(parameters,'(i2,x,4(i5,x))')
     &       n_ap,min_rank,max_rank,ncadiff,iformal
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine cloneop_parameters(rw,parameters,name,dagger)

      implicit none
      
      integer, intent(in) ::
     &     rw
      logical, intent(inout) ::
     &     dagger
      character, intent(inout) ::
     &     parameters*(*), name*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(l,x,a)')
     &       dagger,name
      else
        read (parameters,'(l,x,a)')
     &       dagger,name
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine form_parameters(rw,
     &     parameters,n_par_str,title,inum,mode_str)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      integer, intent(inout) ::
     &     inum
      character*(*), intent(inout) ::
     &     parameters(n_par_str),
     &     title, mode_str

      if (rw.lt.0) then
        write(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       write(parameters(2),'(i8,a)') inum, mode_str
      else
        read(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       read(parameters(2),'(i8,a)') inum, mode_str
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine expand_parameters(rw,
     &     parameters,n_par_str,title,nop,
     &     idx_sv,iblkmin,iblkmax,
     &     connect,nconnect,
     &     avoid,navoid,
     &     inproj,ninproj)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      integer, intent(inout) ::
     &     nop,
     &     idx_sv(*),iblkmin(*),iblkmax(*),
     &     connect(*),nconnect,
     &     avoid(*),navoid,
     &     inproj(*),ninproj
      character*(*), intent(inout) ::
     &     parameters(n_par_str),
     &     title

      integer ::
     &     ii

      if (n_par_str.lt.3)
     &     call quit(1,'expand_parameters','3 strings!')

      if (rw.lt.0) then
        write(parameters(1),'(a)') title
        write(parameters(2),'(50i4,a)') nop,
     &       (idx_sv(ii),iblkmin(ii),iblkmax(ii), ii = 1, nop)
        write(parameters(3),'(50i4,a)') nconnect,navoid,ninproj,
     &       connect(1:2*nconnect), avoid(1:2*navoid),
     &       inproj(1:4*ninproj)
      else
        read(parameters(1),'(a)') title
        read(parameters(2),'(50i4,a)') nop,
     &       (idx_sv(ii),iblkmin(ii),iblkmax(ii), ii = 1, nop)
        read(parameters(3),'(50i4,a)') nconnect,navoid,ninproj,
     &       connect(1:2*nconnect), avoid(1:2*navoid),
     &       inproj(1:4*ninproj)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine form_parameters2(rw,
     &     parameters,n_par_str,title,inum,idxlist)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      integer, intent(inout) ::
     &     inum, idxlist(*)
      character*(*), intent(inout) ::
     &     parameters(n_par_str),
     &     title

      if (rw.lt.0) then
        write(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       write(parameters(2),'(i4,20i4)') inum, idxlist(1:inum)
      else
        read(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       read(parameters(2),'(i4,20i4)') inum, idxlist(1:inum)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine opt_parameters(rw,parameters,ncat,nint)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     ncat,nint
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(2i4,a)') ncat,nint
      else
        read(parameters,'(2i4,a)') ncat,nint
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine me_list_parameters(rw,parameters,
     &     absym,casym,gamma,s2,ms,ms_fix)

      implicit none
      
      integer, intent(inout) ::
     &     rw,absym,casym,gamma,s2,ms
      logical, intent(inout) ::
     &     ms_fix
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(l,8(i5,x))')
     &        ms_fix,absym,casym,gamma,s2,ms
      else
        read(parameters,'(l,8(i5,x))')
     &       ms_fix,absym,casym,gamma,s2,ms
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine import_parameters(rw,parameters,
     &     env_type)

      implicit none
      
      integer, intent(in) ::
     &     rw
      character(*), intent(inout) ::
     &     env_type
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(a)') env_type
      else
        read(parameters,'(a)') env_type
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine solve_parameters(rw,parameters,n_par_str,
     &     nopt,nroots,mode_str)

      implicit none
      
      integer, intent(in) ::
     &     rw,n_par_str
      integer, intent(inout) ::
     &     nopt,nroots
      character(*), intent(inout) ::
     &     parameters(n_par_str), mode_str

      if (n_par_str.lt.2)
     &     call quit(1,'solve_parameters','need >= 2 par-strings!')
      if (rw.lt.0) then
        parameters(1)(1:len(parameters)) = ' '
        parameters(2)(1:len(parameters)) = ' '
        write(parameters(1),'(2i4)') nopt,nroots
        write(parameters(2),'(a)')   trim(mode_str)
      else
        read(parameters(1),'(2i4)') nopt,nroots
        read(parameters(2),'(a)')   mode_str
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine evalprop_parameters(rw,parameters,
     &     ndens,rank,env_type)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     ndens,rank
      character(*), intent(inout) ::
     &     parameters, env_type

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(2i4,a)') ndens,rank,env_type
      else
        read(parameters,'(2i4,a)') ndens,rank,env_type
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine add_parameters(rw,parameters,
     &     nfac,fac,maxfac)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     nfac,maxfac
      real(8), intent(inout) ::
     &     fac(*)
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        if (nfac.gt.12) call quit(1,'add_parameters','too much')
        write(parameters,'(i2,12g20.14)') nfac,fac(1:nfac)
      else
        read(parameters,'(i2,12g20.14)') nfac,fac(1:nfac)
        if (nfac.gt.maxfac)
     &       call quit(1,'add_parameters','too much (>maxfac)')
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine scale_parameters(rw,parameters,
     &     nfac,idxblk,fac,maxfac)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     nfac,maxfac,idxblk(*)
      real(8), intent(inout) ::
     &     fac(*)
      character(*), intent(inout) ::
     &     parameters
      integer ::
     &     ii

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        if (nfac.gt.12) call quit(1,'scale_parameters','too much')
        write(parameters,'(i2,12(i4,g20.14))')
     &       nfac,((idxblk(ii),fac(ii)), ii=1,nfac)
      else
        read(parameters,'(i2,12(i4,g20.14))')
     &       nfac,((idxblk(ii),fac(ii)), ii=1,nfac)
        if (nfac.gt.maxfac)
     &       call quit(1,'scale_parameters','too much (>maxfac)')
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine modify_parameters(rw,parameters,
     &     nterm,idxterm,maxterm)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     nterm,maxterm,idxterm(*)
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        if (nterm.gt.30) call quit(1,'modify_parameters','too much')
        write(parameters,'(i2,30i4)')
     &       nterm,idxterm(1:nterm)
      else
        read(parameters,'(i2,30i4)')
     &       nterm,idxterm(1:nterm)
        if (nterm.gt.maxterm)
     &       call quit(1,'modify_parameters','too much (>maxterm)')
      end if

      return
      end

     
*---------------------------------------------------------------------*

      subroutine ord_parameters(rw,parameters,
     &     iorder)

      implicit none

      integer, intent(inout) ::
     &     rw,iorder
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(i5)')
     &        iorder
      else
        read(parameters,'(i5)')
     &       iorder
      end if

      return
      end

