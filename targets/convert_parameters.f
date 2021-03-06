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
     &     min_rank,max_rank,ncadiff,
     &     min_xrank,max_xrank,iformal,freeze,
     &     hpvx_constr,hpvxca_constr,ngastp,gas_constr,ngas)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      logical, intent(inout) ::
     &     freeze(2)
      integer, intent(inout) ::
     &     min_rank,max_rank,min_xrank,max_xrank,
     &     ncadiff,iformal,ngas,ngastp,
     &     hpvx_constr(*),hpvxca_constr(*),gas_constr(*) 
      character(*), intent(inout) ::
     &     parameters(n_par_str)

      if (n_par_str.ne.3)
     &       call quit(1,'genop_parameters','n_par_str must be 3!')
      if (rw.lt.0) then
        parameters(1)(1:len(parameters(1))) = ' '
        parameters(2)(1:len(parameters(2))) = ' '
        parameters(3)(1:len(parameters(3))) = ' '
        write(parameters(1),'(8(i5,x),2(l2))')
     &        min_rank,max_rank,min_xrank,max_xrank,ncadiff,iformal,
     &        ngastp,ngas,freeze(1:2)
        if (ngastp.gt.8.or.ngas.gt.12)
     &       call quit(1,'genop_parameters',
     &       'ngastp or ngas larger than expected')
        write(parameters(2),'(96(i2))')
     &        hpvx_constr(1:2*ngastp),hpvxca_constr(1:2*ngastp*2)
        write(parameters(3),'(96(i2))')
     &        gas_constr(1:2*ngas*2*2)
        
      else
        read(parameters(1),'(8(i5,x),2(l2))')
     &        min_rank,max_rank,min_xrank,max_xrank,ncadiff,iformal,
     &        ngastp,ngas,freeze(1:2)
        read(parameters(2),'(96(i2))')
     &        hpvx_constr(1:2*ngastp),hpvxca_constr(1:2*ngastp*2)
        read(parameters(3),'(96(i2))')
     &        gas_constr(1:2*ngas*2*2)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine op_from_occ_parameters(rw,parameters,n_par_str,
     &     occ_def,ndef,njoined,freeze,nmax)

      implicit none
      
      include 'stdunit.h'
      include 'opdim.h'

      integer, intent(in) ::
     &     rw, n_par_str, nmax
      integer, intent(inout) ::
     &     occ_def(*), ndef, njoined
      integer, intent(inout) ::
     &     freeze(2,*)
      character(*), intent(inout) ::
     &     parameters(n_par_str)

      if (n_par_str.ne.2)
     &       call quit(1,'op_from_occ_parameters'
     &                  ,'n_par_str must be 2!')
      if (rw.lt.0) then
        parameters(1)(1:len(parameters(1))) = ' '
        parameters(2)(1:len(parameters(2))) = ' '
        write(parameters(1),'(2(i5,x),10(i2))')
     &        njoined, ndef, freeze(1:2,1:njoined)
c        if (2*ngastp*ndef*njoined.gt.120)
c        if (2*ngastp*ndef*njoined.gt.240)
        if (2*ngastp*ndef*njoined.gt.2047) then
           write(lulog,*) 'ndef, njoined: ',ndef,njoined
           call quit(1,'op_from_occ_parameters','2*ngastp*ndef.gt.2047')
        end if
c        write(parameters(2),'(120(i2))')
c        write(parameters(2),'(240(i1))')
c        write(parameters(2),'(480(i1))')
        write(parameters(2),'(2047(i2))')
     &        occ_def(1:2*ngastp*ndef*njoined)        
      else
        read(parameters(1),'(2(i5,x),10(i2))')
     &       njoined, ndef, freeze(1:2,1:njoined)
        if (ndef*njoined.gt.nmax)
     &       call quit(1,'op_from_occ_parameters','nmax too small')
c        read(parameters(2),'(120(i2))')
c        read(parameters(2),'(240(i1))')
c        read(parameters(2),'(480(i1))')
c dbg
c      print *,'"',trim(parameters(2)),'"'
c dbg
        read(parameters(2),'(2047(i2))')
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
        write(parameters,'(3(i5,x),l)')
     &        min_rank,max_rank,iformal,explicit
      else
        read(parameters,'(3(i5,x),l)')
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
      subroutine def_form_parameters(rw,
     &     parameters,n_par_str,form_str,title)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      character*(*), intent(inout) ::
     &     parameters(n_par_str),
     &     title, form_str

      if (rw.lt.0) then
        write(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       write(parameters(2),'(a)') form_str
      else
        read(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       read(parameters(2),'(a)') form_str
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
     &     parameters,n_par_str,title,inum,idxlist,inum2,pop_idx)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      integer, intent(inout) ::
     &     inum, idxlist(*)
      character*(*), intent(inout) ::
     &     parameters(n_par_str),
     &     title
      integer, intent(inout), optional ::
     &     inum2, pop_idx(*)

      if (rw.lt.0) then
        write(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       write(parameters(2),'(i4,20i4)') inum, idxlist(1:inum)
        if (n_par_str.gt.2)
     &       write(parameters(3),'(i4,100i2)') inum2, pop_idx(1:inum2)
      else
        read(parameters(1),'(a)') title
        if (n_par_str.gt.1)
     &       read(parameters(2),'(i4,20i4)') inum, idxlist(1:inum)
        if (n_par_str.gt.2)
     &       read(parameters(3),'(i4,100i2)') inum2, pop_idx(1:inum2)
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
     &     list_type,env_type)

      implicit none
      
      integer, intent(in) ::
     &     rw
      character(*), intent(inout) ::
     &     list_type, env_type
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(a)') trim(list_type)
        write(parameters(33:),'(a)') trim(env_type)
      else
        read(parameters,'(2a32)') list_type,env_type
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
     &     imode,nfac,idxblk,fac,maxfac)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     imode,nfac,maxfac,idxblk(*)
      real(8), intent(inout) ::
     &     fac(*)
      character(*), intent(inout) ::
     &     parameters
      integer ::
     &     ii

      if (nfac.gt.1) call quit(1,'scale_parameter','obsolete. needed?')

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        if (nfac.gt.12) call quit(1,'scale_parameters','too much')
        write(parameters,'(2i2,12(i4,g20.14))')
     &       imode,nfac,idxblk(1),fac(1)
C     &       imode,nfac,((idxblk(ii),fac(ii)), ii=1,nfac)
      else
        read(parameters,'(2i2,12(i4,g20.14))')
     &       imode,nfac,idxblk(1),fac(1)
C     &       imode,nfac,((idxblk(ii),fac(ii)), ii=1,nfac)
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
*----------------------------------------------------------------------*
      subroutine select_parameters(rw,parameters,n_par_str,
     &     ninclude,ninclude_or,nexclude,
     &     iblk_include,iblk_include_or,iblk_exclude)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     n_par_str,
     &     ninclude,ninclude_or,nexclude,
     &     iblk_include(*),iblk_include_or(*),iblk_exclude(*)
      character(len=*), intent(inout) ::
     &     parameters(n_par_str)

      if (rw.lt.0) then
        parameters(1)(1:len(parameters(1))) = ' '
        parameters(2)(1:len(parameters(2))) = ' '
        parameters(3)(1:len(parameters(3))) = ' '
        if (ninclude.gt.30) call quit(1,'select_parameters','too much')
        if (ninclude_or.gt.30)
     &                      call quit(1,'select_parameters','too much')
        if (nexclude.gt.30) call quit(1,'select_parameters','too much')
        if (ninclude.gt.0) then
          write(parameters(1),'(i2,30i4)')
     &       ninclude,iblk_include(1:ninclude)
        else
          write(parameters(1),'(" 0   0")')
        end if
        if (ninclude_or.gt.0) then
          write(parameters(2),'(i2,30i4)')
     &       ninclude_or,iblk_include_or(1:ninclude_or)
        else
          write(parameters(2),'(" 0   0")')
        end if
        if (nexclude.gt.0) then
          write(parameters(3),'(i2,30i4)')
     &       nexclude,iblk_exclude(1:nexclude)
        else
          write(parameters(3),'(" 0   0")')
        end if
      else
        read(parameters(1),'(i2,30i4)')
     &       ninclude,iblk_include(1:ninclude)
        read(parameters(2),'(i2,30i4)')
     &       ninclude_or,iblk_include_or(1:ninclude_or)
        read(parameters(3),'(i2,30i4)')
     &       nexclude,iblk_exclude(1:nexclude)
      end if

      return
      end

     
*---------------------------------------------------------------------*

      subroutine ord_parameters(rw,parameters,
     &     iorder,species,ifreq)

      implicit none

      integer, intent(inout) ::
     &     rw, iorder, species, ifreq(*)
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,*)
     &        iorder, species, ifreq(1:iorder)
      else
        read(parameters,*)
     &       iorder, species, ifreq(1:iorder)
      end if

      return
      end

*---------------------------------------------------------------------*

      subroutine freq_parameters(rw,parameters,freq)

      implicit none

      integer, intent(inout) ::
     &     rw
      real(8), intent(inout) ::
     &     freq
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(g20.14)') freq
      else
        read(parameters,'(g20.14)') freq
      end if

      return
      end
