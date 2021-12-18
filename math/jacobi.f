      subroutine jacobi(a,n,np,d,v,nrot)
C-----------------------------------------------------------------------
C     Subroutine used to find the eigenvalues and eigenvectors of a
C     SYMMETRIC matrix, a, with physical dimension, np, and logical
C     dimension, n.
C     Eigenvalues stored in d.
C     Normalised eigenvectors formed and stored in v.
C     nrot gives the total number of steps.
C     Taken from Numerical Recipes Online book: www.nr.com.
C     Slightly modified to deal with double precision reals.
C-----------------------------------------------------------------------
      integer n,np,nrot,nmax
      real*8 a(np,np),d(np),v(np,np)
      parameter (nmax=600)
      integer i,ip,iq,j
      real*8 c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)
      do ip=1,n
         do iq=1,n
            v(ip,iq)=0.d0
         enddo
         v(ip,ip)=1.d0
      enddo
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.d0
      enddo
      nrot=0
      do i=1,50
         sm=0.d0
         do ip=1,n-1
            do iq=ip+1,n
               sm=sm+abs(a(ip,iq))
            enddo
         enddo
         if(sm.eq.0.d0)return
         if(i.lt.4)then
            tresh=0.2d0*sm/(dble(n))**2
         else
            tresh=0.d0
         endif
         do ip=1,n-1
            do iq=ip+1,n
               g=100.d0*abs(a(ip,iq))
               if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     .            .and.(abs(d(iq))+g.eq.abs(d(iq))))then
                  a(ip,iq)=0.d0
               elseif(abs(a(ip,iq)).gt.tresh)then
                  h=d(iq)-d(ip)
                  if(abs(h)+g.eq.abs(h))then
                     t=a(ip,iq)/h
                  else
                     theta=0.5d0*h/a(ip,iq)
                     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                     if(theta.lt.0.d0)t=-t
                  endif
                  c=1.d0/sqrt(1.d0+t**2)
                  s=t*c
                  tau=s/(1.d0+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=0.d0
                  do j=1,ip-1
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  enddo
                  do j=ip+1,iq-1
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  enddo
                  do j=iq+1,n
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
                  enddo
                  do j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
                  enddo
                  nrot=nrot+1
               endif
            enddo
         enddo
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.d0
         enddo
      enddo
      write(6,'(''Too many iterations in jacobi'')')
      return
      end
