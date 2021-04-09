c ... program ...
        program fmg

        implicit none

C ... variable declarations ...

        integer i,k,capm,capn,nx1,ny1,ncycles,ic(10),gamma

        include 'incs.com'
        
        real h1,u(0:gptmax),f(0:gptmax),ex(0:gptmax)

C ... input and output files ...

        open(unit=1,file='main.in',status='old')
        open(unit=6,file='fmg.out',status='new')

C ... read input data ...

        read(1,*) nx1,ny1,h1
        read(1,*) capm
        read(1,*) gamma
        read(1,*) ncycles

C ... zero out arrays ...

        do 10 i=0,gptmax
          u(i) = 0.
          f(i) = 0.
          ex(i)= 0.
 10     continue

C ... calculate multigrid specifiers ...

        call grdfn(h1,nx1,ny1,capm)

C ... generate rhs and exact solution ...

        do 20 k=1,capm
          call fset(k,f,ex)
 20     continue

C .. full multigrid algorithm ...

        do 35 capn = 2,capm

C       ... solve on coarse grid exactly ... 

          if(capn.eq.2) call relax(1,u,f)

C       ... fmg interpolation ...

          call cubicint(capn-1,capn,u,f)
          
C       ... (2,1) V or W cycles ...

          ic(capn) = ncycles
          k = capn
 25       if(ic(capn).gt.0) then
            if((ic(k).eq.0) .or. (k.eq.1)) then
              if(k.eq.1) then
                call relax(1,u,f)
                ic(1)=ic(1)-1
              endif
              call intadd(k,k+1,u)
              k = k+1
              call relax(k,u,f)
              if(k.eq.capn) call err(k,u,ex,f)
              ic(k) = ic(k) - 1
            else
              call relax(k,u,f)
              call relax(k,u,f)
              call fwrestrn(k,k-1,u,f)
              call putzero(k-1,u)
              k = k-1
              ic(k) = gamma
            endif
            goto 25
          endif

 35     continue

C ... stop ...

        stop
        end
c ... subroutines ...
        subroutine fwrestrn(k,kc,v,f)

        implicit none

C ... variable declarations ...

        integer i,j,ic,jc,imax,jmax,imaxc,jmaxc,k,kc,jeo,istart

        include 'incs.com'

        integer p(0:npmax),pc(0:npmax)
        real h,hc,v(0:gptmax),f(0:gptmax),r(gptmax)

C ... calculate grid specifiers ...

        call key(k,p,imax,jmax,h)
        call key(kc,pc,imaxc,jmaxc,hc)

C ... compute residual at even points ...

        jeo = 1

        do 20 j=1,jmax
          jeo = -jeo
          istart = 1 + (jeo+1)/2
          do 20 i=istart,imax,2
            r(p(i)+j) = f(p(i)+j) + v(p(i)+j+1)
     &           + v(p(i)+j-1) + v(p(i-1)+j)
     &           + v(p(i+1)+j) - 4.*v(p(i)+j)
 20     continue

C ... restrict using full-weighting ...

        do 30 jc=1,jmaxc
          j=2*jc
          do 30 ic=1,imaxc
            i=2*ic
            f(pc(ic)+jc) = 0.25*(4.*r(p(i)+j)
     &                       + r(p(i-1)+j+1) + r(p(i-1)+j-1)
     &                       + r(p(i+1)+j+1) + r(p(i+1)+j-1))
 30     continue

C ... return ...

        return
        end 
        subroutine fset(k,g,ex)

        implicit none

C ... variable declarations ...

        integer i,j,k,imax,jmax

        real x,xx,x6,y,yy,h,h2

        include 'incs.com'

        integer p(0:npmax)
        real g(0:gptmax),ex(0:gptmax)

C ... specifiers for grid k ...

        call key(k,p,imax,jmax,h)

C ... generate rhs and exact solution ...

        h2 = h*h

        do 10 i=1,imax
          x = h*real(i)
          xx = x*x*(x*x - 1.)
          x6 = 6.*x*x - 1.

          do 10 j=1,jmax
            y  = h*real(j)
            yy = y*y*(y*y - 1.)
            ex(p(i)+j) = -xx*yy
            g(p(i)+j) = 2.*h2*(x6*yy + xx*(6.*y*y - 1.))
 10     continue

C ... return ...

        return
        end
        subroutine cubicint(kc,k,v,g)

        implicit none

C ... variable declarations ...

        integer i,j,ic,jc,imax,jmax,imaxc,jmaxc

        include 'incs.com'

        integer pc(0:npmax),p(0:npmax),kc,k
        real hc,h,g(0:gptmax),v(0:gptmax)

C ... call grid specifiers ...

        call key(kc,pc,imaxc,jmaxc,hc)
        call key(k,p,imax,jmax,h)

C ... interpolate solution to finer grid ...

C    ... even points ...

        do 10 jc=1,jmaxc
          j=2*jc
          do 10 ic = 1,imaxc 
            i=2*ic
            v(p(i)+j) = v(pc(ic)+jc)
 10     continue

C       ... rotated 5 point star ...

        do 20 jc=1,jmaxc+1
          j=2*jc
          do 20 ic=1,imaxc+1
            i=2*ic
            v(p(i-1)+j-1) = 0.5*g(p(i-1)+j-1) 
     &                      + 0.25*(v(pc(ic)+jc) + v(pc(ic-1)+jc)
     &                      + v(pc(ic)+jc-1) + v(pc(ic-1)+jc-1))
 20     continue

C    ... odd points ...

        do 30 jc=1,jmaxc
          j=2*jc
          do 30 ic=1,imaxc
            i=2*ic
            v(p(i-1)+j) = 0.25*(g(p(i-1)+j) +v(p(i)+j) +v(p(i-1)+j-1)
     &                          + v(p(i-1)+j+1) + v(p(i-2)+j))
            v(p(i)+j-1) = 0.25*(g(p(i)+j-1) + v(p(i)+j) + v(p(i)+j-2)
     &                          + v(p(i-1)+j-1) + v(p(i+1)+j-1))
 30     continue

C       ... along J+1 boundary ...

        j=2*(jmaxc+1)
        do 40 ic=1,imaxc
          i=2*ic
          v(p(i)+j-1) = 0.25*(g(p(i)+j-1) + v(p(i)+j) + v(p(i)+j-2)
     &                        + v(p(i-1)+j-1) + v(p(i+1)+j-1))
 40     continue

C       ... along I+1 boundary ...

        i=2*(imaxc+1)
        do 50 jc=1,jmaxc
          j=2*jc
          v(p(i-1)+j) = 0.25*(g(p(i-1)+j) +v(p(i)+j) +v(p(i-1)+j-1)
     &                        + v(p(i-1)+j+1) + v(p(i-2)+j))     
 50     continue

C ... zero out present grid approximation ...

        do 60 i=0,imaxc+1
          do 60 j=0,jmaxc+1
            v(pc(i)+j) = 0.
 60     continue

C ... return ...

        return
        end
        subroutine err(k,u,ex,f)

        implicit none

C ... variable declarations ...

        integer i,j,k,imax,jmax

        real emax,esum,rmax,rsum,r,er,eh,rh,h

        include 'incs.com'

        integer p(0:npmax)
        real u(0:gptmax),ex(0:gptmax),f(0:gptmax)

C ... specifiers for grid k ...

        call key(k,p,imax,jmax,h)

C ... compute max norm of global and residual error ...

        emax = 0.
        esum = 0.
        rmax = 0.
        rsum = 0.

        do 10 i= 1,imax
          do 10 j=1,jmax

            er = abs(u(p(i)+j) - ex(p(i) + j))
            esum = esum + er*er
            if(er .gt. emax) emax = er

            r = abs(f(p(i)+j) + u(p(i+1)+j) + u(p(i-1)+j) + 
     &              u(p(i)+j+1) + u(p(i)+j-1) - 4.*u(p(i)+j))
            rsum = rsum + r*r
            r = r/(h*h)
            if(r .gt. rmax) rmax = r

 10     continue

C ... compute discrete L2 norm of global and residual error ...

        eh = sqrt(h*h*esum)
        rh = sqrt(rsum)/h

C ... write to text file ...

        write(6,20) k,1./h,rmax,rh,emax,eh
 20     format('  k = ',i3,5e12.3)

C ... return ...

        return
        end


        subroutine grdfn(h1,nx1,ny1,capm)

        implicit none

C ... variable specifications ...

        integer k,k2,iq,imax,jmax,nx1,ny1,capm
        integer imx(20),jmx(20),nst(20)

        real h1,h(20)

C ... common statements ...

        common/grd/nst,imx,jmx,h

C ... calculate multigrid specifiers ...

        iq = 0

        do 20 k=1,capm
            k2 = 2**(k-1)
          nst(k) = iq
          imx(k) = nx1*k2 - 1
          jmx(k) = ny1*k2 - 1
          h(k)   = h1/k2
            iq = iq + (imx(k)+2)*(jmx(k)+2)
 20     continue

C ... return ...

        return
        end

        subroutine intadd(kc,k,v)

        implicit none

C ... variable declarations ...

        integer ic,jc,i,j,imax,imaxc,jmax,jmaxc,k,kc

        include 'incs.com'

        integer p(0:npmax),pc(0:npmax)
        real h,hc,v(0:gptmax)

C ... specifiers for grid k ...

        call key(kc,pc,imaxc,jmaxc,hc)
        call key(k,p,imax,jmax,h)

C ... interpolate correction and add to fine grid approximation ...

C   ... interior block ...

        do 10 ic = 1,imaxc
          i = 2*ic

          do 10 jc = 1,jmaxc
            j = 2*jc
            v(p(i)+j)     = v(p(i)+j) + v(pc(ic)+jc)
            v(p(i-1)+j)   = v(p(i-1)+j) 
     &                      + 0.5*(v(pc(ic)+jc) + v(pc(ic-1)+jc))
            v(p(i)+j-1)   = v(p(i)+j-1) 
     &                      + 0.5*(v(pc(ic)+jc) + v(pc(ic)+jc-1))
            v(p(i-1)+j-1) = v(p(i-1)+j-1) + 0.25*(v(pc(ic)+jc) 
     &                      + v(pc(ic-1)+jc) + v(pc(ic)+jc-1) 
     &                      + v(pc(ic-1)+jc-1))
 10     continue

C   ... points along J+1 boundary ...

        jc = jmaxc+1
        j = 2*jc

        do 20 ic=1,imaxc
          i=2*ic
          v(p(i)+j-1)   = v(p(i)+j-1) 
     &                      + 0.5*(v(pc(ic)+jc) + v(pc(ic)+jc-1))
          v(p(i-1)+j-1) = v(p(i-1)+j-1) + 0.25*(v(pc(ic)+jc) 
     &                      + v(pc(ic-1)+jc) + v(pc(ic)+jc-1) 
     &                      + v(pc(ic-1)+jc-1))
 20     continue

C   ... points along I+1 boundary ...

        ic = imaxc+1
        i = 2*ic

        do 30 jc=1,jmaxc
          j=2*jc
          v(p(i-1)+j)   = v(p(i-1)+j) 
     &                    + 0.5*(v(pc(ic)+jc) + v(pc(ic-1)+jc))
          v(p(i-1)+j-1) = v(p(i-1)+j-1) + 0.25*(v(pc(ic)+jc) 
     &                    + v(pc(ic-1)+jc) + v(pc(ic)+jc-1) 
     &                    + v(pc(ic-1)+jc-1))
 30     continue

C   ... interior corner point ...

        ic = imaxc+1
        i  = 2*ic
        jc = jmaxc+1
        j  = 2*jc

        v(p(i-1)+j-1) = v(p(i-1)+j-1) + 0.25*(v(pc(ic)+jc) 
     &                  + v(pc(ic-1)+jc) + v(pc(ic)+jc-1) 
     &                  + v(pc(ic-1)+jc-1))

C ... return ...

        return
        end
        subroutine key(k,ist,imax,jmax,hh)

        implicit none

C ... variable declarations ...

        integer i,imax,jmax,is,k
        integer nst(20),imx(20),jmx(20)

        include 'incs.com'

        integer ist(0:npmax)

        real hh,h(20)

C ... common statements ...

        common/grd/nst,imx,jmx,h

C ... determine specifiers for grid k ...

        imax = imx(k)
        jmax = jmx(k)
        hh   = h(k)

        is = nst(k)

        do 10 i=0,imax+1
          ist(i) = is      
          is = is + (jmax+2)
 10     continue

C ... return ...

        return
        end
        subroutine putzero(k,u)

        implicit none

C ... implicit none ...

        integer i,j,k,imax,jmax

        include 'incs.com'

        integer p(0:npmax)
        real u(0:gptmax),h

C ... specifiers for grid k ...

        call key(k,p,imax,jmax,h)

C ... initialize to zero ...

        do 10 i=0,imax+1
          do 10 j=0,jmax+1
            u(p(i)+j) = 0.
 10     continue

C ... return ...

        return
        end
        subroutine relax(k,v,g)

        implicit none

C ... variable declarations ...

        integer i,j,k,krb,irb,j0,jdel,imax,jmax

        include 'incs.com'

        integer p(0:npmax)
        real v(0:gptmax),g(0:gptmax),h

C ... specifiers for grid k ...

        call key(k,p,imax,jmax,h)

C ... red-black gauss-seidel on (n+1)**2 grid ...

        do 10 krb=1,-1,-2

          irb=krb

          do 20 i=1,imax
            irb= -irb
            jdel = (irb+1)/2
            j0 = jdel + 1

            do 20 j=j0,jmax,2
              v(p(i)+j) = (g(p(i)+j) + v(p(i+1)+j) + v(p(i-1)+j)
     &                     + v(p(i)+j+1) + v(p(i)+j-1))/4.
 20       continue

 10     continue

C ... return ...

        return
        end
