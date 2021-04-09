c ... program ...
        program bgetest

        implicit none

C ... variable declarations ...

        integer k,i
        integer nx1,ny1,capm

        include 'incs.com'

        logical phaseone

        real h1,ewcf,nscf,g(0:gptmax),v(0:gptmax),ex(0:gptmax)

C ... common statements ...

        common/psicoefs/ewcf,nscf

C ... input and output files ...

        open(unit=1,file='main.in',status='old')
        open(unit=6,file='bge.out',status='new')

C ... read input data ...

        read(1,*) nx1,ny1,h1
        read(1,*) capm

C ... parameter assignments ...

        ewcf = 1.
        nscf = 1.

C ... zero out arrays ...

        do 5 i=0,gptmax
          v(i) = 0.
          g(i) = 0.
          ex(i) = 0.
 5      continue

C ... grid specifiers ...

        call grdfn(h1,nx1,ny1,capm)

C ... generate rhs ...

        do 10 k=1,capm
          call fset(k,g,ex)
 10     continue

C ... block gaussian elimination ...

        do 20 k=2,capm
          phaseone = .true.
          call bge(k,v,g,phaseone)
          call err(k,v,ex,g)
 20     continue

C ... stop ...

        stop 
        end

c ... subroutines ...
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


        subroutine bge(k,v,g,phaseone)

        implicit none

C ... variable declarations ...

        integer k,icap,jcap

        logical phaseone

        include 'incs.com'

        integer p(0:npmax)

        real ewcf,nscf,v(0:gptmax),g(0:gptmax),h

C ... common statements ...

        common/psicoefs/ewcf,nscf

C ... specifiers for grid k ...

        call key(k,p,icap,jcap,h)

C ... phase one - block LU decomposition ...

        if(phaseone) then
          call blocklud(icap,jcap,nscf,ewcf)
          phaseone = .false.
        endif

C ... phase two - forward & backward substitution ...

        call blockfbs(ewcf,icap,jcap,p,v,g)

C ... return ...

        return 
        end
        subroutine blocklud(icap,jcap,nscf,ewcf)

        implicit none

C ... variable declarations ...

        integer i,icap,jcap,m,n

        include 'incs.com'

        real nscf,ewcf
        real tcap(bnmax,bnmax,bnmax)
        real acap(bnmax,bnmax),ahat(bnmax,bnmax),qq(bnmax,bnmax)

C ... common statements ...

        common/invblock/tcap

C ... set up A matrix ...

        call diagmatx(nscf,ewcf,jcap,acap)

C ... Block LU Decomposition ...

C     ... i = 1 ...

        do 5 m = 1,jcap
          do 5 n = 1,jcap
            ahat(m,n) = acap(m,n)
 5      continue

        call ludecomp(jcap,ahat)
        call fbsubstn(jcap,ahat,qq)

        do 10 m = 1,jcap
          do 10 n = 1,jcap
            tcap(m,n,1) = qq(m,n)
 10     continue

C     ... i = 2,icap ...

        do 20 i=2,icap

          do 25 m = 1,jcap
            do 25 n = 1,jcap
              ahat(m,n) = acap(m,n) - ewcf*ewcf*qq(m,n)
 25       continue

          call ludecomp(jcap,ahat)
          call fbsubstn(jcap,ahat,qq)

          do 30 m = 1,jcap
            do 30 n = 1,jcap
              tcap(m,n,i) = qq(m,n)
 30       continue

 20     continue

C ... return ...

        return 
        end

        subroutine diagmatx(nscf,ewcf,n,acap)

        implicit none

C ... variable declarations ...

        integer i,j,n

        include 'incs.com'

        real diag,nscf,ewcf
        real acap(bnmax,bnmax)

C ... zero out A matrix ...

        do 5 i=i,n
          do 5 j=1,n
            acap(i,j) = 0.
 5      continue

C ... generate A matrix ...

        diag = 2.*(nscf + ewcf)

C     ... i=1 ...

        acap(1,1) = diag
        acap(1,2) = -nscf

C     ... i=2,n-1 ...

        do 10 i=2,n-1
          acap(i,i-1) = -nscf
          acap(i,i)   = diag
          acap(i,i+1) = -nscf
 10     continue

C     ... i=n ...

        acap(n,n)   = diag
        acap(n,n-1) = -nscf

C ... return ...

        return
        end
        subroutine fbsubstn(n,aa,ww)

        implicit none

C ... variable declarations ...

        integer i,j,k,n

        include 'incs.com'

        real sum,krd,aa,ww,yy

        dimension aa(bnmax,bnmax),ww(bnmax,bnmax),yy(bnmax)

C ... calculate inverse from LU decomposition ...

        do 10 j=1,n

C     ... forward substitution ...

          krd = 0.
          if(1 .eq. j) krd = 1.
          yy(1) = krd/aa(1,1)
          do 20 i=2,n
            krd = 0.
            if(i .eq. j) krd = 1.
            sum = 0.
            do 30 k=1,i-1
              sum = sum + aa(i,k)*yy(k)
 30         continue
            yy(i) = (krd - sum)/aa(i,i)
 20       continue

C     ... backward substitution ...

          ww(n,j) = yy(n)
          do 40 i=n-1,1,-1
            sum = 0.
            do 50 k = i+1,n
              sum = sum + aa(i,k)*ww(k,j)
 50         continue
            ww(i,j) = yy(i) - sum
 40       continue

 10     continue

C ... return ...

        return
        end

         
        subroutine ludecomp(n,a)

        implicit none

C  ... variable declarations ...

        integer i,j,k,n

        include 'incs.com'

        real sum,a(bnmax,bnmax)

C ... Crout's LU decomposition ...

C     ... calculate u(1,j) ...

        do 10 j=2,n
          a(1,j) = a(1,j)/a(1,1)
 10     continue

C     ... interior points ...

        do 20 j=2,n-1

C       ... calculate l(i,j) ...

          do 30 i=j,n
            sum = 0.
            do 40 k=1,j-1
              sum = sum + a(i,k)*a(k,j)
 40         continue
            a(i,j) = a(i,j) - sum
 30       continue

C       ... calculate u(j,k) ...

          do 50 k=j+1,n
            sum = 0.
            do 60 i=1,j-1
              sum = sum + a(j,i)*a(i,k)
 60         continue
            a(j,k) = (a(j,k) - sum)/a(j,j)
 50       continue
 20     continue

C     ... calculate l(n,n) ...

        sum = 0.
        do 70 k=1,n-1
          sum = sum + a(n,k)*a(k,n)
 70     continue
        a(n,n) = a(n,n) - sum

C ... return ...

        return
        end
        subroutine blockfbs(ewcf,icap,jcap,p,v,g)

        implicit none

C ... variable declarations ...

        integer i,icap,jcap,m

        include 'incs.com'

        integer p(0:npmax)
        
        real v(0:gptmax),g(0:gptmax)
        real ewcf
        real qq(bnmax),xp(bnmax),tcap(bnmax,bnmax,bnmax)
        real y(bnmax,bnmax),x(bnmax,bnmax)

C ... common statements ...

        common/invblock/tcap

C ... forward substitution ...

        do 10 m=1,jcap
          qq(m) = g(p(1)+m)
 10     continue

        call mulmatx(1,jcap,tcap,qq,y)

        do 20 i=2,icap
          do 25 m=1,jcap
            qq(m) = g(p(i)+m) + ewcf*y(m,i-1)
 25       continue
          call mulmatx(i,jcap,tcap,qq,y)
 20     continue

C ... backward substitution ...

        do 30 m=1,jcap
          xp(m) = y(m,icap)
          v(p(icap)+m) = y(m,icap)
 30     continue

        do 40 i=icap-1,1,-1
          call mulmatx(i,jcap,tcap,xp,x)
          do 45 m=1,jcap
            xp(m) = y(m,i) + ewcf*x(m,i)
            v(p(i)+m) = xp(m)
 45       continue
 40     continue

C ... return ...

        return
        end
        subroutine mulmatx(i,n,tc,q,w)

        implicit none

C ... variable declarations ...

        integer i,n,j,k

        include 'incs.com'

        real sum
        real q(bnmax),w(bnmax,bnmax),tc(bnmax,bnmax,bnmax)

C ... matrix-vector multiply ...

        do 10 j=1,n
          sum = 0.
          do 20 k=1,n
            sum = sum + tc(j,k,i)*q(k)
 20       continue
          w(j,i) = sum
 10     continue

C ... return ...

        return
        end
