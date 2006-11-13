      double precision function htmax(n,n2,k,px,sx,tx,ibin)
      integer n,n2,k
      double precision px(n),sx(n),tx(n2)
      logical ibin

      integer i,j
c      double precision xsum,sx2,x1,x2,rj,rn,tij,xvar,
      double precision sumx,ssqx,xbar,rn,rj,absx,sxmx,bssmx,tss

      rn = dfloat(n)
      sumx = 0.0d0
      ssqx = 0.0d0
      do 10 i = 1,n
         sumx = sumx + px(i)
         ssqx = ssqx + px(i)**2
 10   continue
      xbar = sumx/rn
      tss = ssqx - sumx*xbar
      do 20 i = 1,n
         tx(i) = px(i) - xbar
         tx(n+i) = tx(i)
         sx(i) = tx(i)
 20   continue
      htmax = 0.0d0
      do 40 j = 2,k
         rj = dfloat(j)
         sxmx = 0.0d0
         do 30 i = 1,n
            sx(i) = sx(i) + tx(i+j-1)
            absx = abs(sx(i))
            if (sxmx.lt.absx) sxmx = absx
 30      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5d0)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (htmax.lt.bssmx) htmax = bssmx
 40   continue
      if (ibin) then
         if (tss .le. 0.0001d0) tss = 1.0d0
         htmax = htmax/(tss/rn)
      else
         if (tss .le. htmax+0.0001d0) tss = htmax + 1.0d0
         htmax = htmax/((tss-htmax)/(rn-2.0d0))
      endif

      return
      end

      double precision function tailp(b, delta, m, ngrid, tol)
      double precision b, delta, tol
      integer m, ngrid
c     it1tsq is the integral of 1/(t*(1-t))**2
      double precision nu, it1tsq
      external nu, it1tsq

      double precision t, tl, dincr, bsqrtm, x, nux
      integer i

      dincr = (0.5d0 - delta)/dfloat(ngrid)
      bsqrtm = b/sqrt(dfloat(m))

      tl = 0.5d0 - dincr
      t = 0.5d0 - 0.5d0*dincr
      tailp = 0.0d0
      do 10 i = 1,ngrid
         tl = tl + dincr
         t = t + dincr
         x = bsqrtm/sqrt(t*(1-t))
         nux = nu(x, tol)
         tailp = tailp + (nux**2)*it1tsq(tl, dincr)
 10   continue
      tailp = 9.973557d-2*(b**3)*exp(-b**2/2)*tailp
c     since test is two-sided need to multiply tailp by 2
      tailp = 2.0d0*tailp

      return
      end

c     integral of 1/(t*(1-t))**2 from x to x+a
      double precision function it1tsq(x, a)
      double precision x, a

      double precision y

      y = x + a - 0.5d0
      it1tsq = (8.0d0*y)/(1.0d0 - 4.0d0*y**2) + 
     1     2.0d0*log((1.0d0 + 2.0d0*y)/(1.0d0 - 2.0d0*y))
      y = x - 0.5d0
      it1tsq = it1tsq - (8.0d0*y)/(1.0d0 - 4.0d0*y**2) -
     1     2.0d0*log((1.0d0 + 2.0d0*y)/(1.0d0 - 2.0d0*y))

      return
      end

      double precision function nu(x, tol)
      double precision x, tol

      double precision fpnorm
      external fpnorm

      double precision lnu0, lnu1, dk, xk
      integer i, k

      if (x .gt. 0.01d0) then
         lnu1 = log(2.0d0) - 2*log(x)
         lnu0 = lnu1
         k = 2
         dk = 0
         do 10 i = 1, k
            dk = dk + 1
            xk = -x*sqrt(dk)/2.0d0
            lnu1 = lnu1 - 2.0d0*fpnorm(xk)/dk
 10      continue

         do 50 while (dabs((lnu1-lnu0)/lnu1) .gt. tol)
            lnu0 = lnu1
            do 20 i = 1,k
               dk = dk + 1
               xk = -x*sqrt(dk)/2.0d0
               lnu1 = lnu1 - 2.0d0*fpnorm(xk)/dk
 20         continue
            k = 2*k
 50      enddo
      else
         lnu1 = -0.583d0*x
      endif
      nu = exp(lnu1)

      return
      end
