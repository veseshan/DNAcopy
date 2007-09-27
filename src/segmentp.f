      subroutine bsegp(n, gendat, ostat, pval, ng, tol)
      integer n, ng
      double precision gendat(n), ostat, pval, tol

      double precision btmax, btailp
      external btmax, btailp

      ostat = btmax(n, gendat)
c      call dblepr("Max Stat",8,ostat,1)
      pval = btailp(ostat, n, ng, tol)

      return
      end

      double precision function btmax(n, x)
      integer n
      double precision x(n)

      double precision sumxi, btmaxi, dn, di

      sumxi = x(1)
      ostat = 0.0
      dn = dfloat(n)
      di = 1.0
      do 20 i = 2,n-2
         di = di + 1.0
         sumxi = sumxi + x(i)
         btmaxi = dn*(sumxi**2)/(di*(dn-di))
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
c            ibseg = i
         endif
 20   continue
      btmax = sqrt(ostat)

      return
      end

c     tail probability of binary segmentation statistic
c     from page 387 of Siegmund (1986) paper
      double precision function btailp(b, m, ng, tol)
      integer m, ng
      double precision b, tol

      double precision ll, ul, dincr, nulo, nuhi, x, x1, dm
      integer i, k

      double precision fpnorm, nu
      external fpnorm, nu

      dm = dfloat(m)
      k = 2
      ll = b*sqrt(1.0/dfloat(m-k) - 1.0/dfloat(m))
      ul = b*sqrt(1.0/dfloat(k) - 1.0/dfloat(m))
      dincr = (ul - ll)/dfloat(ng)

      btailp = 0.0
      x = ll
      x1 = x + (b**2)/(dm*x)
      nulo = nu(x1, tol)/x
      do 10 i = 1, ng
         x = x + dincr
         x1 = x + (b**2)/(dm*x)
         nuhi = nu(x1, tol)/x
         btailp = btailp + (nuhi + nulo)*dincr
         nulo = nuhi
 10   continue
      btailp = b*exp(-b**2/2)*btailp/2.506628275

      btailp =  btailp + 2*(1.0-fpnorm(b))

      return
      end
