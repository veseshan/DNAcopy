c     for binary segmentation of exon data re-use the segment.p code 
c     with one-sided p-value
      subroutine esegp(n, exndat, ostat, eloc, pval, ng, tol)
      integer n, eloc, ng
      double precision exndat(n), ostat, pval, tol

      double precision btailp
      external btailp

      integer i
      double precision tss

      tss = 0
      do 10 i = 1, n
         tss = tss + exndat(i)**2
 10   continue
      call etmax(n, exndat, tss, ostat, eloc)
c      call dblepr("Max Stat",8,ostat,1)
      pval = btailp(ostat, n, ng, tol)
      pval = pval/2
      if (pval .gt. 1) pval = 1.0d0

      return
      end

c     looking for increase in expression - use a one-sided t-statistic
c     t_i = S_i*sqrt(n/(i*(n-i)))/sqrt((tss - S_i^2*n/(i*(n-i)))/(n-2))
      subroutine etmax(n, x, tss, ostat, eloc)
      integer n, eloc
      double precision x(n), tss, ostat

      integer i
      double precision sumxi, btmaxi, dn, di

      sumxi = x(1)
      ostat = 0.0
      eloc = -1
      dn = dfloat(n)
      di = 1.0
      do 20 i = 2,n-2
         di = di + 1.0
         sumxi = sumxi + x(i)
         btmaxi = -sumxi/sqrt(di*(dn-di))
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            eloc = i
         endif
 20   continue
      ostat = (ostat/sqrt(tss - dn*ostat**2))*sqrt(dn*(dn-2))

      return
      end
