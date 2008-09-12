      subroutine bsegp(n, gendat, ostat, pval, ng, tol)
      integer n, ng
      double precision gendat(n), ostat, pval, tol

      double precision btmax, btailp
      external btmax, btailp

      ostat = btmax(n, gendat)
c      call dblepr("Max Stat",8,ostat,1)
      pval = btailp(ostat, n, ng, tol)
      if (pval .gt. 1) pval = 1.0d0

      return
      end

      double precision function btmax(n, x)
      integer n
      double precision x(n)

      integer i
      double precision sumxi, btmaxi, dn, di, ostat

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

c     pseudo confidence interval based on permutations
      subroutine bsegci(n, k, sumxk, x, px, sr, vfact, nperm, bsloc)
      integer n, k, sr(2), nperm, bsloc(nperm)
      double precision sumxk, x(n), px(n), vfact(n)

      integer k1, nk, np, ibseg

      call rndstart()
      k1 = k+1
      nk = n-k
      do 10 np = 1, nperm
         call xperm(k,x,px)
         call xperm(nk,x(k1),px(k1))
         call btmxci(n,k,sr,px,vfact,ibseg,sumxk)
         bsloc(np) = ibseg
 10   continue
      call rndend()

      return
      end

      subroutine btmxci(n,k,sr,x,vfact,ibseg,sumxk)
      integer n,k,sr(2),ibseg
      double precision x(n),vfact(n),sumxk

      integer i
      double precision sumxi, ostat, btmaxi

      ostat = vfact(k)*(sumxk**2)
      ibseg = k
      sumxi = sumxk
      do 10 i = k-1,sr(1),-1
         sumxi = sumxi - x(i+1)
         btmaxi = vfact(i)*(sumxi**2)
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            ibseg = i
         endif
 10   continue

      sumxi = sumxk
      do 20 i = k+1,sr(2)
         sumxi = sumxi + x(i)
         btmaxi = vfact(i)*(sumxi**2)
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            ibseg = i
         endif
 20   continue
      ostat = sqrt(ostat)

      return
      end
