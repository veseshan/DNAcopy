c     Ternary segmentation with permutation reference distribution
      subroutine fndcpt(n,twon,x,tss,px,sx,tx,nperm,cpval,ncpt,
     1     icpt,ibin,hybrid,hk,delta,ngrid,sbn,sbdry,tol)
      integer n,twon,nperm,ncpt,icpt(2),hk,ngrid,sbn,sbdry(sbn)
      logical ibin,hybrid
      double precision x(n),tss,px(n),sx(n),tx(twon),cpval,delta,tol

      integer np,nrej,nrejc,iseg(2),n1,n2,n12,l,k
      double precision ostat,ostat1,pstat,tpval,pval1,pval2

      double precision tailp, tmax, htmax, tpermp
      external tailp, tmax, htmax, tpermp

      call rndstart()

      nrej = 0
      ncpt = 0

      call tmax1(n,twon,x,tss,sx,tx,iseg,ostat,ibin)
      ostat1 = sqrt(ostat)
      ostat = ostat * 0.99999
c      call dblepr("Max Stat",8,ostat,1)
c      call intpr("Location",8,iseg,2)

c     if maximal t-statistic is too small (for now use 0.1) don't split
      if (ostat1 .le. 0.1) go to 500
c     if maximal t-statistic is too large (for now use 7.0) split
      if (ostat1 .ge. 7.0) go to 200
c     o.w calculate p-value and decide if & how data are segmented
      if (hybrid) then
         pval1 = tailp(ostat1, delta, n, ngrid, tol)
         if (pval1 .gt. cpval) go to 500
         pval2 = cpval - pval1
         nrejc = int(pval2*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 50 np = 1,nperm
            call xperm(n,x,px)
            pstat = htmax(n,twon,hk,tss,px,sx,tx,ibin)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
            if (nrej.gt.nrejc) go to 500
            if (np .ge. sbdry(k)) go to 200
 50      continue
      else
         nrejc = int(cpval*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 100 np = 1,nperm
            call xperm(n,x,px)
            pstat = tmax(n,twon,tss,px,sx,tx,ibin)
c     call dblepr("Perm Max Stat",13,pstat,1)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
c     call intpr("num rej",7,nrej,1)
            if (nrej.gt.nrejc) go to 500
            if (np .ge. sbdry(k)) go to 200
 100     continue
      endif
 200  if (iseg(2).eq.n) then
         ncpt = 1
         icpt(1) = iseg(1)
      else
         if(iseg(1).eq.0) then
            ncpt = 1
            icpt(1) = iseg(2)
         else
            l = 1
            n1 = iseg(1)
            n12 = iseg(2)
            n2 = n12 - n1
            tpval = tpermp(n1,n2,n12,x(l),px,nperm)
c            call dblepr("binseg p-value",14,tpval,1)
            if (tpval.le.cpval) then
               ncpt = 1
               icpt(1) = iseg(1)
            endif
            l = iseg(1) + 1
            n12 = n - iseg(1)
            n2 = n - iseg(2)
            n1 = n12 - n2
            tpval = tpermp(n1,n2,n12,x(l),px,nperm)
c            call dblepr("binseg p-value",14,tpval,1)
            if (tpval.le.cpval) then
               ncpt = ncpt + 1
               icpt(ncpt) = iseg(2)
            endif
         endif
      endif

 500  call rndend()

      return
      end

      double precision function tmax(n,twon,tss,px,sx,tx,ibin)
      integer n,twon
      double precision tss,px(n),sx(n),tx(twon)
      logical ibin

      integer i,j,k,l
c      double precision xsum,sx2,x1,x2,rj,rn,tij,xvar,
      double precision sumx,ssqx,xbar,rn,rj,absx,sxmx,bssmx,wtmax

      rn = dfloat(n)
      do 20 i = 1,n
         tx(i) = px(i)
         tx(n+i) = tx(i)
         sx(i) = tx(i)
 20   continue
      tmax = 0.0
c     compute the max statistic for segments of length j
      do 40 j = 2,(n-1)/2
         rj = dfloat(j)
         sxmx = 0.0
         do 30 i = 1,n
            sx(i) = sx(i) + tx(i+j-1)
            absx = abs(sx(i))
            if (sxmx.lt.absx) sxmx = absx
 30      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (tmax.lt.bssmx) tmax = bssmx
 40   continue
c     compute the max statistic for segments of length n/2 (if integer)
      if (n.eq.2*(n/2)) then
         j = n/2
         rj = dfloat(j)
         sxmx = 0.0
         do 50 i = 1,n/2
            sx(i) = sx(i) + tx(i+j-1)
            absx = abs(sx(i))
            if (sxmx.lt.absx) sxmx = absx
 50      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (tmax.lt.bssmx) tmax = bssmx
      endif
      if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         tmax = tmax/(tss/rn)
      else
         if (tss.le.tmax+0.0001) tss = tmax + 1.0
         tmax = tmax/((tss-tmax)/(rn-2.0))
      endif

      return
      end

      subroutine tmax1(n,twon,x,tss,sx,tx,iseg,ostat,ibin)
      integer n,twon,iseg(2)
      double precision x(n),tss,sx(n),tx(twon),ostat
      logical ibin

      integer i,j,k,l, sxmxi, tmaxi, tmaxj
c      double precision xsum,sx2,x1,x2,rij,rw,tij,xvar
      double precision rn,rj,absx,sxmx,bssmx,tmax

      rn = dfloat(n)
      ostat = -0.5
      do 20 i = 1,n
         tx(i) = x(i)
         tx(n+i) = tx(i)
         sx(i) = tx(i)
 20   continue
      tmax = -0.5
c     compute the max statistic for segments of length j
      do 40 j = 2,(n-1)/2
         rj = dfloat(j)
c     sxmx is the maximum of abs partial sums x[i] + ... + x[i+j-1]
         sxmx = -0.5
         do 30 i = 1,n
            sx(i) = sx(i) + tx(i+j-1)
            absx = abs(sx(i))
            if (sxmx.lt.absx) then 
               sxmx = absx
c     i is the start point for the current maximum
               sxmxi = i
            endif
 30      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (tmax.lt.bssmx) then 
            tmax = bssmx
            tmaxi = sxmxi
            tmaxj = j
         endif
 40   continue
c     compute the max statistic for segments of length n/2 (if integer)
      if (n.eq.2*(n/2)) then
         j = n/2
         rj = dfloat(j)
         sxmx = -0.5
         do 50 i = 1,n/2
            sx(i) = sx(i) + tx(i+j-1)
            absx = abs(sx(i))
            if (sxmx.lt.absx) then 
               sxmx = absx
               sxmxi = i
            endif
 50      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (tmax.lt.bssmx) then
            tmax = bssmx
            tmaxi = sxmxi
            tmaxj = j
         endif
      endif
c     convert statistic to t^2 form
      if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         tmax = tmax/(tss/rn)
      else
         if (tss.le.tmax+0.0001) tss = tmax + 1.0
         tmax = tmax/((tss-tmax)/(rn-2.0))
      endif

      ostat = tmax
      if (tmaxi+tmaxj-1.le.n) then
         iseg(1) = tmaxi-1
         iseg(2) = tmaxi-1+tmaxj
      else
         iseg(1) = tmaxi-1+tmaxj-n
         iseg(2) = tmaxi-1
      endif

      return
      end

      subroutine xperm(n,x,px)
      integer n
      double precision x(n),px(n)

      integer i,j
      double precision cc,tmpx

      double precision dunif
      external dunif

      do 10 i = 1,n
         px(i) = x(i)
 10   continue

      do 20 i = n,1,-1
         cc = dunif()
         j = int(cc*dfloat(i))+1
         tmpx = px(i)
         px(i) = px(j)
         px(j) = tmpx
 20   continue
      return
      end

      double precision function tpermp(n1,n2,n,x,px,nperm)
      integer n1,n2,n,nperm
      double precision x(n),px(n)
      integer np,i,m1
      double precision xsum1,xsum2,xbar,ostat,pstat,rn1,rn2,rm1

      rn1 = dfloat(n1)
      rn2 = dfloat(n2)
      if (n1.eq.1 .or. n2.eq.1) then
         nrej = nperm
         go to 110
      endif
      xsum1 = 0.0
      do 10 i=1,n1
         xsum1 = xsum1 + x(i)
 10   continue
      xsum2 = 0.0
      do 20 i=n1+1,n
         xsum2 = xsum2 + x(i)
 20   continue
      xbar = (xsum1 + xsum2)/(rn1+rn2)
      if (n1.le.n2) then
         m1 = n1
         rm1 = rn1
         ostat = 0.99999*abs(xsum1/rn1 - xbar)
      else
         m1 = n2
         rm1 = rn2
         ostat = 0.99999*abs(xsum2/rn2 - xbar)
      endif
c      call dblepr("O-Stat",6,ostat,1)
      nrej = 0
      do 100 np = 1,nperm
         call xperm(n,x,px)
         xsum1 = 0.0
         do 30 i=1,m1
            xsum1 = xsum1 + px(i)
 30      continue
         pstat = abs(xsum1/rm1 - xbar)
c         call dblepr("P-Stat",6,pstat,1)
         if (ostat.le.pstat) nrej = nrej + 1
 100  continue
 110  tpermp = dfloat(nrej)/dfloat(nperm)

      return
      end
