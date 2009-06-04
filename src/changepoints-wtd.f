c     Ternary segmentation with permutation reference distribution
c     probes have weights due to differences in variances

      subroutine wfindcpt(n,x,tss,wts,rwts,cwts,px,sx,nperm,cpval,ncpt,
     1     icpt,hybrid,al0,hk,mncwt,delta,ngrid,sbn,sbdry,tol)
      integer n,nperm,ncpt,icpt(2),al0,hk,ngrid,sbn,sbdry(sbn)
      logical hybrid
      double precision x(n),tss,wts(n),rwts(n),cwts(n),px(n),sx(n),
     1     cpval,mncwt(hk),delta,tol

      integer np,nrej,nrejc,iseg(2),n1,n2,n12,l,k
      double precision ostat,ostat1,pstat,tpval,pval1,pval2

c     new functions to replace tmax and htmax (also tmaxo replaces tmax1)
      double precision tailp, wtmaxp, hwtmaxp, wtpermp
      external tailp, wtmaxp, hwtmaxp, wtpermp

      call rndstart()

      nrej = 0
      ncpt = 0

c     call the observed statistic routine
      call wtmaxo(n,x,wts,tss,sx,cwts,iseg,ostat,al0)
      ostat1 = sqrt(ostat)
      ostat = ostat * 0.99999

c     if maximal t-statistic is too small (for now use 0.1) don't split
      if (ostat1 .le. 0.1) go to 500
c     if maximal t-statistic is too large (for now use 7.0) split
c     also make sure it's not affected by outliers i.e. small seglength
      l = min(iseg(2) - iseg(1), n - iseg(2) + iseg(1))
      if ((ostat1 .ge. 7.0) .and. (l .ge. 10)) go to 200
c     o.w calculate p-value and decide if & how data are segmented
      if (hybrid) then
         call getmncwt(n, cwts, hk, mncwt, delta)
c     delta is a function of arc lengths
         pval1 = tailp(ostat1, delta, n, ngrid, tol)
         if (pval1 .gt. cpval) go to 500
         pval2 = cpval - pval1
         nrejc = int(pval2*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 50 np = 1,nperm
c     call permutation code for data with weights
            call wxperm(n,x,px,rwts)
c     call the small arc permutation statistic function
            pstat = hwtmaxp(n,hk,px,wts,sx,cwts,mncwt,al0)
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
c     call permutation code for data with weights
            call wxperm(n,x,px,rwts)
c     call full data permutation statistic function
            pstat = wtmaxp(n,px,wts,sx,cwts,al0)
            if (ostat.le.pstat) then
               nrej = nrej + 1
               k = k + 1
            endif
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
            tpval = wtpermp(n1,n2,n12,x(l),px,wts(l),rwts(l),nperm)
            if (tpval.le.cpval) then
               ncpt = 1
               icpt(1) = iseg(1)
            endif
            l = iseg(1) + 1
            n12 = n - iseg(1)
            n2 = n - iseg(2)
            n1 = n12 - n2
            tpval = wtpermp(n1,n2,n12,x(l),px,wts(l),rwts(l),nperm)
            if (tpval.le.cpval) then
               ncpt = ncpt + 1
               icpt(ncpt) = iseg(2)
            endif
         endif
      endif

 500  call rndend()

      return
      end

c     *******  code to permute the data vector with weights  ********
c     since variance of probe i is inversely proportional to weights 
c     multiply by square root, permute and then divide by square root
      subroutine wxperm(n,x,px,rwts)
      integer n
      double precision x(n),px(n),rwts(n)

      integer i,j
      double precision cc,tmpx

      double precision dunif
      external dunif

      do 10 i = 1,n
         px(i) = x(i)*rwts(i)
 10   continue

      do 20 i = n,1,-1
         cc = dunif()
         j = int(cc*dfloat(i))+1
         tmpx = px(i)
         px(i) = px(j)/rwts(i)
         px(j) = tmpx
 20   continue

      return
      end

c     function for the p-value of t-statistics for removing edge effects
      double precision function wtpermp(n1,n2,n,x,px,wts,rwts,nperm)
      integer n1,n2,n,nperm
      double precision x(n),px(n),wts(n),rwts(n)

      integer np,i,m1,j,nrej
      double precision xsum1,xsum2,xbar,ostat,pstat,rn1,rn2,rm1,
     1     tstat, tss, rn, cc, tmpx

      double precision dunif
      external dunif

      if (n1.eq.1 .or. n2.eq.1) then
         nrej = nperm
         go to 110
      endif
      xsum1 = 0.0
      tss = 0.0
      rn1 = 0.0
      do 10 i=1,n1
         px(i) = x(i)*rwts(i)
         xsum1 = xsum1 + wts(i)*x(i)
         tss = tss + wts(i)*x(i)**2
         rn1 = rn1 + wts(i)
 10   continue
      xsum2 = 0.0
      rn2 = 0.0
      do 20 i=n1+1,n
         px(i) = x(i)
         xsum2 = xsum2 + wts(i)*x(i)
         tss = tss + wts(i)*x(i)**2
         rn2 = rn2 + wts(i)
 20   continue
      rn = rn1 + rn2
      xbar = (xsum1 + xsum2)/rn
      tss = tss - rn*(xbar**2)
      if (n1.le.n2) then
         m1 = n1
         rm1 = rn1
         ostat = 0.99999*abs(xsum1/rn1 - xbar)
         tstat = (ostat**2)*rn1*rn/rn2
      else
         m1 = n2
         rm1 = rn2
         ostat = 0.99999*abs(xsum2/rn2 - xbar)
         tstat = (ostat**2)*rn2*rn/rn1
      endif
      nrej = 0
      tstat = tstat/((tss-tstat)/(dfloat(n)-2.0))
c     if observed t is large (> 5) don't bother with permutation p-value
c     also make sure there are enough observations i.e. m1 >= 10
      if ((tstat .gt. 25) .and. (m1 .ge. 10)) go to 110
      do 100 np = 1,nperm
         xsum1 = 0
         do 30 i = n,n-m1+1,-1
            cc = dunif()
            j = int(cc*dfloat(i))+1
            tmpx = px(i)
            px(i) = px(j)
            px(j) = tmpx
c     the observation should be divided by sqrt(wts(i)) to get the correct 
c     probe variance.  But should be multiplied by wts(i) for statistic
            xsum1 = xsum1 + px(i)*rwts(i)
 30      continue
         pstat = abs(xsum1/rm1 - xbar)
         if (ostat.le.pstat) nrej = nrej + 1
 100  continue
 110  wtpermp = dfloat(nrej)/dfloat(nperm)

      return
      end
