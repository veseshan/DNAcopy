c     Ternary segmentation with permutation reference distribution
      subroutine fndcpt(n,x,tss,px,sx,nperm,cpval,ncpt,icpt,ibin,
     1     hybrid,al0,hk,delta,ngrid,sbn,sbdry,tol)
      integer n,nperm,ncpt,icpt(2),al0,hk,ngrid,sbn,sbdry(sbn)
      logical ibin,hybrid
      double precision x(n),tss,px(n),sx(n),cpval,delta,tol

      integer np,nrej,nrejc,iseg(2),n1,n2,n12,l,k
      double precision ostat,ostat1,pstat,tpval,pval1,pval2

c     new functions to replace tmax and htmax (also tmaxo replaces tmax1)
      double precision tailp, tmaxp, htmaxp, tpermp
      external tailp, tmaxp, htmaxp, tpermp

      call rndstart()

      nrej = 0
      ncpt = 0

c      call tmax1(n,twon,x,tss,sx,tx,iseg,ostat,ibin)
      call tmaxo(n,x,tss,sx,iseg,ostat,al0,ibin)
      ostat1 = sqrt(ostat)
      ostat = ostat * 0.99999
c      call dblepr("Max Stat",8,ostat,1)
c      call intpr("Location",8,iseg,2)

c     if maximal t-statistic is too small (for now use 0.1) don't split
      if (ostat1 .le. 0.1) go to 500
c     if maximal t-statistic is too large (for now use 7.0) split
c     also make sure it's not affected by outliers i.e. small seglength
      l = min(iseg(2) - iseg(1), n - iseg(2) + iseg(1))
      if ((ostat1 .ge. 7.0) .and. (l .ge. 10)) go to 200
c     o.w calculate p-value and decide if & how data are segmented
      if (hybrid) then
         pval1 = tailp(ostat1, delta, n, ngrid, tol)
         if (pval1 .gt. cpval) go to 500
         pval2 = cpval - pval1
         nrejc = int(pval2*dfloat(nperm))
         k=nrejc*(nrejc+1)/2 + 1
         do 50 np = 1,nperm
            call xperm(n,x,px)
c            pstat = htmax(n,twon,hk,tss,px,sx,tx,ibin)
            pstat = htmaxp(n,hk,tss,px,sx,al0,ibin)
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
c            pstat = tmax(n,twon,tss,px,sx,tx,ibin)
            pstat = tmaxp(n,tss,px,sx,al0,ibin)
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

c     code to permute the data vector
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

c     function for the p-value of t-statistics for removing edge effects
      double precision function tpermp(n1,n2,n,x,px,nperm)
      integer n1,n2,n,nperm
      double precision x(n),px(n)

      integer np,i,m1,j,nrej
      double precision xsum1,xsum2,xbar,ostat,pstat,rn1,rn2,rm1,
     1     tstat, tss, rn, cc, tmpx

      double precision dunif
      external dunif

      rn1 = dfloat(n1)
      rn2 = dfloat(n2)
      rn = rn1 + rn2
      if (n1.eq.1 .or. n2.eq.1) then
         nrej = nperm
         go to 110
      endif
      xsum1 = 0.0
      tss = 0.0
      do 10 i=1,n1
         px(i) = x(i)
         xsum1 = xsum1 + x(i)
         tss = tss + x(i)**2
 10   continue
      xsum2 = 0.0
      do 20 i=n1+1,n
         px(i) = x(i)
         xsum2 = xsum2 + x(i)
         tss = tss + x(i)**2
 20   continue
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
c      call dblepr("O-Stat",6,ostat,1)
      nrej = 0
      tstat = tstat/((tss-tstat)/(rn-2.0))
c      call dblepr("T-square",8,tstat,1)
c     if observed t is large (> 5) don't bother with permutation p-value
c     also make sure there are enough observations i.e. m1 >= 10
      if ((tstat .gt. 25) .and. (m1 .ge. 10)) go to 110
      do 100 np = 1,nperm
c*******************************************
c     the following is very inefficient
c*******************************************
c         call xperm(n,x,px)
c         xsum1 = 0.0
c         do 30 i=1,m1
c            xsum1 = xsum1 + px(i)
c 30      continue
c*******************************************
c     changed to the following: instead of
c     full permutation sample m1 w.o. repl
c******************************************* 
         xsum1 = 0
         do 30 i = n,n-m1+1,-1
            cc = dunif()
            j = int(cc*dfloat(i))+1
            tmpx = px(i)
            px(i) = px(j)
            px(j) = tmpx
            xsum1 = xsum1 + px(i)
 30      continue
         pstat = abs(xsum1/rm1 - xbar)
c         call dblepr("P-Stat",6,pstat,1)
         if (ostat.le.pstat) nrej = nrej + 1
 100  continue
 110  tpermp = dfloat(nrej)/dfloat(nperm)

      return
      end
