      subroutine getbdry(eta, m, nperm, mb, ibdry, etastr, tol)
      integer m, nperm, mb, ibdry(mb)
      double precision eta, etastr(m), tol

      double precision eta0, etalo, etahi, plo, phi, pexcd
      integer j, l

      l = 1
      ibdry(1) = nperm-int(dfloat(nperm)*eta)
      etastr(1) = eta
      eta0 = eta
      do 20 j = 2,m
         etahi = eta0*1.1
         call etabdry(nperm, etahi, j, ibdry(l+1))
         call pexceed(nperm, j, ibdry(l+1), phi)
         etalo = eta0*0.25
         call etabdry(nperm, etalo, j, ibdry(l+1))
         call pexceed(nperm, j, ibdry(l+1), plo)
         do 10 while ((etahi-etalo)/etalo .gt. tol)
            eta0 = etalo + (etahi-etalo)*(eta-plo)/(phi-plo)
            call etabdry(nperm, eta0, j, ibdry(l+1))
            call pexceed(nperm, j, ibdry(l+1), pexcd)
            if (pexcd .gt. eta) then
               etahi = eta0
               phi = pexcd
            else
               etalo = eta0
               plo = pexcd
            endif
 10      continue
         etastr(j) = eta0
         l = l+j
 20   continue

      return
      end

      subroutine etabdry(nperm, eta0, n1s, ibdry)
      integer nperm, n1s, ibdry(n1s)
      double precision eta0

      double precision fphypr
      external fphypr

      integer i, k
      double precision di, dn, dn1s, dk, tprob

      dn1s = dfloat(n1s)
      dn = dfloat(nperm-n1s)
      
      k = 0
      dk = 0.0d0
      do 10 i = 1, nperm
         di = dfloat(i)
         tprob = fphypr(dk, dn1s, dn, di)
         if (tprob .le. eta0) then
            k = k+1
            dk = dk + 1.0d0
            ibdry(k) = i
         endif
 10   continue

      return
      end

      subroutine pexceed(nperm, n1s, ibdry, pexcd)
      integer nperm, n1s, ibdry(n1s)
      double precision pexcd

      double precision dn, dk, dn1, dk1, dn2, dk2, dn3, dk3, dlcnk
      integer i

      double precision flchoose
      external flchoose

      dn = dfloat(nperm)
      dk = dfloat(n1s)
      dn1 = dfloat(nperm-ibdry(1))
      dlcnk = flchoose(dn, dk)

      pexcd = exp(flchoose(dn1, dk) - dlcnk)

      if (n1s .ge. 2) then
         dn1 = dfloat(ibdry(1))
         dn = dfloat(nperm-ibdry(2))
         dk = dfloat(n1s-1)
         pexcd = pexcd + exp(log(dn1) + flchoose(dn, dk) - dlcnk)
      endif

      if (n1s .ge. 3) then
         dn1 = dfloat(ibdry(1))
         dn2 = dfloat(ibdry(2))
         dn = dfloat(nperm-ibdry(3))
         dk = dfloat(n1s-2)
         pexcd = pexcd + 
     1        exp(log(dn1) + log(dn1-1.0) - log(2.0) + 
     2                        flchoose(dn, dk) - dlcnk) +
     3        exp(log(dn1) + log(dn2-dn1) + flchoose(dn, dk) - dlcnk)
      endif

      if (n1s .gt. 3) then
         do 10 i = 4, n1s
            dn1 = dfloat(ibdry(i-3))
            dk1 = dfloat(i-1)
            dk2 = dfloat(i-2)
            dk3 = dfloat(i-3)
            dn2 = dfloat(ibdry(i-2))
            dn3 = dfloat(ibdry(i-1))
            dn = dfloat(nperm-ibdry(i))
            dk = dfloat(n1s-i+1)
            pexcd = pexcd + 
     1           exp(flchoose(dn1, dk1) + flchoose(dn, dk) - dlcnk) +
     2           exp(flchoose(dn1, dk2) + log(dn3-dn1) + 
     3                       flchoose(dn, dk) - dlcnk) +
     4           exp(flchoose(dn1, dk3) + log(dn2-dn1) + log(dn3-dn2) +
     3                        flchoose(dn, dk) - dlcnk) +
     5           exp(flchoose(dn1, dk3) + log(dn2-dn1) - log(2.0) + 
     6                    log(dn2-dn1-1.0) + flchoose(dn, dk) - dlcnk)
 10      continue
      endif

      return
      end
