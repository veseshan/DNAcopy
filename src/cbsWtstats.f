c     these are the subroutines to do the weigthed version of CBS
c     which is useful in order to merge data from multiple platforms
c     --------------------------------------------------------------
c         This is relevant only for log-ratio not binary data
c     --------------------------------------------------------------
c     function for calculating the full max weighted t-statistic
c     new approach to maximizing t-statistic

      subroutine wtmaxo(n,x,wts,tss,sx,cwts,iseg,ostat,al0)
      integer n,iseg(2),al0
      double precision x(n),wts(n),tss,sx(n),cwts(n),ostat
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, ihi1, jlo1, jhi1,
     2     tmaxi, tmaxj, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, sij1, sij2, sijmx0, bijbss, awtmax, psrnov2,
     2     psdiff, psrj, psrn, psrnj, awtlo, awthi, awt1
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j
      double precision, allocatable :: bssbij(:), bssijmax(:), awt(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:)

c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
      if (n .ge. 50) then
         nb = nint(sqrt(dfloat(n)))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb), ibmin(nb), ibmax(nb))
      allocate(bssbij(nb2), bssijmax(nb2), awt(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2))

c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      psmin0 = 0
      psmax0 = 0
      ipsmin0 = n
      ipsmax0 = n
      do 20 j = 1, nb
         sx(ilo) = psum + x(ilo)*wts(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + x(i)*wts(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         ibmin(j) = ipsmin
         ibmax(j) = ipsmax
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     adjust global min, max and locations
         if (psmin .lt. psmin0) then
            psmin0 = psmin
            ipsmin0 = ipsmin
         endif
         if (psmax .gt. psmax0) then
            psmax0 = psmax
            ipsmax0 = ipsmax
         endif
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
 20   continue

c     calculate bss for max s_i - min s_i
      psdiff = psmax0 - psmin0
c     if the segment is all constant then psdiff = 0 and so bssmax = 0
      if (psdiff .le. 0) then
         bssmax = 0
         go to 120
      endif
      psrn = cwts(n)
      psrj = abs(cwts(ipsmax0) - cwts(ipsmin0))
      psrnj = psrj*(psrn-psrj)
      bssmax = (psdiff**2)/psrnj
      tmaxi = min(ipsmax0, ipsmin0)
      tmaxj = max(ipsmax0, ipsmin0)

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      psrnov2 = psrn/2
      l = 0
      nal0 = n - al0
      do 40 i = 1, nb
         do 30 j = i, nb
c     calculate bsslim
            if (i .eq. 1) then
               ilo = 1
            else
               ilo = bb(i-1) + 1
            endif
            ihi = bb(i)
            if (j .eq. 1) then
               jlo = 1
            else
               jlo = bb(j-1) + 1
            endif
            jhi = bb(j)
c     for wCBS calculated hi and lo arc weights instead of lengths
            awthi = cwts(jhi) - cwts(ilo)
            if (jhi - ilo .gt. nal0) then
               awthi = 0
               do 35 k = 1, al0
                  awthi = max(awthi, cwts(nal0+k) - cwts(k))
 35            continue
            endif
            if (i .eq. j) then
               awtlo = cwts(ilo+al0) - cwts(ilo)
               do 36 k = ilo + 1, ihi - al0
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 36            continue
            else if (i+1 .eq. j) then
               awtlo = cwts(jlo) - cwts(jlo-al0)
               do 37 k = jlo - al0 + 1, ihi
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 37            continue
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            psrnj = min(awtlo*(psrn-awtlo), awthi*(psrn-awthi))
            bsslim = (sijmx0**2)/psrnj
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  awt(l) = abs(cwts(ibmax(j)) - cwts(ibmin(i)))
                  bssbij(l) = (sij1**2)/(awt(l)*(psrn-awt(l)))
               else
                  awt(l) = abs(cwts(ibmin(j)) - cwts(ibmax(i)))
                  bssbij(l) = (sij2**2)/(awt(l)*(psrn-awt(l)))
               endif
            endif
 30      continue
 40   continue

      nb1 = l

c     Now sort the t-statistics by their magnitude
      call qsort4(bssbij, loc, 1, nb1)

c     now go through the blocks in reverse order (largest down)
      do 100 l = nb1, 1, -1
         k = loc(l)
c     need to check a block only if it has potential to increase bss
c     rjlo is the smalllest (j-i) in the block and rjhi is the largest
         bsslim = bssijmax(k)
         if (bssmax .le. bsslim) then
c     bi, bj give the block location
            bi = bloci(k)
            bj = blocj(k)
            awtmax = awt(k)
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1) + 1
            endif
            ihi = bb(bi)
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
            awthi = cwts(jhi) - cwts(ilo)
            if (bi .eq. bj) then
               awtlo = 0
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c
c     if arc wt is larger than half total wt (psrn/2) make is psrn - arc wt
c
            if (awtmax .gt. psrn - awtmax) awtmax = psrn - awtmax
c
c     if awtlo <= psrn/2 start from (ihi, jlo) and go up
c     if awthi >= psrn/2 start from (ilo, jhi) and go down
c
            if (awtlo .le. psrnov2) then
               if (bi .eq.bj) then 
                  ihi1 = ihi - al0
               else
                  ihi1 = ihi
               endif
               do 60 i = ihi1, ilo, -1
                  jlo1 = max(i + al0, jlo)
                  do 55 j = jlo1, jhi
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .le. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) then
                           bssmax = bijbss
                           tmaxi = i
                           tmaxj = j
                        endif
                     endif
 55               continue
 60            continue
            endif
c
c     make arc wt  psrn - arc wt
c
            awtmax = psrn - awtmax
            if (awthi .ge. psrnov2) then
               do 70 i = ilo, ihi
                  if ((bi .eq. 1) .and. (bj .eq. nb)) then 
                     jhi1 = min(jhi, jhi - al0 + i)
                  else
                     jhi1 = jhi
                  endif
                  do 65 j = jhi1, jlo, -1
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .ge. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) then
                           bssmax = bijbss
                           tmaxi = i
                           tmaxj = j
                        endif
                     endif
 65               continue
 70            continue
            endif
         endif
 100  continue

 120  if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
      bssmax = bssmax/((tss-bssmax)/(rn-2.0))

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, awt)

      ostat = bssmax
      iseg(1) = tmaxi
      iseg(2) = tmaxj

      return
      end

c     function for calculating the full max wtd t-statistic on permuted data
c     using a new approach to maximizing t-statistic
      double precision function wtmaxp(n,px,wts,sx,cwts,al0)
      integer n,al0
      double precision px(n),wts(n),sx(n),cwts(n)
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, ihi1, jlo1, jhi1, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, sij1, sij2, sijmx0, bijbss, awtmax, psrnov2,
     2     psdiff, psrj, psrn, psrnj, awtlo, awthi, awt1, ssq, tss
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j
      double precision, allocatable :: bssbij(:), bssijmax(:), awt(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:)

c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
      if (n .ge. 50) then
         nb = nint(sqrt(dfloat(n)))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb), ibmin(nb), ibmax(nb))
      allocate(bssbij(nb2), bssijmax(nb2), awt(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2))

c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      psmin0 = 0
      psmax0 = 0
      ipsmin0 = n
      ipsmax0 = n
      ssq = 0
      do 20 j = 1, nb
         sx(ilo) = psum + px(ilo)*wts(ilo)
         ssq = ssq + (px(ilo)**2)*wts(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)*wts(i)
            ssq = ssq + (px(i)**2)*wts(i)
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         ibmin(j) = ipsmin
         ibmax(j) = ipsmax
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     adjust global min, max and locations
         if (psmin .lt. psmin0) then
            psmin0 = psmin
            ipsmin0 = ipsmin
         endif
         if (psmax .gt. psmax0) then
            psmax0 = psmax
            ipsmax0 = ipsmax
         endif
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
 20   continue

c     calculate bss for max s_i - min s_i
      psdiff = psmax0 - psmin0
      psrn = cwts(n)
      tss = ssq - (sx(n)/psrn)**2
      psrj = abs(cwts(ipsmax0) - cwts(ipsmin0))
      psrnj = psrj*(psrn-psrj)
      bssmax = (psdiff**2)/psrnj

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      psrnov2 = psrn/2
      l = 0
      nal0 = n - al0
      do 40 i = 1, nb
         do 30 j = i, nb
c     calculate bsslim
            if (i .eq. 1) then
               ilo = 1
            else
               ilo = bb(i-1) + 1
            endif
            ihi = bb(i)
            if (j .eq. 1) then
               jlo = 1
            else
               jlo = bb(j-1) + 1
            endif
            jhi = bb(j)
c     for wCBS calculated hi and lo arc weights instead of lengths
            awthi = cwts(jhi) - cwts(ilo)
            if (jhi - ilo .gt. nal0) then
               awthi = 0
               do 35 k = 1, al0
                  awthi = max(awthi, cwts(nal0+k) - cwts(k))
 35            continue
            endif
            if (i .eq. j) then
               awtlo = cwts(ilo+al0) - cwts(ilo)
               do 36 k = ilo + 1, ihi - al0
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 36            continue
            else if (i+1 .eq. j) then
               awtlo = cwts(jlo) - cwts(jlo-al0)
               do 37 k = jlo - al0 + 1, ihi
                  awtlo = min(awtlo, cwts(k+al0) - cwts(k))
 37            continue
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            psrnj = min(awtlo*(psrn-awtlo), awthi*(psrn-awthi))
            bsslim = (sijmx0**2)/psrnj
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  awt(l) = abs(cwts(ibmax(j)) - cwts(ibmin(i)))
                  bssbij(l) = (sij1**2)/(awt(l)*(psrn-awt(l)))
               else
                  awt(l) = abs(cwts(ibmin(j)) - cwts(ibmax(i)))
                  bssbij(l) = (sij2**2)/(awt(l)*(psrn-awt(l)))
               endif
            endif
 30      continue
 40   continue
      nb1 = l

c     Now sort the t-statistics by their magnitude
      call qsort4(bssbij, loc, 1, nb1)

c     now go through the blocks in reverse order (largest down)
      do 100 l = nb1, 1, -1
         k = loc(l)
c     need to check a block only if it has potential to increase bss
c     rjlo is the smalllest (j-i) in the block and rjhi is the largest
         bsslim = bssijmax(k)
         if (bssmax .le. bsslim) then
c     bi, bj give the block location
            bi = bloci(k)
            bj = blocj(k)
            awtmax = awt(k)
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1) + 1
            endif
            ihi = bb(bi)
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
            awthi = cwts(jhi) - cwts(ilo)
            if (bi .eq. bj) then
               awtlo = 0
            else
               awtlo = cwts(jlo) - cwts(ihi)
            endif
c
c     if arc wt is larger than half total wt (psrn/2) make is psrn - arc wt
c
            if (awtmax .gt. psrn - awtmax) awtmax = psrn - awtmax
c
c     if awtlo <= psrn/2 start from (ihi, jlo) and go up
c     if awthi >= psrn/2 start from (ilo, jhi) and go down
c
            if (awtlo .le. psrnov2) then
               if (bi .eq.bj) then 
                  ihi1 = ihi - al0
               else
                  ihi1 = ihi
               endif
               do 60 i = ihi1, ilo, -1
                  jlo1 = max(i + al0, jlo)
                  do 55 j = jlo1, jhi
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .le. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) bssmax = bijbss
                     endif
 55               continue
 60            continue
            endif
c
c     make arc wt  psrn - arc wt
c
            awtmax = psrn - awtmax
            if (awthi .ge. psrnov2) then
               do 70 i = ilo, ihi
                  if ((bi .eq. 1) .and. (bj .eq. nb)) then 
                     jhi1 = min(jhi, jhi - al0 + i)
                  else
                     jhi1 = jhi
                  endif
                  do 65 j = jhi1, jlo, -1
                     awt1 = cwts(j) - cwts(i)
                     if (awt1 .ge. awtmax) then
                        bijbss = (sx(j) - sx(i))**2/(awt1*(psrn-awt1))
                        if (bijbss .gt. bssmax) bssmax = bijbss
                     endif
 65               continue
 70            continue
            endif
         endif
 100  continue

      if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
      wtmaxp = bssmax/((tss-bssmax)/(rn-2.0))

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, awt)

      return
      end

c     function for the max (over small arcs) wtd t-statistic on permuted data
c     new code to speed up this part 4/1/2010
      double precision function hwtmaxp(n,k,px,wts,sx,cwts,mncwt,al0)
      integer n,k,al0
      double precision px(n),wts(n),sx(n),cwts(n),mncwt(k)

      integer i, j, nmj, ipj, ipnmj
      double precision rn, rj, rnj, bssmax, bssij, psmin, psmax, psdiff,
     1     bsslim, ssq, tss

c     create blocks of size k (or k+1) to span 1 thru n
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:)
c     variables to work on block specific data
      integer nb, ilo, ihi, l
      double precision psum, psdiffsq

      rn = dfloat(n)
c     number of blocks
      nb = int(rn/dfloat(k))
c     allocate memory
      allocate(bpsmax(nb), bpsmin(nb))
      allocate(bb(nb))
c     block boundaries
      do 110 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 110  continue

c     don't need global min and max
c     find the max, min of partial sums and their locations within blocks
      ilo = 1
      psum = 0
      ssq = 0.0d0
      bssmax = 0.0d0
      rn = cwts(n)
      do 20 j = 1, nb
         sx(ilo) = psum + px(ilo)*wts(ilo)
         ssq = ssq + wts(ilo)*px(ilo)**2
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)*wts(i)
            ssq = ssq + wts(i)*px(i)**2
            if (sx(i) .lt. psmin) then 
               psmin = sx(i)
               ipsmin = i
            endif
            if (sx(i) .gt. psmax) then 
               psmax = sx(i)
               ipsmax = i
            endif
 10      continue
c     store the block min, max and locations
         bpsmin(j) = psmin
         bpsmax(j) = psmax
c     reset ilo to be the block boundary + 1
         psum = sx(bb(j))
         ilo = bb(j) + 1
c     calculate the bss at the block max & min pr
         i = abs(ipsmin - ipsmax)
         if ((i .le. k) .and. (i .ge. al0)) then
            rj = abs(cwts(ipsmax) - cwts(ipsmin))
            rnj = rj*(rn-rj)
            bssij = (bpsmax(j) - bpsmin(j))**2/rnj
            if (bssmax .lt. bssij) bssmax = bssij
         endif
 20   continue
      tss = ssq - (sx(n)/rn)**2

c     check the first block
      ilo = 1
      ihi = bb(1)
      psdiff = bpsmax(1) - bpsmin(1)
      psdiffsq = psdiff**2
      do 40 j = al0,k
         rj = mncwt(j)
         bsslim = psdiffsq/(rj*(rn-rj))
         if (bsslim .lt. bssmax) go to 50
         sxmx = 0.0d0
         do 30 i = ilo,ihi-j
            ipj = i+j
            rj = cwts(ipj) - cwts(i)
            bssij = (sx(ipj) - sx(i))**2/(rj*(rn-rj))
            if (bssij .gt. bssmax) bssmax = bssij
 30      continue
 40   continue

c     now the minor arcs spanning the end (n)
 50   psdiff = max(abs(bpsmax(1)-bpsmin(nb)), abs(bpsmax(nb)-bpsmin(1)))
      psdiffsq = psdiff**2
      do 70 j = al0,k
         rj = mncwt(j)
         bsslim = psdiffsq/(rj*(rn-rj))
         if (bsslim .lt. bssmax) go to 100
         nmj = n-j
         do 60 i = 1,j
            ipnmj = i + nmj
            rj = cwts(ipnmj) - cwts(i)
            bssij = (sx(ipnmj) - sx(i))**2/(rj*(rn-rj))
            if (bssij .gt. bssmax) bssmax = bssij
 60      continue
 70   continue

c     now the other blocks
 100  do 200 l = 2,nb
         ilo = bb(l-1)+1
         ihi = bb(l)
         psdiff = bpsmax(l) - bpsmin(l)
         psdiffsq = psdiff**2
         do 140 j = al0,k
            rj = mncwt(j)
            bsslim = psdiffsq/(rj*(rn-rj))
            if (bsslim .lt. bssmax) go to 150
            sxmx = 0.0d0
            do 130 i = ilo,ihi-j
               ipj = i+j
               rj = cwts(ipj) - cwts(i)
               bssij = (sx(ipj) - sx(i))**2/(rj*(rn-rj))
               if (bssij .gt. bssmax) bssmax = bssij
 130        continue
 140     continue
 150     psdiff = max(abs(bpsmax(l)-bpsmin(l-1)), 
     1        abs(bpsmax(l-1)-bpsmin(l)))
         psdiffsq = psdiff**2
         do 170 j = al0,k
            rj = mncwt(j)
            bsslim = psdiffsq/(rj*(rn-rj))
            if (bsslim .lt. bssmax) go to 200
            do 160 i = ilo-j,ilo-1
               ipj = i+j
               rj = cwts(ipj) - cwts(i)
               bssij = (sx(ipj) - sx(i))**2/(rj*(rn-rj))
               if (bssij .gt. bssmax) bssmax = bssij
 160        continue
 170     continue
 200  continue

c      call dblepr("bss max", 7, bssmax, 1)

      if (tss .le. bssmax+0.0001d0) tss = bssmax + 1.0d0
      hwtmaxp = bssmax/((tss-bssmax)/(dfloat(n)-2.0d0))

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb)

      return
      end

c     the new statistic routine doesn't compute mncwt
      subroutine getmncwt(n, cwts, k, mncwt, delta)
      integer n, k
      double precision cwts(n), mncwt(k), delta

      integer i, j, nmj
      double precision rj, rn

      rn = cwts(n)
      do 30 j = 1,k
         mncwt(j) = cwts(j)
         nmj = n-j
         do 10 i = 1,nmj
            rj = cwts(i+j) - cwts(i)
            mncwt(j) = min(mncwt(j), rj)
 10      continue
         do 20 i = 1, j
            rj = cwts(i+nmj) - cwts(i)
            mncwt(j) = min(mncwt(j), rn-rj)
 20      continue
 30   continue

      j = k+1
      nmj = n-j
      delta = cwts(j)
      do 40 i = 1,nmj
         rj = cwts(i+j) - cwts(i)
         delta = min(delta, rj)
 40   continue
      do 50 i = 1, j
         rj = cwts(i+nmj) - cwts(i)
         delta = min(delta, rn-rj)
 50   continue

      delta = delta/cwts(n)

      return
      end
