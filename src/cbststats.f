c     new approach to maximizing t-statistic
c     dynamic memory allocation using allocatable arrays 
      subroutine tmaxo(n,x,tss,sx,iseg,ostat,al0,ibin)
      integer n,iseg(2),al0
      double precision x(n),tss,sx(n),ostat
      logical ibin
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, alenmax, i2j, sxmxi,
     2     alenlo, alenhi, tmaxi, tmaxj, ixlo, ixhi, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1, sij2, sijmx0,
     2     absx, sxmx, bijbss, rnov2, psdiff
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j (and max possible)
      double precision, allocatable :: bssbij(:), bssijmax(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:), alen(:)

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
      allocate(bssbij(nb2), bssijmax(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2), alen(nb2))

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
         sx(ilo) = psum + x(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + x(i)
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
      rj = dfloat(abs(ipsmax0 - ipsmin0))
      rnjov1 = rn/(rj*(rn-rj))
      if (ibin) then
         bssmax = rnjov1*(psdiff-0.5)**2
      else
         bssmax = rnjov1*psdiff**2
      endif
      tmaxi = min(ipsmax0, ipsmin0)
      tmaxj = max(ipsmax0, ipsmin0)

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      rnov2 = rn/2
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
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (i .eq. j) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            rjlo = dfloat(alenlo)
            rnjov1 = rn/min(rjlo*(rn-rjlo), rjhi*(rn-rjhi))
            if (ibin) then
               bsslim = rnjov1*(sijmx0-0.5)**2
            else
               bsslim = rnjov1*(sijmx0**2)
            endif
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  alen(l) = abs(ibmax(j) - ibmin(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij1-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij1**2)
                  endif
               else
                  alen(l) = abs(ibmin(j) - ibmax(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij2-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij2**2)
                  endif
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
c     max arc length of interest in block
            alenmax = alen(k)
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
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (bi .eq. bj) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
            rjlo = dfloat(alenlo)
c
c     if arc length is larger than n/2 make is n - arc length
c
            if (alenmax .gt. n - alenmax) alenmax = n - alenmax
c
c     if alenlo <= n/2 start from (ihi, jlo) and go up
c     if alenhi >= n/2 start from (ilo, jhi) and go down
c
            if ((rjlo .le. rnov2) .and. (alenlo .le. alenmax)) then
               do 60 i2j = alenlo, alenmax
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 55 i = ilo + ixlo, ihi - ixhi
                     j = i+i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) then
                        sxmx = absx
                        sxmxi = i
                     endif
 55               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) then
                     bssmax = bijbss
                     tmaxi = sxmxi
                     tmaxj = sxmxi + i2j
                  endif
 60            continue
            endif
c
c     make arclength n - arc length
c
            alenmax = n - alenmax
            if ((rjhi .ge. rnov2) .and. (alenhi .ge. alenmax)) then
               do 70 i2j = alenhi, alenmax, -1
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 65 i = ilo + ixlo, ihi - ixhi
                     j = i + i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) then
                        sxmx = absx
                        sxmxi = i
                     endif
 65               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) then
                     bssmax = bijbss
                     tmaxi = sxmxi
                     tmaxj = sxmxi + i2j
                  endif
 70            continue
            endif
         endif
 100  continue

      if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         bssmax = bssmax/(tss/rn)
      else
         if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
         bssmax = bssmax/((tss-bssmax)/(rn-2.0))
      endif

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, alen)

      ostat = bssmax
      iseg(1) = tmaxi
      iseg(2) = tmaxj

      return
      end

c     function for calculating the full max t-statistic on permuted data
c     new approach to maximizing t-statistic using allocatable arrays 
      double precision function tmaxp(n,tss,px,sx,al0,ibin)
      integer n,al0
      double precision tss,px(n),sx(n)
      logical ibin
c     
c     look at the partial sums in blocks of size sqrt(n)
c     
      integer ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
     1     nb2, bi, bj, ilo, ihi, jlo, jhi, alenmax, i2j, alenlo,
     2     alenhi, ixlo, ixhi, nal0
      double precision psum, psmin, psmax, psmin0, psmax0, bssmax,
     1     bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1, sij2, sijmx0,
     2     absx, sxmx, bijbss, rnov2, psdiff
c     
c     use local arrays for working within blocks
c     block partial sum max and min 
      double precision, allocatable :: bpsmax(:), bpsmin(:)
c     location of the max and min
      integer, allocatable :: bb(:), ibmin(:), ibmax(:)

c     t statistic corresponding to max for block i,j (and max possible)
      double precision, allocatable :: bssbij(:), bssijmax(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:), alen(:)

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
      allocate(bssbij(nb2), bssijmax(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2), alen(nb2))

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
         sx(ilo) = psum + px(ilo)
         psmin = sx(ilo)
         ipsmin = ilo
         psmax = sx(ilo)
         ipsmax = ilo
         do 10 i = ilo+1, bb(j)
            sx(i) = sx(i-1) + px(i)
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
      rj = dfloat(abs(ipsmax0 - ipsmin0))
      rnjov1 = rn/(rj*(rn-rj))
      if (ibin) then
         bssmax = rnjov1*(psdiff-0.5)**2
      else
         bssmax = rnjov1*psdiff**2
      endif

c     for a pair of blocks (i,j) calculate the max absolute t-statistic
c     at the (min_i, max_j) and (max_i, min_j) locations 
c     for other indices the t-statistic can be bounded using this
c
c     if a block doesn't have the potential to exceed bssmax ignore it
c     calculate the bsslim for each block and include ones >= bssmax

      rnov2 = rn/2
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
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (i .eq. j) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
c     max S_k over block j - min S_k over block i
            sij1 = abs(bpsmax(j) - bpsmin(i))
c     max S_k over block i - min S_k over block j
            sij2 = abs(bpsmax(i) - bpsmin(j))
c     if i = j then sij1 and sij2 are the same
            sijmx0 = max(sij1, sij2)
            rjlo = dfloat(alenlo)
            rnjov1 = rn/min(rjlo*(rn-rjlo), rjhi*(rn-rjhi))
            if (ibin) then
               bsslim = rnjov1*(sijmx0-0.5)**2
            else
               bsslim = rnjov1*(sijmx0**2)
            endif
c     if its as large as bssmax add block
            if (bssmax .le. bsslim) then
               l = l+1
               loc(l) = l
               bloci(l) = i
               blocj(l) = j
               bssijmax(l) = bsslim
c     max sij in the (i,j) block, t-statistic etc
               if (sij1 .gt. sij2) then
                  alen(l) = abs(ibmax(j) - ibmin(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij1-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij1**2)
                  endif
               else
                  alen(l) = abs(ibmin(j) - ibmax(i))
                  rj = dfloat(alen(l))
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bssbij(l) = rnjov1*(sij2-0.5)**2
                  else
                     bssbij(l) = rnjov1*(sij2**2)
                  endif
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
c     max arc length of interest in block
            alenmax = alen(k)
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
            alenhi = jhi - ilo
            if (alenhi .gt. nal0) alenhi = nal0
            rjhi = dfloat(alenhi)
            if (bi .eq. bj) then
               alenlo = 1
            else
               alenlo = jlo - ihi
            endif
            if (alenlo .lt. al0) alenlo = al0
            rjlo = dfloat(alenlo)
c
c     if arc length is larger than n/2 make is n - arc length
c
            if (alenmax .gt. n - alenmax) alenmax = n - alenmax
c
c     if alenlo <= n/2 start from (ihi, jlo) and go up
c     if alenhi >= n/2 start from (ilo, jhi) and go down
c
            if ((rjlo .le. rnov2) .and. (alenlo .le. alenmax)) then
               do 60 i2j = alenlo, alenmax
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 55 i = ilo + ixlo, ihi - ixhi
                     j = i+i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) sxmx = absx
 55               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) bssmax = bijbss
 60            continue
            endif
c
c     make arclength n - arc length
c
            alenmax = n - alenmax
            if ((rjhi .ge. rnov2) .and. (alenhi .ge. alenmax)) then
               do 70 i2j = alenhi, alenmax, -1
c     excess calcultaions to set range of i
                  ixlo = max(0, jlo - ilo - i2j)
                  ixhi = max(0, ihi + i2j - jhi)
                  sxmx = 0
                  do 65 i = ilo + ixlo, ihi - ixhi
                     j = i + i2j
                     absx = abs(sx(j) - sx(i))
                     if (sxmx .lt. absx) sxmx = absx
 65               continue
                  rj = dfloat(i2j)
                  rnjov1 = rn/(rj*(rn-rj))
                  if (ibin) then
                     bijbss = rnjov1*(sxmx-0.5)**2
                  else
                     bijbss = rnjov1*(sxmx**2)
                  endif
                  if (bijbss .gt. bssmax) bssmax = bijbss
 70            continue
            endif
         endif
 100  continue

      if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         tmaxp = bssmax/(tss/rn)
      else
         if (tss.le.bssmax+0.0001) tss = bssmax + 1.0
         tmaxp = bssmax/((tss-bssmax)/(rn-2.0))
      endif

c     deallocate memory
      deallocate(bpsmax, bpsmin, bb, ibmin, ibmax)
      deallocate(bssbij, bssijmax, bloci, blocj, loc, alen)

      return
      end

c     function for the max (over small arcs) t-statistic on permuted data
      double precision function htmaxp(n,k,tss,px,sx,al0,ibin)
      integer n,k,al0
      double precision tss,px(n),sx(n)
      logical ibin

      integer i, j, nmj, nmji
      double precision rn, rj, absx, sxmx, bssmx, psmin, psmax, psdiff,
     1     bsslim, rnjov1

      rn = dfloat(n)
      sx(1) = px(1)
      psmin = sx(1)
      psmax = sx(1)
      do 20 i = 2,n
         sx(i) = sx(i-1) + px(i)
         psmin = min(psmin, sx(i))
         psmax = max(psmax, sx(i))
c         if (sx(i) .lt. psmin) psmin = sx(i)
c         if (sx(i) .gt. psmax) psmax = sx(i)
 20   continue
      psdiff = psmax - psmin

      htmaxp = 0.0d0
      do 50 j = al0,k
         rj = dfloat(j)
         rnjov1 = rn/(rj*(rn-rj))
         if (ibin) then
            bsslim = rnjov1*(psdiff-0.5)**2
         else
            bsslim = rnjov1*psdiff**2
         endif
         if (bsslim .lt. htmaxp) go to 60
         sxmx = 0.0d0
         do 30 i = 1,n-j
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 30      continue
         nmj = n - j
         do 40 i = 1, j
            nmji = nmj + i
            absx = abs(sx(nmji) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 40      continue
         if (ibin) then
            bssmx = rnjov1*(abs(sxmx)-0.5)**2
         else
            bssmx = rnjov1*sxmx**2
         endif
         if (htmaxp.lt.bssmx) htmaxp = bssmx
 50   continue
 60   if (ibin) then
         if (tss .le. 0.0001d0) tss = 1.0d0
         htmaxp = htmaxp/(tss/rn)
      else
         if (tss .le. htmaxp+0.0001d0) tss = htmaxp + 1.0d0
         htmaxp = htmaxp/((tss-htmaxp)/(rn-2.0d0))
      endif

      return
      end
