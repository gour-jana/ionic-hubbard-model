program mapping_levels_clus_sys
implicit none
integer,allocatable::pi(:),pj(:),pci(:),pcj(:)
double precision, allocatable::n_total(:),n_total_c(:)

integer::ls,rf,q,c,p,i,j,pbcx,x,y,k,s,pbcy,yc,xc,l,d,dc
d=4
dc=2
allocate(n_total(d**2),n_total_c(dc**2))
allocate(pi(d**2),pj(d**2),pci(dc**2),pcj(dc**2))
!!***************************************************************************************
!   Levels for the system
!!***************************************************************************************
       do l=1,d**2                     !l -> levels of the lattice points
         do y=1,d
          do x=1,d
           if (d*(y-1)+x.eq.l) then 
            pi(l)=x                   !x co-ordinate of the system
            pj(l)=y                   !y co-ordinate of the system
           endif
          enddo
         enddo
        enddo
!________________________________________________________________________________________

!*****************************************************************************************
!  Levels for the cluster
!*****************************************************************************************
       do c=1,dc**2
         do yc=1,dc
          do xc=1,dc
           if (dc*(yc-1)+xc.eq.c) then
            pci(c)=xc                            ! c for cluster
            pcj(c)=yc
           endif
          enddo
         enddo
        enddo

do rf=1,1!d**2
ls=rf
q=0
i=0
x=pi(rf)               !x=pi(rf) > x levels in the system
y=pj(rf)                !y=pj(rf) > y levels in the system
!write(44,*)rf,pi(rf),pj(rf)
do c=1,dc**2
 xc=pci(c)
 yc=pcj(c)
 !write(45,*)c,pci(c),pcj(c)
q=q+1
if(((d-(x-1)).lt.dc).and.((x-1+xc).gt.d))then
pbcx=-d                                        !pbcx > periodic boundary condition along x
else
pbcx=0
endif
if(((d-(y-1)).lt.dc).and.((y-1+yc).gt.d))then
pbcy=-d**2                                      !pbcy > periodic boundary condition along y
else
pbcy=0
endif
p=ls+(q-1)+pbcx+pbcy                       ! p > variable for mapping levels from system to cluster including PBC
n_total(p)=n_total_c(c)
!print*,c,p,m_c(c),th_c(c),ph_c(c)
if (mod(c,dc).eq.0)then
q=0
i=i+1
ls=rf+d*i
endif
write(32,*)p,c
enddo
enddo

end

