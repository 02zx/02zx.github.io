program main
 IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 character fname*216, dum*216
 PARAMETER (NT=8000000)
 real(8),dimension(NT)    :: rx,ry,rz,qtet
 integer,dimension(NT,4)     :: nei
 real(8),dimension(NT,4)     :: drn
 

 call getarg(1,fname)
 open(10,file=fname,status='OLD',form='FORMATTED',action='READ')

 nl=0
 do 
   nl=nl+1
   read(10,*,end=900)
   read(10,*) num
   do i=1,num
     read(10,"(A20,3F8.3)") dum, rx(i), ry(i), rz(i)
   enddo
   read(10,*) bx,by,bz

!$OMP PARALLEL PRIVATE(i,j,dr,nj,nk,dx_ij,dy_ij,dz_ij,dr_ij,dx_ik,dy_ik,dz_ik,dr_ik,cosijk,sumc) SHARED(rx,ry,rz,nei, drn,qtet)  
!$OMP DO SCHEDULE(dynamic)
   do i = 1, num
     nei(i,:)=-10
     drn(i,:)=1000
     do j = 1, num
       dr = compute_distance(rx(i),ry(i),rz(i),rx(j),ry(j),rz(j),bx,by,bz)
       if(dr>10) cycle
       if(i==j) cycle
       if(dr<drn(i,1)) then 
        drn(i,4)=drn(i,3)
        drn(i,3)=drn(i,2)
        drn(i,2)=drn(i,1)
        drn(i,1)=dr
        nei(i,4)=nei(i,3)
        nei(i,3)=nei(i,2)
        nei(i,2)=nei(i,1)
        nei(i,1)=j
       elseif(dr<drn(i,2)) then
         drn(i,4)=drn(i,3)
         drn(i,3)=drn(i,2)
         drn(i,2)=dr
         nei(i,4)=nei(i,3)
         nei(i,3)=nei(i,2)
         nei(i,2)=j
       elseif(dr<drn(i,3)) then
         drn(i,4)=drn(i,3)
         drn(i,3)=dr
         nei(i,4)=nei(i,3)
         nei(i,3)=j
       elseif(dr<drn(i,4)) then
        drn(i,4)=dr
        nei(i,4)=j
       endif      
     enddo

     sumc=0.d0
     do nj=1,3
       do nk=nj+1,4
         dx_ij = rx(nei(i,nj)) - rx(i); dx_ij = dx_ij - bx*dnint(dx_ij/bx)
         dy_ij = ry(nei(i,nj)) - ry(i); dy_ij = dy_ij - by*dnint(dy_ij/by)
         dz_ij = rz(nei(i,nj)) - rz(i); dz_ij = dz_ij - bz*dnint(dz_ij/bz)
         dr_ij = dsqrt(dx_ij**2+dy_ij**2+dz_ij**2)

         dx_ik = rx(nei(i,nk)) - rx(i); dx_ik = dx_ik - bx*dnint(dx_ik/bx)
         dy_ik = ry(nei(i,nk)) - ry(i); dy_ik = dy_ik - by*dnint(dy_ik/by)
         dz_ik = rz(nei(i,nk)) - rz(i); dz_ik = dz_ik - bz*dnint(dz_ik/bz)
         dr_ik = dsqrt(dx_ik**2+dy_ik**2+dz_ik**2)

         cosjik=(dx_ij*dx_ik + dy_ij*dy_ik + dz_ij*dz_ik)/dr_ij/dr_ik

         sumc=sumc+(cosjik+1.d0/3.d0)**2
       enddo
     enddo
     qtet(i)=1.d0-3.d0/8.d0*sumc
   enddo
!$OMP END DO
!$OMP END PARALLEL

do n=1,num
     write (6, "(F8.3)") qtet(n)
enddo
  enddo
  


900  continue
 return 
contains

function compute_distance(x1,y1,z1,x2,y2,z2,bx,by,bz) result(dr)
real(8),intent(in) ::x1, y1, z1, x2, y2, z2
real(8) :: distance
       dx = x1 - x2; dx = dx - bx*dnint(dx/bx)
       dy = y1 - y2; dy = dy - by*dnint(dy/by)
       dz = z1 - z2; dz = dz - bz*dnint(dz/bz)
       dr = dsqrt(dx*dx+dy*dy+dz*dz)
end function compute_distance
end
