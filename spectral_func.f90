	module array
	 complex*16,allocatable::H(:,:),work1(:)
	 integer,allocatable::pi(:),pj(:)
	 double precision, allocatable::m_s(:),th_s(:),ph_s(:),rwork1(:),evl_s(:),op_fl(:,:)
	 double precision, allocatable::n_total(:),ion_dis(:),spec(:,:)
	end module array
	
	module global
	 integer::d,MCSW,temp_max,intrvl,dc,seed
	 double precision::T,strnth
	 double precision ::t1,t2,U1,filling,gama,gama_m
	end module global

    program spectral_fun
     use array
     use global
     implicit none
     integer ::l,x,y,temp,count1,config,mc_count,q,i,p,ink,iw
     double precision ::ds,inp,get_mu_s,w,kx,ky,pp,omg,mu_sys
     character(len=8 )::grnsfn
	 open (5,file='input.dat',status='unknown')
	 
     do i=1,13
         read (5,*)inp
		 if (i.eq.1)d=int(inp)
		 if (i.eq.2)dc=int(inp)
	     if (i.eq.3)t1=dble(inp)
		 if (i.eq.4)t2=dble(inp)
		 if (i.eq.5)temp_max=int(inp)
		 if (i.eq.6)mcsw=int(inp)
		 if (i.eq.7)intrvl=int(inp)
		 if (i.eq.8)U1=dble(inp)
		 if (i.eq.9)filling=dble(inp)
		 if (i.eq.10)gama=dble(inp)
		 if (i.eq.11)gama_m=dble(inp)
		 if (i.eq.12)seed=int(inp)
		 if (i.eq.13)strnth=dble(inp)
	 end do      
    
    
		print*,"system size=",d
		print*,"cluster size=",dc
		print*,"nn hopping parameter=",t1
		print*,"nnn hopping parameter=",t2
		print*,"maximum no of temperature points=",temp_max
		print*,"total no of system sweeps MCSW=",MCSW
		print*,"interval in system sweeps steps to cal. observables=",intrvl
		print*,"the interaction value U1=",U1
		print*,"filling in the system=",filling
		print*,"the broadening of the lorenzian for cal. DOS=",gama
		print*,"the broadening of the lorentzian for cal. of distributin p(m) of m=",gama_m
		print*,"strnth of the ionic disorder=",strnth       


	      print*,"print gour"

        config=(MCSW/(2*intrvl))+1
        
		allocate(m_s(d**2),th_s(d**2),ph_s(d**2),n_total(d**2),ion_dis(d**2))
		allocate(pi(d**2),Pj(d**2),op_fl(temp_max*config*d**2,5))
		allocate(evl_s(2*d**2),work1(2*(2*d**2)-1),rwork1(3*(2*d**2)-2),H(2*d**2,2*d**2))
		allocate(spec(3*d+3,1000))
		12 format('fort.',I3, '')
		
		ds=filling*d**2 
		pp=acos(-1.0d0)
		 !__________________________________________________________________
     !   Levels for the system

       do l=1,d**2
         do y=1,d
          do x=1,d
           if (d*(y-1)+x.eq.l) then 
            pi(l)=x
            pj(l)=y
           endif
          enddo
         enddo
        enddo
     !__________________________________________________________________
     
         l=1
		do y=1,d
		 do x=1,d
		  ion_dis(l)=strnth*(-1)**(x+y)  !ion_dis(l) -> site dependent ionic disorder with strength "strnth"
          !write(31,*)l,ion_dis(l)
          l=l+1
         enddo
        enddo
  !_____________________________________________________________________
	 
	    do temp=1,temp_max
         write(grnsfn,12) 300+temp
          open(unit=7,file=grnsfn,status='unknown')
           do q=1,config*d**2
            read(7,*)op_fl(q+((temp-1)*config*d**2),1:5)
           enddo
         close(7)
        enddo
	 
	       T=0.1750d0
        do temp=1,temp_max
	      if(temp.le.3)T=T-0.025
	      if((temp.gt.3).and.(temp.le.12))T=T-0.01
		  if(temp.eq.13)T=0.005
          if(temp.eq.14)T=0.001
          mc_count=0
          spec=0.0d0
            
          
	     do count1=1,1!config
	        mc_count=mc_count+1
	        1919 format (1x,i4,4f16.8)
	      do l=1,d**2
            p=(temp-1)*config*d**2+(count1-1)*d**2+l
            m_s(l)=op_fl(p,2)
            th_s(l)=op_fl(p,3)
            ph_s(l)=op_fl(p,4)
            n_total(l)=op_fl(p,5)
          enddo !  end of the site loop 
          
            call system_mat_gen
             mu_sys=get_mu_s(dble(ds)) 
            print*,mu_sys
            call spectrl_func
            
           
             print*,T
          
          
         enddo  !system sweep loop
         
         
            do ink=1,3*d+3
        
               if(ink.le.(d+1))kx=(pp/dble(d))*(ink-1);ky=0.0d0
               if((ink.gt.(d+1)).and.(ink.le.(2*d+2)))kx=0.0d0;ky=(pp/dble(d))*(ink-(d+2))
               if((ink.gt.(2*d+2)).and.(ink.le.(3*d+3)))then
                kx=pp-(pp/dble(d))*(ink-(2*d+3))
                ky=pp-(pp/dble(d))*(ink-(2*d+3))
         
               endif
                w=-30.0d0
                omg=0.060d0
             do iw=1,100
               w=w+omg
 write(6000+temp,*)ink,kx,ky,(w-mu_sys),spec(ink,iw)/dble(mc_count)                                                                                                                                                                                                                                                                                                                               
           
             enddo     ! iw loop
            enddo  !ink loop  
    
         
        enddo   !temperature loop
         
    
    end 
    
    
!********************************************************************************************
!subroutine for matrix formation and diagonalization for the system
!********************************************************************************************
	subroutine system_mat_gen
 !	use input
  	use array
  	use global
 	implicit none
 	integer :: l,ii,id,ji,jd,k,a,b,i,j,info
 	double precision r
	H=cmplx(0.0d0,0.0d0)

  	do l=1,d**2
         ii=1
         id=-1
         ji=1
         jd=-1
         i=pi(l)
         j=pj(l)
 !        write(47,*)l,pi(l),pj(l)
         if (i.eq.1) id=-1+d
         if (i.eq.d) ii=1-d
         if (j.eq.1) jd=-1+d
         if (j.eq.d) ji=1-d
         do k=1,d**2
          if(l.eq.k)then
           a=2*k-1
           b=2*k-1
           H(a,b)=(-abs(m_s(k))*cos(th_s(k))*(U1/2.0d0))+((n_total(k))*U1/2.0d0)+ ion_dis(k)
           H(a+1,b+1)=(abs(m_s(k))*cos(th_s(k))*(U1/2.0d0))+((n_total(k))*U1/2.0d0)+ ion_dis(k)
           H(a,b+1)=-abs(m_s(k))*sin(th_s(k))*cmplx(cos(ph_s(k)),-sin(ph_s(k)))*(U1/2.0d0)
           H(a+1,b)=conjg(H(a,b+1))
          endif
          if (((pi(k).eq.(i+ii)) .and. (pj(k).eq.j))&
          &.or. ((pi(k).eq.i) .and. (pj(k).eq.(j+ji)))) then
           a=2*l-1
           b=2*k-1
           H(a,b)=t1
           H(a+1,b+1)=t1
          endif
          if (((pi(k).eq.(i+id)) .and. (pj(k).eq.j))&
          &.or. ((pi(k).eq.i) .and. (pj(k).eq.(j+jd)))) then
           a=2*l-1
           b=2*k-1
           H(a,b)=t1
           H(a+1,b+1)=t1
          endif
          if (((pi(k).eq.(i+ii)).and.(pj(k).eq.(j+ji)))&
          &.or. ((pi(k).eq.(i+id)).and.(pj(k).eq.(j+jd)))) then
           a=2*l-1
           b=2*k-1
           H(a,b)=t2
           H(a+1,b+1)=t2
          endif
         if (((pi(k).eq.(i+id)).and.(pj(k).eq.(j+ji)))&
          &.or. ((pi(k).eq.(i+ii)).and.(pj(k).eq.(j+jd)))) then
           a=2*l-1
           b=2*k-1
           H(a,b)=t2
           H(a+1,b+1)=t2
          endif
         enddo
        enddo
        
       do l=1,d**2
   !     write(3500,*)m_s(l)
       enddo
       
     call zheev ('V','U',2*d**2,H,2*d**2,evl_s,work1,2*(2*d**2)-1,rwork1,info)
     
     print*,info
    	return
  	end


!****************************************************************************
!subroutine for calculating spectral function A(K,W)
!****************************************************************************

	subroutine spectrl_func
     use array
     use global
     implicit none
     integer::ink,iw,j,l,a
     double precision::pp,kx,ky,w,omg,arg,prd,lrn,sp
     pp=acos(-1.0d0)
     
     
     
        sp=0.0d0
 
        do ink=1,3*d+3
        
         if(ink.le.(d+1))kx=(pp/dble(d))*(ink-1);ky=0.0d0
         if((ink.gt.(d+1)).and.(ink.le.(2*d+2)))kx=0.0d0;ky=(pp/dble(d))*(ink-(d+2))
         if((ink.gt.(2*d+2)).and.(ink.le.(3*d+3)))then
          kx=pp-(pp/dble(d))*(ink-(2*d+3))
          ky=pp-(pp/dble(d))*(ink-(2*d+3))
         
         endif
        
        
           w=-30.0d0
           omg=0.060d0
         do iw=1,100
           w=w+omg
           
          do j=1,d**2
           do l=1,d**2
             arg=(kx*(pi(j)-pi(l))+ky*(pj(j)-pj(l)))
             
             do a=1,2*d**2                                !a--> eigen vector variable
              lrn=((gama/pp)/((w-evl_s(a))**2+(gama**2)))
              prd=lrn*cmplx(cos(arg),sin(arg))
              sp=sp+(H(2*l-1,a)*conjg(H(2*j-1,a))+H(2*l,a)*conjg(H(2*j,a)))*prd
              
             enddo   ! a loop 
              
                
              
            enddo   ! l loop
          
          enddo    ! j loop
         spec(ink,iw)=spec(ink,iw)+sp*(1.0d0/(pp*dble(d**4)))
 
 
 
         enddo     ! iw loop
          
          
          
        enddo  !ink loop  
 
     
     
     
     
     return
	end






!****************************************************************************
!calculation of chemical potential for the system
!***************************************************************************
double precision function get_mu_s(fill)
    use array
!    use input
    use global
	implicit none
	double precision f0, f, fL2, fR, mR, mL, rtmp,m_d
	integer i
	double precision fill
	mR = maxval(evl_s)       !right-side chemical potential
	fr=0.0d0
	do i=1,2*d**2
	fr=fr+(1.0d0/(exp((evl_s(i)-mR)/T)+1.0d0))
	end do
	 mL = minval(evl_s)       !left-side chemical potential
	 fL2=0.0d0
 	do i=1,2*d**2
 	fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/T)+1.0d0))
	end do
	m_d = 0.5d0*(mL+mR)    !middle chemical potential
	f=0.0d0
	do i=1,2*d**2
	f=f+(1.0d0/(exp((evl_s(i)-m_d)/T)+1.0d0))
	end do
	!print*,f,fill
	do while(abs(f-fill).ge.1e-8)
	 m_d = 0.5d0*(mL+mR)
	f=0.0d0
	do i=1,2*d**2
	f=f+(1.0d0/(exp((evl_s(i)-m_d)/T)+1.0d0))
	end do
	 if(f.gt.fill)then
	  !if middle filling is above target, make it the new right bound.
	  mR = m_d
	  fR = f
	 elseif(f.lt.fill)then
	  !if middle filling is below target, make it the new left bound.
	  mL = m_d
	  fR = f
	 endif
	enddo
	!Return the middle value
	get_mu_s = m_d
	return
	end function get_mu_s
























