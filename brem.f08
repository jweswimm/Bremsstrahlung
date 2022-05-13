!--------------------Functions----------------------
function sigma(Ee,Evv) result(out)
	real, intent(in) :: Ee, Evv	
	real		 :: out
	if (Ee-Evv.lt.0) then 
		out=0
	else 
	out=((7.9E-26)/(Ee*Evv))*log((((sqrt(Ee-Evv))+sqrt(Ee))**2)/Evv)
	end if

end function
!---------------------------------------------------


program brem
implicit none
real 			      :: sigma,dummy,tot,tot2
real, allocatable :: n(:),z(:), phi(:,:), Ee(:),P(:,:),Ev(:),PI(:)
integer                       :: i,ii,iii,ibin,j,jj,jjj,k

!-----------------------Ev Array-----------------------
!This creates an Ev array from 0 -> 10,000 eV with increments of 10

allocate (Ev(1000))
Ev(1)=10 !start array at 10 eV
do i=2,1000
	Ev(i)=Ev(i-1)+10
end do	




!-----------------------Neutral Densities---------------
!This section reads in the Neutral Densities and creates arrays for them.

!Find out how many altitudes there are.
call execute_command_line('wc -l <neutdens.txt > wc.txt' ) !uses the wordcounter command line function to determine how many rows there are in the file
open(unit=12, file='wc.txt') !open dummy txtfile that has the row count
read(12,*) ii !reads in the total number of rows from the dummy text file
ii=ii-1+100 !accounts for the first row + 100 manufactured altitudes (more on these later)


!Allocate the Neutral Density and Altitude arrays
allocate (n(ii))
allocate (z(ii))



!manufacture 100 altitudes first, since the phi file starts at ~0km and neut dens starts at 200km with increments of 2km
!Initialize altitude and density
z(1)=0.0
n(1)=2.552E-21 !this value was given by Dr. Cravens as a good estimate

do j=2,100 !go until 200km
		z(j)=z(j-1)+2E3 !increment the altitude by 2km
		n(j)=2.552E-21
end do


!Our arrays are complete, now we can
!grab the densities from the file

open(unit=10, file='neutdens.txt') !open neutdens file
read(10,*) !discard first line
do i=101, ii !Start from 202km and go until max altitude
	read(10,*) z(i),n(i)
	z(i)=z(i)*(10**3) ! km -> m (fixes units)
	n(i)=n(i)*(10**6) ! cm^-3 -> m^-3 (fixes units)
end do 
!-----------------------Phi-Up------------------------
!This section reads in the Phi-Up data and creates arrays for them.
open(unit=13, file='phi.txt') !opens the phiup file
open(unit=14, file='wc1.txt') !open dummy txtfile

call execute_command_line('wc -l <phi.txt > wc1.txt' ) !use bash command word count
read(14,*) iii !this is the total number of rows in the Phi file

ibin=260 !this is specific to our file, this needs to be generalized with a counter for any file
!the value of 260 from above is the amount of different incident electron energies there are

!How many Flux values are there per incident electron energy level?
iii=(iii-(2*ibin))/ibin !this is specific to our file, this needs to be generalized for any file


!since phi file starts at -88km and it increments by 2km, we must get rid of the first 44 data points
!this allows us to start the altitude at 0 km. Then our arrays will be matched up 
iii=iii-44


!allocate phi and electron energy arrays
allocate (phi(ibin,ii))
allocate (Ee(ibin))






!Grab Electron Energy and PhiUp(Ee,z)
do j=ibin, 1, -1 !since the phi file starts at the highest incident electron energy
		 !count down
	read (13,*) Ee(j) !read in the incident electron energy (the phi file has been edited to only have Ee in this line)
	read (13,*) !this skips the info line (r/s, phiup/ etc..)



	do i=1, 44 !ignore the first 44 points
	read (13,*)
	end do

	do jj=1, ii !start assigning values to the phi array, ignore the first column since we already have the altitudes
		read (13,*) dummy, phi(j,jj)
		phi(j,jj)=phi(j,jj)*1E4
		
	end do

	do i=1, iii-ii !difference between total number of altitudes and total fluxes per energy
	read (13,*) !skip the rest of the rows (because the neutdens file has less data than the phi file)
	end do

end do





!-------------------Main Phase------------------------
!This section creates productions P(Ev,z) and integrates production over altitude for each energy to obtain intensity 
open (unit=100, file='production.txt') 
open (unit=106, file='intensity.txt') 
!open (unit=111, file='testproduction.txt')

allocate (P(1000,ii))
allocate (PI(1000))

  write(100,*) 'z', (Ev(I), I = 1, 1000) !cravens


!PRODUCTION
!E integral (determines prod by integrated cross section and phi over Electron energy)
do i=1, 1000

	do jj=1, ii !z loop 
	tot=0.0
	
		do j=1,ibin !start integral do loop
		if (j.eq.1) then
		tot=sigma(Ee(j),Ev(i))*phi(j,jj)*Ee(j)
		else
		tot=sigma(Ee(j),Ev(i))*phi(j,jj)*(Ee(j)-Ee(j-1))+tot
		end if
	
		end do
	P(i,jj)=4*4*ATAN(1.d0)*n(jj)*tot !4*ATAN(1.d0)=pi

	!write to the file
!	if (i.eq.2) then !we want just 20eV to write to the output document !cravens
!	write (111,*)P(i,jj),z(jj) !cravens
!	else !cravens
!	end if !cravens


	end do !z loop

end do !Ev loop

    do jj=1,ii !cravens
          write(100,*) z(jj), (P(I,jj), I = 1, 1000)
    End do


!INTENSITY
!z integral (takes prod and integrates over altitude) to obtain intensity
do i=1, 1000
	PI(i)=0.0
	do k=1, ii
	if (k.eq.1) then
		PI(i)=P(i,k)*z(k)
			
	else 
		PI(i)=((P(i,k-1)+P(i,k))*(z(k)-z(k-1))/2.0)+PI(i)
			
	end if
	end do !z integral loop

write (106,*) Ev(i), PI(i) !this is the value of the prod integrated over altitude
end do



!--------------------Make Test Documents-------------------------
!
!
!
!
!!Phi vs. Altitude at given energies (50.5 eV)
!open (unit=101, file='phi_alt50.txt')
!
!
!!do j=ibin, 1, -1
!!	write (101,*) Ee(200), 'eV'
!
!	do i=1,ii
!		write(101,*) phi(200,i),z(i)	
!
!
!	end do
!!end do
!
!
!
!!Phi vs. Altitude at given energies (1050 eV)
!open (unit=99, file='phi_alt1050.txt')
!
!
!!do j=ibin, 1, -1
!!	write (99,*) Ee(100), 'eV'
!
!	do i=1,ii
!*****		write(99,*) n(i), phi(100,i),z(i)	
!
!
!	end do
!!end do
!
!
!
!
!
!
!!Flux vs. Electron Energy at given altitudes (400km)
!open (unit=102, file='phi_ee400.txt')
!
!write (102,*) '#', z(101), 'm'
!do i=ibin-1, 1, -1 !the ibin-1 is to get rid of the outlier
!	write(102,*) phi(i,102), Ee(i)
!end do
!
!!Flux vs. Electron Energy at given altitudes (600km)
!open (unit=103, file='phi_ee600.txt')
!
!write (103,*) '#', z(201), 'm'
!do i=ibin-1, 1, -1 
!	write(103,*) phi(i,201), Ee(i)
!end do
!
!
!
!
!-------------------Graphing------------------------------
!The follow commands execute gnuplot with each .plt
!These files are similar to scripts, but don't need to be executed independently
!!call execute_command_line('gnuplot -p production.plt')
!call execute_command_line('gnuplot -p phi_ee600.plt')
!call execute_command_line('gnuplot -p phi_ee400.plt')
!call execute_command_line('gnuplot -p phi_alt1050.plt')
!call execute_command_line('gnuplot -p phi_alt50.plt')
!!call execute_command_line('gnuplot -p intensity.plt')







!-------------------Test File------------------------------
!This file will create a mass file with each array lined up to ensure they are correctly indexed
open (unit=105, file='test.txt')

write (105,*) 'This is a test file for ',Ee(100),' eV'
write (105,*) '       z        ','        Flux        ','     nn        ','  Production        ' !this is clunky, need fix

do i=1, ii
	write (105,*) z(i), phi(100,i),n(i),p(100,i)
end do

!--------------------Housekeeping-----------------------------


!Close all opened documents
!close(99)
close(103)
!close(102)
!close(101)
close(100)
close(10)
close(13)
close(12)
close(14)
close (105)
close (106)
!close (111)

!deallocate arrays
deallocate (n)
deallocate (phi)
deallocate (Ee)
deallocate (P)
deallocate (PI)


!delete wc.txt and wc1.txt from current directory
call execute_command_line('rm wc.txt')
call execute_command_line('rm wc1.txt')

end program brem
