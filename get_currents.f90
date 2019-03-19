
program get_currents

  use ModKind, ONLY: Real8_
  use ModTimeConvert
  use ModNumConst
  
  implicit none

  real, parameter :: Re = 6372.0
  
  integer, parameter :: nLons = 48
  integer, parameter :: nLats = 30
  integer, parameter :: nAlts = 20
    
  real :: dlat = 1.0
  real :: dlon = 7.5

  real :: StartLat = 60.0 ! deg
  real :: StartLon = 1.0 ! deg
  real :: StartAlt = 100.0 ! km
  real :: EndAlt   = 150.0 ! km

  real :: NormalPeakAlt = 105.0
  real :: NormalPeakWidth = 10.0
  real :: NormalIntegral
  
  real :: GeoLats(nLats), cosLats(nLats)
  real :: GeoLons(nLons)
  real :: GeoAlts(nAlts)
  real :: Normalization(nAlts)

  real :: dalt

  integer :: iLon, iLat, iAlt, iVar
  integer :: iOutputUnit_ = 31

  real :: B0(nLons, nLats, nAlts,4)
  real :: bDir(nLons, nLats, nAlts,3)
  real :: ApexLat(nLons, nLats, nAlts)
  real :: ApexLon(nLons, nLats, nAlts)
  real :: mlt(nLons, nLats, nAlts)

  real :: Date
  real :: alat, alon, xmag, ymag, zmag, bmag, magpot, lshell

  integer :: iError

  integer      :: iTime(7), jDay
  real(Real8_) :: CurrentTime, StartTime
  real :: SubSolarLatitude, SubSolarLongitude
  real :: MagneticPoleColat, MagneticPoleLon, MagneticPoleStrength
  
  integer, parameter :: iCharLen_ = 100
  integer, parameter :: iCharVarLen_ = 40
  
  character (len=iCharLen_), dimension(10) :: Lines

  ! Here I am setting the file names of the IE files to read:
  character (len=iCharLen_) :: cAMIEFileSouth = "b20090722n.swmf"
  character (len=iCharLen_) :: cAMIEFileNorth = "b20090722s.swmf"

  character (len=iCharLen_) :: outfilePre = "test", outfileNum, outfile

  character (len=iCharVarLen_) :: VarName
  real :: Version = 0.1
  integer, parameter :: nVars = 7
  real :: AllData(nLons, nLats, nAlts, nVars)

  real :: TempPotential(nLons, nLats)
  real :: TempHall(nLons, nLats)
  real :: TempPed(nLons, nLats)

  ! These are the state variables:
  real :: Potential(nLons,nLats,nAlts)
  real :: Efield(nLons,nLats,nAlts,3)
  real :: ExB(nLons,nLats,nAlts,3)
  real :: ExBdir(nLons,nLats,nAlts,3)
  real :: Hall(nLons,nLats,nAlts)
  real :: Ped(nLons,nLats,nAlts)

  real :: jPed(nLons,nLats,nAlts,3)
  real :: jHall(nLons,nLats,nAlts,3)
  real :: jHor(nLons,nLats,nAlts,3)
  real :: jDiv(nLons,nLats,nAlts,3)

  real :: edotb, ExBmag, EFieldMag, div
  real :: eMax, jMax

  ! --------------------------------------------------------------
  ! Here I am setting the number of iterations to take and the dt
  ! between the iterations
  ! --------------------------------------------------------------
  
  integer :: nTimes = 3, iter
  real :: dt = 60.0

  ! Set the initial time (obviously, you have to change the time
  ! if you want to change the time to investigate with new files):

  itime(1) = 2009
  itime(2) = 07
  itime(3) = 22
  itime(4) = 06
  itime(5) = 00
  itime(6) = 00
  itime(7) = 00
  call time_int_to_real(itime, StartTime)
  jDay = n_day_of_year(iTime(1), iTime(2), iTime(3))
  Date = float(iTime(1)) + float(jDay)/365.0
  
  ! Set the delta-altitude
  
  dalt = (EndAlt-StartAlt)/nAlts
  
  ! This is setting the values for the IE reader:
  Lines(1) = "#AMIEFILES"
  Lines(2) = cAMIEFileNorth
  Lines(3) = cAMIEFileSouth
  Lines(4) = ""
  Lines(5) = "#DEBUG"
  Lines(6) = "10"
  Lines(7) = "0"
  Lines(8) = ""
  Lines(9) = "#END"
  Lines(10) = ""

  call EIE_set_inputs(Lines,10)
  call EIE_Initialize(iError)
  
  ! Set the geographic grid to investigate:

  do iLon = 1, nLons
     GeoLons(iLon) = StartLon + dlon*(iLon-1)
  enddo
  do iLat = 1, nLats
     GeoLats(iLat) = StartLat + dlat*(iLat-1)
     cosLats(iLat) = cos(cosLats(iLat)*cDegToRad)
  enddo

  ! This is a total, 100% hack - it just sets the altitude distribution
  ! of the conductances.  This is really not a great thing, and I just made
  ! it up.  We should look up in a book for both a Pedersen and Hall
  ! given a certain energy input and take down the distribution.  But, I
  ! didn't do that.  I made a total hack:

  NormalIntegral = 0.0
  do iAlt = 1, nAlts
     GeoAlts(iAlt) = StartAlt + dalt*(iAlt-1)
     Normalization(iAlt) = dalt*(iAlt-1) * &
          exp(-(GeoAlts(iAlt)-NormalPeakAlt)**2/ &
          NormalPeakWidth**2)
     NormalIntegral = NormalIntegral + &
          Normalization(iAlt) * dAlt * 1000.0
  enddo
  Normalization = Normalization/NormalIntegral

  ! Need a dummy call to initialize the APEX stuff:
  call APEX(DATE, 65.0, 250.0, 100.0, &
       LShell, alat,alon,bmag,xmag,ymag,zmag,MagPot)

  ! Calculate the location of the subsolar point
  call SUBSOLR(iTime(1),jDay,iTime(4),iTime(5),iTime(6), &
       SubsolarLatitude, SubsolarLongitude)

  ! Calculate the location of the dipole pole
  call dypol( MagneticPoleColat, MagneticPoleLon, MagneticPoleStrength)

  ! Fill in magnetic coordinates for each of the geographic grid centers:

  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           alat = 0.0
           alon = 0.0
           xmag = 0.0
           ymag = 0.0
           zmag = 0.0
           MagPot = 0.0
           LShell = 0.0

           call APEX(DATE, &
                GeoLats(iLat), &
                GeoLons(iLon), &
                GeoAlts(iAlt), &
                LShell, alat,alon,bmag,xmag,ymag,zmag,MagPot)

           ! 1 = North
           ! 2 = East
           ! 3 = vertical
           ! 4 = magnitude
           
           B0(iLon,iLat,iAlt,1) = xmag * 1.0e-9
           B0(iLon,iLat,iAlt,2) = ymag * 1.0e-9
           B0(iLon,iLat,iAlt,3) = -zmag * 1.0e-9
           B0(iLon,iLat,iAlt,4) = sqrt( &
                B0(iLon,iLat,iAlt,1) ** 2 + &
                B0(iLon,iLat,iAlt,2) ** 2 + &
                B0(iLon,iLat,iAlt,3) ** 2)

           ! This defines a b-hat (direction of B)
           
           bDir(iLon,iLat,iAlt,1) = &
                B0(iLon,iLat,iAlt,1) / B0(iLon,iLat,iAlt,4)
           bDir(iLon,iLat,iAlt,2) = &
                B0(iLon,iLat,iAlt,2) / B0(iLon,iLat,iAlt,4)
           bDir(iLon,iLat,iAlt,3) = &
                B0(iLon,iLat,iAlt,3) / B0(iLon,iLat,iAlt,4)
           
           ! Store the latitude and londitude

           ApexLat(iLon,iLat,iAlt) = alat
           ApexLon(iLon,iLat,iAlt) = alon

        enddo
     enddo
  enddo

  ! This sets the number of longitudes and latitudes that IE will be expecting:

  call UA_SetnMLTs(nLons)
  call UA_SetnLats(nLats)
  call UA_SetNorth

  ! -----------------------------------------------------------------------
  ! Looping can start here, since this stuff is all dependent on the time:
  ! -----------------------------------------------------------------------

  do iter = 0,nTimes-1

     CurrentTime = StartTime + dt * iter
     call time_real_to_int(CurrentTime, iTime)
     write(*,*) "Working on Time : ", iTime
     jDay = n_day_of_year(iTime(1), iTime(2), iTime(3))
     
     ! Calculate the location of the subsolar point
     call SUBSOLR(iTime(1),jDay,iTime(4),iTime(5),iTime(6), &
          SubsolarLatitude, SubsolarLongitude)
     
     ! Calculate the magnetic local time of each grid point
     do iLon = 1, nLons
        do iLat = 1, nLats
           do iAlt = 1, nAlts
              call magloctm( &
                   ApexLon(iLon,iLat,iAlt), &
                   SubsolarLatitude,   &
                   SubsolarLongitude,  &
                   MagneticPoleColat, &
                   MagneticPoleLon,   &
                   MLT(iLon,iLat,iAlt))
           enddo
        enddo
     enddo

     call get_AMIE_values(CurrentTime)
  
     ! Loop through altitudes, since the magnetic coordinates are dependent on
     ! Altitude:
  
     do iAlt=1,nAlts

        ! Set the grid for IE:
        call UA_SetGrid(MLT(:,:,iAlt), ApexLat(:,:,iAlt), iError)

        ! Get the potential:
        TempPotential = 0.0
        call UA_GetPotential(TempPotential, iError)
        potential(:,:,iAlt) = TempPotential

        ! Get height-integrated Hall Conductance:
        call UA_GetHall(TempHall, iError)
        ! Now spread it out over altitude
        hall(:,:,iAlt) = TempHall * Normalization(iAlt)

        ! Get height-integrated Pedersen Conductance:
        call UA_GetPed(TempPed, iError)
        ! Now spread it out over altitude
        ped(:,:,iAlt) = TempPed * Normalization(iAlt)

     enddo

     eMax = 0.0
     jMax = 0.0
  
     do iLon = 2, nLons-1
        do iLat = 2, nLats-1
           do iAlt = 2, nAlts-1

              ! Longitude Direction
              EField(iLon,iLat,iAlt,1) = - &
                   (Potential(iLon+1,iLat,iAlt) - Potential(iLon-1,iLat,iAlt)) / &
                   (2*dlon*cDegToRad*((Re+GeoAlts(iAlt))*1000.0)*cosLats(iLat))
           
              ! Latitude Direction
              EField(iLon,iLat,iAlt,2) = - &
                   (Potential(iLon,iLat+1,iAlt) - Potential(iLon,iLat-1,iAlt)) / &
                   (2*dlat*cDegToRad*((Re+GeoAlts(iAlt))*1000.0))
           
              ! Altitude Direction
              EField(iLon,iLat,iAlt,3) = - &
                   (Potential(iLon,iLat,iAlt+1) - Potential(iLon,iLat,iAlt-1)) / &
                   (2*dAlt*1000.0)

              ! Remove any E in the B-field direction
           
              edotb = &
                   EField(iLon,iLat,iAlt,1) * bDir(iLon,iLat,iAlt,1) + &
                   EField(iLon,iLat,iAlt,2) * bDir(iLon,iLat,iAlt,2) + &
                   EField(iLon,iLat,iAlt,3) * bDir(iLon,iLat,iAlt,3)           
           
              EField(iLon,iLat,iAlt,1) = EField(iLon,iLat,iAlt,1) - &
                   edotb * bDir(iLon,iLat,iAlt,1)
              EField(iLon,iLat,iAlt,2) = EField(iLon,iLat,iAlt,2) - &
                   edotb * bDir(iLon,iLat,iAlt,2)
              EField(iLon,iLat,iAlt,3) = EField(iLon,iLat,iAlt,3) - &
                   edotb * bDir(iLon,iLat,iAlt,3)

              ExB(iLon,iLat,iAlt,1) = &
                   EField(iLon,iLat,iAlt,2) * B0(iLon,iLat,iAlt,3) - &
                   EField(iLon,iLat,iAlt,3) * B0(iLon,iLat,iAlt,2)

              ExB(iLon,iLat,iAlt,2) = - ( &
                   EField(iLon,iLat,iAlt,1) * B0(iLon,iLat,iAlt,3) - &
                   EField(iLon,iLat,iAlt,3) * B0(iLon,iLat,iAlt,1))

              ExB(iLon,iLat,iAlt,3) = &
                   EField(iLon,iLat,iAlt,1) * B0(iLon,iLat,iAlt,2) - &
                   EField(iLon,iLat,iAlt,2) * B0(iLon,iLat,iAlt,1)

              EFieldMag = sqrt( &
                   EField(iLon,iLat,iAlt,1)**2 + &
                   EField(iLon,iLat,iAlt,2)**2 + &
                   EField(iLon,iLat,iAlt,3)**2)

              eMax = max(EFieldMag, eMax)
           
              ExBmag = sqrt( &
                   ExB(iLon,iLat,iAlt,1)**2 + &
                   ExB(iLon,iLat,iAlt,2)**2 + &
                   ExB(iLon,iLat,iAlt,3)**2)

              ExBdir(iLon,iLat,iAlt,:) = &
                   ExB(iLon,iLat,iAlt,:)/ExBmag

              jPed(iLon,iLat,iAlt,:) = &
                   EField(iLon,iLat,iAlt,:) * &
                   Ped(iLon,iLat,iAlt)
                      
              jHall(iLon,iLat,iAlt,:) = &
                   EFieldMag * &
                   ExBdir(iLon,iLat,iAlt,:) * &
                   Hall(iLon,iLat,iAlt)

              jHor(iLon,iLat,iAlt,:) = &
                   jPed(iLon,iLat,iAlt,:) + &
                   jHall(iLon,iLat,iAlt,:)

              jMax = max(sqrt(sum(jHor(iLon,iLat,iAlt,:)**2)), jMax)
          
           enddo
        enddo
     enddo

     write(*,*) eMax, jMax
  
     do iLon = 3, nLons-2
        do iLat = 3, nLats-2
           do iAlt = 3, nAlts-2

              ! Longitude Direction
              div = &
                   (jHor(iLon+1,iLat,iAlt,1) - jHor(iLon-1,iLat,iAlt,1)) / &
                   (2*dlon*cDegToRad*((Re+GeoAlts(iAlt))*1000.0)*cosLats(iLat))
                
              ! Latitude Direction
              div = div + &
                   (jHor(iLon,iLat+1,iAlt,2) - jHor(iLon,iLat-1,iAlt,2)) / &
                   (2*dlat*cDegToRad*((Re+GeoAlts(iAlt))*1000.0))
           
              ! Altitude Direction
              div = div + &
                   (jHor(iLon,iLat,iAlt+1,3) - jHor(iLon,iLat,iAlt-1,3)) / &
                   (2*dAlt*1000.0)

              jDiv(iLon,iLat,iAlt,:) = div * bDir(iLon,iLat,iAlt,:)
              
           enddo
        enddo
     enddo

     ! Output data to file

     write(outfileNum,'(i4.4)') iter
     outfile = trim(outfilePre)//'_'//trim(outfileNum)//'.bin'
  
     open(iOutputUnit_, file=trim(outfile), status="unknown", form="unformatted")

     write(iOutputUnit_) Version
     write(iOutputUnit_) nLons, nLats, nAlts
     write(iOutputUnit_) nVars

     VarName = "Longitude (deg)"
     write(iOutputUnit_) VarName
     do iLat = 1,nLats
        do iAlt = 1,nAlts
           AllData(:,iLat,iAlt,1) = GeoLons * cDegToRad
        enddo
     enddo
  
     VarName = "Latitude (deg)"
     write(iOutputUnit_) VarName
     do iLon = 1,nLons
        do iAlt = 1,nAlts
           AllData(iLon,:,iAlt,2) = GeoLats * cDegToRad
        enddo
     enddo
  
     VarName = "Altitude (deg)"
     write(iOutputUnit_) VarName
     do iLon = 1,nLons
        do iLat = 1,nLats
           AllData(iLon,iLat,:,3) = GeoAlts * 1000.0
        enddo
     enddo
  
     VarName = "Potential (V)"
     write(iOutputUnit_) VarName
     AllData(:,:,:,4) = Potential

     VarName = "jDiv-East (A/m)"
     write(iOutputUnit_) VarName
     AllData(:,:,:,5) = jDiv(:,:,:,1)

     VarName = "jDiv-North (A/m)"
     write(iOutputUnit_) VarName
     AllData(:,:,:,6) = jDiv(:,:,:,2)

     VarName = "jDiv-Vertical (A/m)"
     write(iOutputUnit_) VarName
     AllData(:,:,:,7) = jDiv(:,:,:,3)

     write(iOutputUnit_) iTime
  
     do iVar = 1, nVars
        write(iOutputUnit_) AllData(:,:,:,iVar)
     enddo

     close(iOutputUnit_)

  enddo
     
end program get_currents
