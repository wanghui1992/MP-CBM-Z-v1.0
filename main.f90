
      program main
      use mod_cbmz
      use gas_data
      implicit none

      integer :: lun_input,lun_output1,nchar
      real :: t1,t2,Elapsed_time
      integer :: comp_stime, comp_etime, clock_rate, clock_max

      !dir$ attributes align:64 :: cbmzobj
      type(cbmztype) cbmzobj
      !dir$ attributes align:64 :: cppb
      real ::  cppb(VLEN,ngas_max)
      !dir$ attributes align:64 :: fcld
      real,dimension(VLEN)   :: fcld          ! the coffi of cloud to photo

      character*40 inputfilename, outputfilename, dumword

      real, dimension(:,:), allocatable :: ggas
      real, dimension(:), allocatable :: grlat,grlon,gte,gRH,gpr_atm,gzalt_m
      integer :: i,i2,ig,ilen,test_n
      logical :: linit
      
      allocate(ggas(ngas_max,npt))
      allocate(grlat(npt),grlon(npt),gte(npt),&
               gRH(npt),gpr_atm(npt),gzalt_m(npt))
 

      write(6,*)'   '
      write(6,*)'******************************************************'
      write(6,*)'             Optimized CBM-Z'
      write(6,*)'   '
      write(6,*)'******************************************************'
      write(6,*)'   '
      write(6,*)'   '


!!----- input ----------------------
      lun_input = 10
      write(6,*)'Input filename: cases.input'
      inputfilename="cases.input"

!!----- output ----------------------
      lun_output1 = 20
      nchar = nbllen(inputfilename) - 6
      outputfilename = inputfilename(1:nchar)//'.output'
      open(lun_output1, file = outputfilename)

!!!!!debugging
      do test_n=1,10,1 !repeat ten times

      open(lun_input, file = inputfilename)
        call ReadInputFile(lun_input,ggas,grlat,grlon,gte,gRH,gpr_atm,gzalt_m)
      close(lun_input)

      call system_clock(comp_stime,clock_rate,clock_max)
      call cpu_time(t1) ! time counting


      linit = .false.
      fcld(:) = 1.0
      emission(:,:) = 0.

gas_chemistry:   do i=1,npt,VLEN! npt,VLEN
      ilen = min(16-i+1, VLEN)
      if (.not.linit) then
          cbmzobj%p_het(:,:) = 0.           !init only once
          cbmzobj%p_com(:,:) = 0.           !init only once
          cbmzobj%rl_com(:,:) = 0.           !init only once
          cbmzobj%r2_com(:,:) = 0.           !init only once
          cbmzobj%p_urb(:,:) = 0.           !init only once
          cbmzobj%rl_het(:,:) = 0.           !init only once
          cbmzobj%rl_urb(:,:) = 0.           !init only once
          cbmzobj%p_bio(:,:) = 0.           !init only once
          cbmzobj%rl_bio(:,:) = 0.           !init only once
          cbmzobj%p_mar(:,:) = 0.           !init only once
          cbmzobj%rl_mar(:,:) = 0.           !init only once
          linit = .true.
      endif
      do ig=1,ngas_max
        do i2 = 1, ilen
          cnn(i2,ig) = ggas(ig,i+i2-1)
        enddo !ig
      enddo !ig

      cbmzobj%pmask(:)  =  .false.             ! enabled flag
      do i2 = 1, ilen
        rlon(i2)   =  grlon(i+i2-1)       !rlat the box lat
        rlat(i2)   =  grlat(i+i2-1)       ! as the above ,but lon
        zalt_m(i2) =  gzalt_m(i+i2-1)     ! the altitude(asl) of box(m)
        cbmzobj%pmask(i2)  =  .true.             ! enabled flag
        cbmzobj%rh(i2)     =  gRH(i+i2-1)        ! the rh
        cbmzobj%te(i2)     =  gte(i+i2-1)        ! the temp
        cbmzobj%pr_atm(i2) =  gpr_atm(i+i2-1)  ! the pressure but as the atm
      enddo !i2

      call cbmz(0,cbmzobj,cppb,fcld)  !! need to modify

  enddo gas_chemistry !i


    call cpu_time(t2)
    call system_clock(comp_etime,clock_rate,clock_max)

!!------------------------------------------------------------------

      !write(6,*)'   '
      !write(6,*)'   '
      !write(6,*)'         End of Simulation'
      !write(6,*)'   '
      !write(6,*)'   '
      !write(6,*)'******************************************************'

      Elapsed_time=t2-t1
      write(6,*) test_n,'Elapsed_time=',Elapsed_time
      print *, test_n,'Elapsed wall time:', real(comp_etime-comp_stime)/real(clock_rate)

!!!!!debugging
      end do
      deallocate(ggas)
      deallocate(grlat,grlon,gte,gRH,gpr_atm,gzalt_m)
      stop
      end            

!!**********************************************************************
!!*********************  Subroutines ***********************************


      subroutine ReadInputFile(lin,ggas,grlat,grlon,gte,gRH,gpr_atm,gzalt_m)
      use gas_data

      real, dimension(ngas_max,npt), intent(out) :: ggas
      real, dimension(npt), intent(out) :: grlat,grlon,gte,gRH,gpr_atm,gzalt_m
      character*40 dword
      integer :: ig,i,j,k, x,y,z
      real :: dummy
      
      write(6,*)'Begin reading inputs.'

      read(lin,*)dword ! PARAMETERS

!!----------- begin time from 12:00 (noon) March 21 [min]
      read(lin,*)dword, tbeg_dd, tbeg_hh, tbeg_mm, tbeg_ss, dword
      read(lin,*)dword, trun_dd, trun_hh, trun_mm, trun_ss, dword
      ! forced to 10 hours
      !trun_hh = 10
      read(lin,*)dword, dt_min,dword ! transport time-step [min]
      read(lin,*)dword, dummy,  dword ! longitude [deg]
      read(lin,*)dword, dummy,  dword ! latitude [deg]
      read(lin,*)dword, dummy,dword ! altitude  above mean sea level [m]
      read(lin,*)dword, dummy,    dword ! relative humidity [%]
      read(lin,*)dword, dummy,    dword ! temperature [K]
      read(lin,*)dword, dummy,dword ! pressure [atm]
      read(lin,*)dword, msolar,dword ! msolar flag
      read(lin,*)dword, mphoto,dword ! mphoto flag
      read(lin,*)dword, iprint,dword ! freq of output

!!----------- read lat,lon,te,RH and gas
      read(lin,*)dword ! GAS

      do i=1,npt !, reduced to 16
      read(lin,*)dword, x ! X index
      read(lin,*)dword, y ! Y index
      read(lin,*)dword, z ! Z index
      !if (i.ne.x .or. j.ne.y .or. k.ne.z) then
      !  write(6,*)'inputs error at: ',i,j,k,x,y,z
      !endif
      !write(6,*)'X:',x,' Y:',y,' Z:',z
      read(lin,*)dword, grlat(i) ! LATITUDE
      !write(6,*)'rlat:',grlat(i)
      read(lin,*)dword, grlon(i) ! LONGITUDE
      !write(6,*)'rlon:',grlon(i)
      read(lin,*)gte(i) ! temperature [K]
      !write(6,*)'te:',gte(i)
      read(lin,*)gRH(i) ! relative humidity [%]
      !write(6,*)'RH:',gRH(i)
      read(lin,*)gpr_atm(i) ! pressure [atm]
      !write(6,*)'pr_atm:',gpr_atm(i)
      read(lin,*)gzalt_m(i) ! pressure [atm]
      !write(6,*)'zalt_m:',gzalt_m(i)
      do ig=1, 67!ngas_max
        read(lin,*)ggas(ig,i)
        !write(6,*)'gas(',ig,'):',ggas(ig,i)
      enddo
      ggas(68:ngas_max,i) = 0.
      enddo ! i

      !extend workload size
!      do i=1,3
!        grlat(npt*i/4+1:npt*(i+1)/4) = grlat(1:npt/4)
!        grlon(npt*i/4+1:npt*(i+1)/4) = grlon(1:npt/4)
!        gte(npt*i/4+1:npt*(i+1)/4) = gte(1:npt/4)
!        gRH(npt*i/4+1:npt*(i+1)/4) = gRH(1:npt/4)
!        gpr_atm(npt*i/4+1:npt*(i+1)/4) = gpr_atm(1:npt/4)
!        gzalt_m(npt*i/4+1:npt*(i+1)/4) = gzalt_m(1:npt/4)
!        ggas(:,npt*i/4+1:npt*(i+1)/4) = ggas(:,1:npt/4)
!      enddo
      write(6,*)'Finished reading all inputs...'
      return
      end

