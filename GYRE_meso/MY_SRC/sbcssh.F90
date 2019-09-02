MODULE sbcssh
   !!======================================================================
   !!                       ***  MODULE  sbcssh  ***
   !! Surface module :   ssh restoring
   !!======================================================================
   !! History :  4.0  !   2019-09  B. Deremble
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   sbc_ssh        : read ssh in netcdf files 
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition
   USE phycst          ! physical constants
   !
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE lib_fortran     ! distribued memory computing library
   USE iom             ! IOM library
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssh       ! routine called in istate
   PUBLIC   sbc_ssh_init  ! routine called in istate
   
   !                                !!* namsbc_ssh namelist *
   LOGICAL, PUBLIC ::   ln_ssh_init   !: initial condition
   LOGICAL, PUBLIC ::   ln_ssh_dmp    !: ssh damping?
   
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_ssh   ! structure of input fields (file informations, fields read)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssh_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssh  ***
      !!
      !! ** Purpose :   read ssh field in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_ssh
      !!              - Read ssh fields in netcdf files 
      !! ** action  :   ssh      : ssh at kt
      !!---------------------------------------------------------------------
      INTEGER            ::   ierror  ! local integer 
      INTEGER            ::   ios     ! Local integer output status for namelist read
      !!
      CHARACTER(len=100) ::  cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::  sn_ssh   ! informations about the fields to be read
      LOGICAL            ::  lrxios   ! read restart using XIOS?
      !!
      NAMELIST/namsbc_ssh/ cn_dir, sn_ssh, ln_ssh_init, ln_ssh_dmp
      !!----------------------------------------------------------------------

      !BD: not reading in ref for now because requires a model change but it would be better to do it in the future
      !BD: *** this also requires that all fields are in the cfg namelist ***

!       REWIND( numnam_ref )              ! Namelist namsbc_ssh in reference namelist : File for ssh
!       READ  ( numnam_ref, namsbc_ssh, IOSTAT = ios, ERR = 901)
! 901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_ssh in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsbc_ssh in configuration namelist : File for ssh
      READ  ( numnam_cfg, namsbc_ssh, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_ssh in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsbc_ssh )
      !
      ALLOCATE( sf_ssh(1), STAT=ierror )           !* allocate and fill sf_sst (forcing structure) with sn_sst
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssh: unable to allocate sf_ssh structure' )
      !
      CALL fld_fill( sf_ssh, (/ sn_ssh /), cn_dir, 'sbc_ssh', 'Sea surface height', 'namsbc_ssh' )
                                ALLOCATE( sf_ssh(1)%fnow(jpi,jpj,1)   )
      IF( sn_ssh%ln_tint )   ALLOCATE( sf_ssh(1)%fdta(jpi,jpj,1,2) )
      !
      IF( lwp )THEN                                 !* control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namsbc_ssh : SSH forcing'
         WRITE(numout,*) '      Initialization of ocean SSH  with input data   ln_ssh_init   = ', ln_ssh_init
         WRITE(numout,*) '      Damping of ocean ssh toward input data        ln_ssh_dmp = ', ln_ssh_dmp
      ENDIF
      !
   END SUBROUTINE sbc_ssh_init

   SUBROUTINE sbc_ssh( kt , pssh)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssh  ***
      !!
      !! ** Purpose :   read ssh field in netcdf files.
      !!
      !! ** Method  : - Read ssh field in netcdf files 
      !! ** action  :   ssh      : ssh at kt
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)::   kt   ! ocean time step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pssh   ! ssh
      !
      !!----------------------------------------------------------------------
      CALL fld_read( kt, 1, sf_ssh )      !==   read ssh at kt time step   ==!

      pssh(:,:) = sf_ssh(1)%fnow(:,:,1)     

   END SUBROUTINE sbc_ssh
      
   !!======================================================================
END MODULE sbcssh
