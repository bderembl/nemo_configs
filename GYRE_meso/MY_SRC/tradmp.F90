MODULE tradmp
   !!======================================================================
   !!                       ***  MODULE  tradmp  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!======================================================================
   !! History :  OPA  ! 1991-03  (O. Marti, G. Madec)  Original code
   !!                 ! 1992-06  (M. Imbard)  doctor norme
   !!                 ! 1998-07  (M. Imbard, G. Madec) ORCA version
   !!            7.0  ! 2001-02  (M. Imbard)  add distance to coast, Original code
   !!            8.1  ! 2001-02  (G. Madec, E. Durand)  cleaning
   !!  NEMO      1.0  ! 2002-08  (G. Madec, E. Durand)  free form + modules
   !!            3.2  ! 2009-08  (G. Madec, C. Talandier)  DOCTOR norm for namelist parameter
   !!            3.3  ! 2010-06  (C. Ethe, G. Madec) merge TRA-TRC 
   !!            3.4  ! 2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!            3.6  ! 2015-06  (T. Graham)  read restoring coefficient in a file
   !!            3.7  ! 2015-10  (G. Madec)  remove useless trends arrays
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_dmp_alloc : allocate tradmp arrays
   !!   tra_dmp       : update the tracer trend with the internal damping
   !!   tra_dmp_init  : initialization, namlist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean: variables
   USE dom_oce        ! ocean: domain variables
   USE c1d            ! 1D vertical configuration
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE zdf_oce        ! ocean: vertical physics
   USE phycst         ! physical constants
   USE dtatsd         ! data: temperature & salinity
   USE zdfmxl         ! vertical physics: mixed layer depth
   !
   USE in_out_manager ! I/O manager
   USE iom            ! XIOS
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing

   USE lbclnk         ! for filter

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_dmp        ! called by step.F90
   PUBLIC   tra_dmp_init   ! called by nemogcm.F90

   !                                           !!* Namelist namtra_dmp : T & S newtonian damping *
   LOGICAL            , PUBLIC ::   ln_tradmp   !: internal damping flag
   LOGICAL            , PUBLIC ::   ln_filter   !: filter flag
   LOGICAL            , PUBLIC ::   ln_gaussian_flt   !: gaussian filter
   LOGICAL            , PUBLIC ::   ln_shapiro_flt   !: shapiro filter
   INTEGER            , PUBLIC ::   nn_zdmp     !: = 0/1/2 flag for damping in the mixed layer
   CHARACTER(LEN=200) , PUBLIC ::   cn_resto    !: name of netcdf file containing restoration coefficient field
   INTEGER , PUBLIC :: kiter_max
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   resto    !: restoring coeff. on T and S (s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  filter_len    !: filtering length scale
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  filter_iter    !: filtering length scale

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: tradmp.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( resto(jpi,jpj,jpk), STAT= tra_dmp_alloc )
      !
      CALL mpp_sum ( 'tradmp', tra_dmp_alloc )
      IF( tra_dmp_alloc > 0 )   CALL ctl_warn('tra_dmp_alloc: allocation of arrays failed')

      ! BD: TODO: add flag for errors
      ALLOCATE( filter_len(jpi,jpj) )
      ALLOCATE( filter_iter(jpi,jpj) )

      !
   END FUNCTION tra_dmp_alloc


   SUBROUTINE gaussian_filter(ptabin, ptabout) !GIG
     !!----------------------------------------------------------------------
     !!                  ***  routine Gaussian Filter  ***
     !!
     !! ** Purpose :  Multiple application (kiter) of a box filter
     !!               on ptabin to produce ptabout.
     !!
     !! ** Method  :   
     !!
     !! ** Action  :   ptabout filtered output from ptabin
     !!
     !!----------------------------------------------------------------------
     REAL(wp), DIMENSION(:,:), INTENT(IN)  :: ptabin     ! input array
     REAL(wp), DIMENSION(:,:), INTENT(OUT) :: ptabout    ! output array

     ! * Local variable
     INTEGER                               :: ji, jj, jn ! dummy loop index
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zvarout    ! working array
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zvarout2    ! working array
     REAL(wp)                              :: znum       ! numerator
     REAL(wp)                              :: z1_9       ! 1/9
     REAL(wp)                              :: dx         ! dx spacing
     !
     !------------------------------------------------------------------------------
     !
     ! assume uniform grid
     ! dx = e1t(1,1)
     ! sigma^2 = n(w^2 - 1)/12
     ! -> n = 12*sigma^2/(w^2 - 1)
     ! length_scale = 3*sigma

     ! kiter = floor(12.0*(len_cut/3/dx)**2/9.0)

     IF( ln_timing )  CALL timing_start('Gaussian Filter')
     !
     z1_9 = 1._wp/9._wp

     ALLOCATE( zvarout(jpi,jpj) )
     ALLOCATE( zvarout2(jpi,jpj) )

     ptabout(:,:) = ptabin(:,:)
     zvarout(:,:) = ptabout(:,:) 
     zvarout2(:,:) = ptabout(:,:) 

     ! BD: adjust level for tmask ? bug
     DO jn = 1,kiter_max
       DO jj = 2,jpjm1
         DO ji = 2,jpim1
           znum = zvarout(ji,jj  )
           znum = ( zvarout(ji-1,jj-1)*tmask(ji-1,jj-1,1) &
              + zvarout(ji-1,jj  )*tmask(ji-1,jj  ,1) &
              + zvarout(ji-1,jj+1)*tmask(ji-1,jj+1,1) &
              + zvarout(ji,jj-1  )*tmask(ji,jj-1  ,1) &
              + zvarout(ji,jj    )*tmask(ji,jj    ,1) &
              + zvarout(ji,jj+1  )*tmask(ji,jj+1  ,1) &
              + zvarout(ji+1,jj-1)*tmask(ji+1,jj-1,1) &
              + zvarout(ji+1,jj  )*tmask(ji+1,jj  ,1) &
              + zvarout(ji+1,jj+1)*tmask(ji+1,jj+1,1) &
              )/ (tmask(ji-1,jj-1,1) + tmask(ji-1,jj  ,1) &
              + tmask(ji-1,jj+1,1) + tmask(ji,jj-1  ,1) &
              + tmask(ji,jj    ,1) + tmask(ji,jj+1  ,1) &
              + tmask(ji+1,jj-1,1) + tmask(ji+1,jj  ,1) &
              + tmask(ji+1,jj+1,1))
             ptabout(ji,jj)=znum*tmask(ji,jj,1) + ptabin(ji,jj)*(1.-tmask(ji,jj,1))
         ENDDO  ! end loop jj
       ENDDO  ! end loop ji
       !
       CALL lbc_lnk('gaussian_filter', ptabout, 'T', 1.) ! Boundary condition
       zvarout(:,:) = ptabout(:,:)

       DO jj = 2,jpjm1
         DO ji = 2,jpim1
           if (jn == filter_iter(ji,jj)) then
             zvarout2(ji,jj) = zvarout(ji,jj)
           ENDIF
         ENDDO
       ENDDO
     ENDDO  ! end loop jn

     ptabout(:,:) = zvarout2(:,:)
     CALL lbc_lnk('gaussian_filter', ptabout, 'T', 1.) ! Boundary condition

     DEALLOCATE( zvarout )
     DEALLOCATE( zvarout2 )
     IF( ln_timing )  CALL timing_stop('Gaussian Filter')
     !
   END SUBROUTINE gaussian_filter



   SUBROUTINE shapiro_filter(ptabin, ptabout) !GIG
     !!----------------------------------------------------------------------
     !!                  ***  routine Shapiro Filter  ***
     !!
     !! ** Purpose :  Multiple application (kiter) of a box filter
     !!               on ptabin to produce ptabout.
     !!
     !! ** Method  :   
     !!
     !! ** Action  :   ptabout filtered output from ptabin
     !!
     !!----------------------------------------------------------------------
     REAL(wp), DIMENSION(:,:), INTENT(IN)  :: ptabin     ! input array
     REAL(wp), DIMENSION(:,:), INTENT(OUT) :: ptabout    ! output array

     ! * Local variable
     INTEGER                               :: ji, jj, jn ! dummy loop index
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zvarout    ! working array
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zvarout2    ! working array
     REAL(wp)                              :: znum       ! numerator
     REAL(wp)                              :: z1_9       ! 1/9
     REAL(wp)                              :: dx         ! dx spacing
     REAL(wp)                              :: zalphax    ! weight coeficient (x direction)
     REAL(wp)                              :: zalphay    ! weight coeficient (y direction)
     !------------------------------------------------------------------------------
     !
     IF( ln_timing )  CALL timing_start('Shapiro Filter')
     !
     z1_9 = 1._wp/9._wp
     zalphax=1./2.
     zalphay=1./2.

     ALLOCATE( zvarout(jpi,jpj) )
     ALLOCATE( zvarout2(jpi,jpj) )

     ptabout(:,:) = ptabin(:,:)
     zvarout(:,:) = ptabout(:,:) 
     zvarout2(:,:) = ptabout(:,:) 

     ! BD: adjust level for tmask ? bug
     DO jn = 1,kiter_max
       DO jj = 2,jpjm1
         DO ji = 2,jpim1
           znum = zvarout(ji,jj  )
           znum = zvarout(ji,jj)   &
              + 0.25*zalphax*(zvarout(ji-1,jj  )-zvarout(ji,jj))*tmask(ji-1,jj  ,1)  &
              + 0.25*zalphax*(zvarout(ji+1,jj  )-zvarout(ji,jj))*tmask(ji+1,jj  ,1)  &
              + 0.25*zalphay*(zvarout(ji  ,jj-1)-zvarout(ji,jj))*tmask(ji  ,jj-1,1)  &
              + 0.25*zalphay*(zvarout(ji  ,jj+1)-zvarout(ji,jj))*tmask(ji  ,jj+1,1)
           ptabout(ji,jj)=znum*tmask(ji,jj,1) + ptabin(ji,jj)*(1.-tmask(ji,jj,1))
         ENDDO  ! end loop jj
       ENDDO  ! end loop ji
       !
       CALL lbc_lnk('shapiro_filter', ptabout, 'T', 1.) ! Boundary condition
       zvarout(:,:) = ptabout(:,:)

       DO jj = 2,jpjm1
         DO ji = 2,jpim1
           if (jn == filter_iter(ji,jj)) then
             zvarout2(ji,jj) = zvarout(ji,jj)
           ENDIF
         ENDDO
       ENDDO
     ENDDO  ! end loop jn

     ptabout(:,:) = zvarout2(:,:)
     CALL lbc_lnk('shapiro_filter', ptabout, 'T', 1.) ! Boundary condition

     DEALLOCATE( zvarout )
     DEALLOCATE( zvarout2 )
     IF( ln_timing )  CALL timing_stop('Shapiro Filter')
     !
   END SUBROUTINE shapiro_filter



   SUBROUTINE tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp  ***
      !!                  
      !! ** Purpose :   Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - tsa: tracer trends updated with the damping trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts)     ::  zts_dta
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ztrdts
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tsb_flt
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_dmp')
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdts(jpi,jpj,jpk,jpts) ) 
         ztrdts(:,:,:,:) = tsa(:,:,:,:) 
      ENDIF

      ALLOCATE( tsb_flt(jpi,jpj,jpk,jpts) ) 
      tsb_flt(:,:,:,:) = tsb(:,:,:,:) 

      IF( ln_filter )   THEN   
        DO jn = 1, jpts
          DO jk = 1, jpkm1
            IF (ln_gaussian_flt) CALL gaussian_filter(tsb(:,:,jk,jn), tsb_flt(:,:,jk,jn))
            IF (ln_shapiro_flt) CALL shapiro_filter(tsb(:,:,jk,jn), tsb_flt(:,:,jk,jn))
          ENDDO
        ENDDO
      ENDIF
      !                           !==  input T-S data at kt  ==!
      CALL dta_tsd( kt, zts_dta )            ! read and interpolates T-S data at kt
      !
      SELECT CASE ( nn_zdmp )     !==  type of damping  ==!
      !
      CASE( 0 )                        !*  newtonian damping throughout the water column  *!
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     tsa(ji,jj,jk,jn) = tsa(ji,jj,jk,jn) + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jn) - tsb_flt(ji,jj,jk,jn) )
                  END DO
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                       !*  no damping in the turbocline (avt > 5 cm2/s)  *!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( avt(ji,jj,jk) <= avt_c ) THEN
                     tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb_flt(ji,jj,jk,jp_tem) )
                     tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb_flt(ji,jj,jk,jp_sal) )
                  ENDIF
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                       !*  no damping in the mixed layer   *!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( gdept_n(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_tem) - tsb_flt(ji,jj,jk,jp_tem) )
                     tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)   &
                        &                 + resto(ji,jj,jk) * ( zts_dta(ji,jj,jk,jp_sal) - tsb_flt(ji,jj,jk,jp_sal) )
                  ENDIF
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( l_trdtra )   THEN       ! trend diagnostic
         ztrdts(:,:,:,:) = tsa(:,:,:,:) - ztrdts(:,:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_dmp, ztrdts(:,:,:,jp_tem) )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_dmp, ztrdts(:,:,:,jp_sal) )
         DEALLOCATE( ztrdts ) 
      ENDIF

      DEALLOCATE(tsb_flt)

      !                           ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' dmp  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_dmp')
      !
   END SUBROUTINE tra_dmp


   SUBROUTINE tra_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the namtra_dmp namelist and check the parameters
      !!----------------------------------------------------------------------
      INTEGER ::   ji,jj, ioptio, ios, imask   ! local integers 
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  filter_iter_real
      !
      NAMELIST/namtra_dmp/ ln_tradmp, ln_filter, ln_gaussian_flt, ln_shapiro_flt,&
         nn_zdmp, cn_resto
      !!----------------------------------------------------------------------
      !
      ALLOCATE( filter_iter_real(jpi,jpj) )
      ! temporary
      ln_filter = .FALSE.
      ln_gaussian_flt = .FALSE.
      ln_shapiro_flt = .FALSE.

      REWIND( numnam_ref )   ! Namelist namtra_dmp in reference namelist : T & S relaxation
      READ  ( numnam_ref, namtra_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_dmp in reference namelist', lwp )
      !
      REWIND( numnam_cfg )   ! Namelist namtra_dmp in configuration namelist : T & S relaxation
      READ  ( numnam_cfg, namtra_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_dmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_dmp )
      !
      IF(lwp) THEN                  ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dmp_init : T and S newtonian relaxation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_dmp : set relaxation parameters'
         WRITE(numout,*) '      Apply relaxation   or not       ln_tradmp   = ', ln_tradmp
         WRITE(numout,*) '      relax large scale only          ln_filter   = ', ln_filter
         WRITE(numout,*) '         mixed layer damping option      nn_zdmp  = ', nn_zdmp
         WRITE(numout,*) '         Damping file name               cn_resto = ', cn_resto
         WRITE(numout,*) '      Gaussian filter       ln_gaussian   = ', ln_gaussian_flt
         WRITE(numout,*) '      shapiro filter       ln_shapiro   = ', ln_shapiro_flt
         WRITE(numout,*)
      ENDIF
      !
      IF( ln_tradmp ) THEN
         !                          ! Allocate arrays
         IF( tra_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'tra_dmp_init: unable to allocate arrays' )
         !
         SELECT CASE (nn_zdmp)      ! Check values of nn_zdmp
         CASE ( 0 )   ;   IF(lwp) WRITE(numout,*) '   tracer damping as specified by mask'
         CASE ( 1 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixing layer (kz > 5 cm2/s)'
         CASE ( 2 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed  layer'
         CASE DEFAULT
            CALL ctl_stop('tra_dmp_init : wrong value of nn_zdmp')
         END SELECT
         !
         !!TG: Initialisation of dtatsd - Would it be better to have dmpdta routine
         !    so can damp to something other than intitial conditions files?
         !!gm: In principle yes. Nevertheless, we can't anticipate demands that have never been formulated.
         IF( .NOT.ln_tsd_dmp ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout, *)  '   read T-S data not initialized, we force ln_tsd_dmp=T'
            CALL dta_tsd_init( ld_tradmp=ln_tradmp )        ! forces the initialisation of T-S data
         ENDIF
         !                          ! Read in mask from file
         CALL iom_open ( cn_resto, imask)
         CALL iom_get  ( imask, jpdom_autoglo, 'resto', resto )
         CALL iom_close( imask )
         !                          ! Read 
         ioptio = 0
         if (ln_filter) then
           if (ln_gaussian_flt) ioptio = ioptio + 1
           if (ln_shapiro_flt)  ioptio = ioptio + 1
           IF( ioptio /= 1 )    CALL ctl_stop( 'tra_dmp_init: one and only one filter type has to be defined ' )
         endif

         if (ln_filter) then
           CALL iom_open ( cn_resto, imask)
           CALL iom_get  ( imask, jpdom_autoglo, 'filter_len',  filter_len)
           CALL iom_close( imask )
         endif

         ! construct iteration array from length scale
         ! assume uniform grid X/Y
         if (ln_gaussian_flt) then
           DO jj = 1, jpj
             DO ji = 1, jpi
               filter_iter(ji,jj) = floor(12.0*(filter_len(ji,jj)/3/e1t(ji,jj))**2/9.0)
               filter_iter_real(ji,jj) = floor(12.0*(filter_len(ji,jj)/3/e1t(ji,jj))**2/9.0)
             ENDDO
           ENDDO
         else if(ln_shapiro_flt) then
           DO jj = 1, jpj
             DO ji = 1, jpi
               ! BD: BUG: update this for shapiro filter
               filter_iter(ji,jj) = floor(20.0*(filter_len(ji,jj)/3/e1t(ji,jj))**2/9.0)
               filter_iter_real(ji,jj) = floor(20.0*(filter_len(ji,jj)/3/e1t(ji,jj))**2/9.0)
             ENDDO
           ENDDO
         endif

         kiter_max = MAXVAL( filter_iter_real(:,:) )
         CALL mpp_max( 'tra_dmp_init', kiter_max )                 ! max over the global domain

      ENDIF
      DEALLOCATE(filter_iter_real)
      !
   END SUBROUTINE tra_dmp_init

   !!======================================================================
END MODULE tradmp
