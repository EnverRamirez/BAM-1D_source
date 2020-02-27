MODULE ModAddTend

  USE Constants, ONLY :	  &
       r8,i8,r4,i4,  &
       rk                ! R/cp=gasr/cp

  USE Sizes,  ONLY : &
       ibMax,   & ! intent(in)
       jbMax,        & ! intent(in)
       kMax ,   &      ! intent(in)
       sl              ! half sigma levels
  USE FieldsDynamics, ONLY : &
       fguo, &  !variables for Nudging
       fgvo, &  !variables for Nudging
       fgto, &  !variables for Nudging
       fgqo, &  !variables for Nudging
       fgug, &  !variables for Geostrophy
       fgvg, &  !variables for Geostrophy
       fgyum,        & ! intent(in)
       fgyvm,        & ! intent(in)
       fgtdm,        & ! intent(in)
       fgvdlnpm,     & ! intent(in)
       fgum,         & ! intent(in)
       fgvm,         & ! intent(in)
       fgtmpm,       & ! intent(in)
       fgqm,         & ! intent(in)
       fglnpm,       & ! intent(in)
       fgyu,         & ! intent(inout)
       fgyv,         & ! intent(inout)
       fgtd,         & ! intent(inout)
       fgqd,         & ! intent(inout)
       fgvdlnp,      &         ! intent(inout)
       fgu,          &
       fgv,          &
       fgtmp,        &
       fgq  ,        &
       fglnps ,      &
       fgice ,       &
       fgicem ,      &
       fgicet,       &
       fgliq,        &
       fgliqm,       &
       fgliqt,       &
       fgvar,        &
       fgvarm,       &
       fgvart,       & 
       omg,          &  ! vertical velocity @ time T (cb/s)
       fgps             ! surface pressure @ time T (cb)

  USE Options, ONLY :       &
      nClass  ,microphys           , slhum 

!Enver (start)
  USE Utils, ONLY: &
       fcor
!Enver (end)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: AddTend
CONTAINS

  ! addtend: finish tendency computations, adding contributions from
  !          old and current time step.

  SUBROUTINE AddTend(dt, nlnminit,grid)
    REAL(KIND=r8)            , INTENT(IN) :: dt
    LOGICAL(KIND=i8)         , INTENT(IN) :: nlnminit    
    CHARACTER(len=*), INTENT(IN) :: grid     
    INTEGER :: ib, jb, k,kk
    REAL(KIND=r8) :: dt2
!Enver (start)
    REAL(KIND=r8)            :: ug(kMax)   
    REAL(KIND=r8)            :: vg(kMax) 
    REAL(KIND=r8)            :: dtndg
    dtndg = 2.0*3600.0  ! two hours
!Enver (end)

!print*, 'kappa=',rk
!print*, 'ps=',fgps, 'size fgps', size(fgps)
!print*, 'sl=', sl, 'size sl', size(sl)
!print*, 'kMax=',kMax
!print*, 'omega=', omg
!stop


    IF(grid=='1D  ') THEN
      jbMax=1
      ibMax=1
    END IF
    !Enver (start)
     ug(1:kMax) = fgug(1,:,1) ! 0.0_r8
     vg(1:kMax) = fgvg(1,:,1) ! 0.0_r8
    !Enver (end)

    print*, 'fguo,vo,to,qo',fguo(1,1,1),fgvo(1,1,1),fgto(1,1,1),fgqo(1,1,1)
    print*, 'fgu,v,t,q',fgum(1,1,1),fgvm(1,1,1),fgtmpm(1,1,1),fgqm(1,1,1)
    print*, 'fgu,v,t,q',fgu(1,1,1),fgv(1,1,1),fgtmp(1,1,1),fgq(1,1,1)
!   
!  
    print*, 'fguo,vo,to,qo',fguo(1,kMax-4,1),fgvo(1,kMax-4,1),fgto(1,kMax-4,1),fgqo(1,kMax-4,1)
    print*, 'fgu,v,t,q',fgum(1,kMax-4,1),fgvm(1,kMax-4,1),fgtmpm(1,kMax-4,1),fgqm(1,kMax-4,1)
    print*, 'fgu,v,t,q',fgu(1,kMax-4,1),fgv(1,kMax-4,1),fgtmp(1,kMax-4,1),fgq(1,kMax-4,1)
    print*, 'fgu ================'

    dt2 = dt + dt
     DO jb = 1, jbMax
      DO k = 1, kMax
       DO ib = 1, ibMax

!n02
!09102017        fgu  (ib,k,jb)   = fgum(ib,k,jb) * 1.0   !Enver082017 + dt2 * fcor*(fgvm(ib,k,jb)-vg(k)) + dt2 * ( fgyu(ib,k,jb)) + 1.0*dt2/dtndg*(fguo(ib,k,jb)-fgu(ib,k,jb)) 
!09102017        fgv  (ib,k,jb)   = fgvm(ib,k,jb)  * 1.0  !Enver082017 - dt2 * fcor*(fgum(ib,k,jb)-ug(k)) + dt2 * ( fgyv(ib,k,jb)) + 1.0*dt2/dtndg*(fgvo(ib,k,jb)-fgv(ib,k,jb))
!09102017        fgtmp(ib,k,jb)   = fgtmpm(ib,k,jb) * ( 1.0 + dt2 * rk * omg(ib,k,jb) / (10.0*fgps(ib,jb)*sl(k)) ) + 0.0* dt2 * ( fgtd(ib,k,jb)) +0.0* dt2/dtndg*((fgto(ib,k,jb)+273.-300.)-1.0*fgtmp(ib,k,jb))
!testEnver16082017        fgq  (ib,k,jb)   = fgqm(ib,k,jb) !testEnver082017 + dt2 * ( fgqd(ib,k,jb)) + 0.0*dt2/dtndg*(fgqo(ib,k,jb)-fgq(ib,k,jb))
!09102017        fgq  (ib,k,jb)   = fgqm(ib,k,jb) + dt2 * ( fgqd(ib,k,jb) ) + 0.0*dt2/dtndg*(fgqo(ib,k,jb)-fgq(ib,k,jb))

!new
!        fgu  (ib,k,jb)   = fgum(ib,k,jb) * 0.9   + dt2 * fcor*(fgvm(ib,k,jb)-vg(k)) + dt2 * ( fgyu(ib,k,jb)) + 1.0*dt2/dtndg*(fguo(ib,k,jb)-fgu(ib,k,jb)) 
!        fgv  (ib,k,jb)   = fgvm(ib,k,jb)  * 0.9  - dt2 * fcor*(fgum(ib,k,jb)-ug(k)) + dt2 * ( fgyv(ib,k,jb)) + 1.0*dt2/dtndg*(fgvo(ib,k,jb)-fgv(ib,k,jb))
!        fgtmp(ib,k,jb)   = fgtmpm(ib,k,jb) * ( 1.0 + dt2 * rk * omg(ib,k,jb) / (10.0*fgps(ib,jb)*sl(k)) ) +  dt2 * ( fgtd(ib,k,jb)) +0.0* dt2/dtndg*((fgto(ib,k,jb)+273.-300.)-fgtmp(ib,k,jb))
!        fgq  (ib,k,jb)   = fgqm(ib,k,jb)   + dt2 * ( fgqd(ib,k,jb)) + 0.0*dt2/dtndg*(fgqo(ib,k,jb)-fgq(ib,k,jb))

!old
                fgu  (ib,k,jb)   = fgum(ib,k,jb)   !+ dt2 * ( fgyu(ib,k,jb))
                fgv  (ib,k,jb)   = fgvm(ib,k,jb)   !+ dt2 * ( fgyv(ib,k,jb))
                fgtmp(ib,k,jb)   = fgtmpm(ib,k,jb) + dt2 * ( fgtd(ib,k,jb))
                fgq  (ib,k,jb)   = fgqm(ib,k,jb)   + dt2 * ( fgqd(ib,k,jb))



       END DO
      END DO
     END DO
       IF (microphys) THEN
          DO jb = 1, jbMax
             DO k = 1, kMax
                DO ib = 1, ibMax
                   fgice  (ib,k,jb)   = fgicem(ib,k,jb)   + dt2 * ( fgicet(ib,k,jb)) 
                   fgliq  (ib,k,jb)   = fgliqm(ib,k,jb)   + dt2 * ( fgliqt(ib,k,jb)) 
                END DO
             END DO
             DO kk=1,nClass
                DO k = 1, kMax
                   DO ib = 1, ibMax
                      fgvar  (ib,k,jb,kk)   = fgvarm(ib,k,jb,kk)   + dt2 * ( fgvart(ib,k,jb,kk)) 
                   END DO
                END DO
             END DO
          END DO
       ENDIF
       
       DO jb = 1, jbMax
          DO ib = 1, ibMax
             fglnps(ib,jb) = fglnpm(ib,jb) !+ dt2 * ( fgvdlnp(ib,jb) + fgvdlnpm(ib,jb) )
          END DO
       END DO
  END SUBROUTINE AddTend
END MODULE ModAddTend
