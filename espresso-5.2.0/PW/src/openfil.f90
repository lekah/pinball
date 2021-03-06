!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE openfil()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens some files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  ! ... All units are set in Modules/io_files.f90
  !
  USE kinds,            ONLY : DP
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level
  USE io_global,        ONLY : stdout,&
                               ionode !aris 
  USE basis,            ONLY : natomwfc, starting_wfc
  USE wvfct,            ONLY : nbnd, npwx
  USE fixed_occ,        ONLY : one_atom_occupations
  USE ldaU,             ONLY : lda_plus_U, U_projection, nwfcU
  USE io_files,         ONLY : prefix, iunpun, iunsat, iunigk,  &
                               iunhub, nwordwfcU, nwordwfc, nwordatwfc,&
                               iunefield, iunefieldm, iunefieldp, seqopn,&
                               iunevp,iunpos,iunvel,tmp_dir, &               !aris
                               iunfor ! LEONID: File for forces
  USE noncollin_module, ONLY : npol
  USE bp,               ONLY : lelfield
  USE wannier_new,      ONLY : use_wannier
  !
  IMPLICIT NONE
  !
  LOGICAL            :: exst
  !
  ! ... Files needed for LDA+U
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunhub as above, only wfcs with a U correction
  !
  ! ... nwordwfc is the record length (IN COMPLEX WORDS)
  ! ... for the direct-access file containing wavefunctions
  ! ... nwordatwfc as above (IN REAL WORDS) for atomic wavefunctions
  !
  if (ionode) then !aris
      open(iunpos,file=trim(tmp_dir)//'../verlet.pos',position='append') !aris
      open(iunevp,file=trim(tmp_dir)//'../verlet.evp',position='append') !aris
      open(iunfor,file=trim(tmp_dir)//'../verlet.for',position='append') ! LEONID
      open(iunvel,file=trim(tmp_dir)//'../verlet.vel',position='append') !aris
      ! write(iunevp,'("NSTEP , TIME(ps), EKIN , ETOT , EKIN+ETOT , TEMP")') !aris
  end if !aris

  nwordwfc  = nbnd*npwx*npol
  nwordatwfc= npwx*natomwfc*npol
  nwordwfcU = npwx*nwfcU*npol
  !
  IF ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) &
     CALL open_buffer ( iunhub, 'hub',    nwordwfcU, io_level, exst )
  IF ( use_wannier .OR. one_atom_occupations ) &
     CALL open_buffer ( iunsat, 'satwfc', nwordatwfc, io_level, exst )
  !
  ! ... iunigk contains the number of PW and the indices igk
  !
  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )
  !
  ! ... open units for electric field calculations
  !
  IF ( lelfield ) THEN
      CALL open_buffer( iunefield , 'ewfc' , nwordwfc, io_level, exst )
      CALL open_buffer( iunefieldm, 'ewfcm', nwordwfc, io_level, exst )
      CALL open_buffer( iunefieldp, 'ewfcp', nwordwfc, io_level, exst )
  END IF
  !
  RETURN
  !
END SUBROUTINE openfil
