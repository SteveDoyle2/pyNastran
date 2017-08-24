@ECHO off



ECHO.
ECHO %INDENT_ECHO_05%   **************************************************************
ECHO %INDENT_ECHO_05%   *** Chargement Variable Environnement LOCAL SATELLITE V02  ***
ECHO %INDENT_ECHO_05%   **************************************************************
ECHO.

SET Satellite_V02_CAS_LOCAL=1

IF "%Satellite_V02_CAS_LOCAL%" == "0" (
  REM CAS RESEAU
  ECHO.
  ECHO VARIABLE ENVIRONEMENT : CAS RESEAU
  SET Satellite_V02_base=\\castore03\save_tsr\ACA\Satellite_V02
)


IF "%Satellite_V02_CAS_LOCAL%" == "1" (
  REM CAS LOCAL
  ECHO.
  ECHO VARIABLE ENVIRONEMENT : CAS LOCAL
  SET Satellite_V02_base=M:\ACA\Satellite_V02
)
  
  
ECHO %INDENT_ECHO_05%      Satellite_V02_base = %Satellite_V02_base%
  
  
REM
REM ===============================================================
REM ---------- Base Materiaux  -----------
REM ===============================================================
REM
SET Satellite_V02_bddm=%Satellite_V02_base%/BULK/MATERIAUX
  ECHO %INDENT_ECHO_05%      Satellite_V02_bddm = %Satellite_V02_bddm%
  
  
REM
REM =================================
REM ----------- FEM INCLUUDE -----------
REM =================================
REM
SET Satellite_V02_INCLUDE=%Satellite_V02_base%/INCLUDE
  ECHO %INDENT_ECHO_05%      Satellite_V02_INCLUDE = %Satellite_V02_INCLUDE%

  
REM
REM =================================
REM ----------- FEM BULK -----------
REM =================================
REM  
REM
SET Satellite_V02_BULK=%Satellite_V02_base%/BULK
  ECHO %INDENT_ECHO_05%      Satellite_V02_BULK = %Satellite_V02_BULK%
  



ECHO.
ECHO.
ECHO %INDENT_ECHO_02%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ECHO.
  
REM PAUSE




