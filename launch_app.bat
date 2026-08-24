@echo off
setlocal enabledelayedexpansion

:: 1. Check if Rscript is already in system PATH
where Rscript >nul 2>&1
if %errorlevel% equ 0 (
    set "RSCRIPT_PATH=Rscript"
    goto RUN
)

:: 2. Search standard Program Files (64-bit)
set "R_DIR=C:\Program Files\R"
if exist "%R_DIR%" (
    for /f "delims=" %%I in ('dir /b /o:-n "%R_DIR%\R-*" 2^>nul') do (
        if exist "%R_DIR%\%%I\bin\Rscript.exe" (
            set "RSCRIPT_PATH=%R_DIR%\%%I\bin\Rscript.exe"
            goto RUN
        )
        if exist "%R_DIR%\%%I\bin\x64\Rscript.exe" (
            set "RSCRIPT_PATH=%R_DIR%\%%I\bin\x64\Rscript.exe"
            goto RUN
        )
    )
)

:: 3. Search Program Files (x86) (32-bit or legacy installs)
set "R_DIR_X86=C:\Program Files (x86)\R"
if exist "%R_DIR_X86%" (
    for /f "delims=" %%I in ('dir /b /o:-n "%R_DIR_X86%\R-*" 2^>nul') do (
        if exist "%R_DIR_X86%\%%I\bin\Rscript.exe" (
            set "RSCRIPT_PATH=%R_DIR_X86%\%%I\bin\Rscript.exe"
            goto RUN
        )
        if exist "%R_DIR_X86%\%%I\bin\i386\Rscript.exe" (
            set "RSCRIPT_PATH=%R_DIR_X86%\%%I\bin\i386\Rscript.exe"
            goto RUN
        )
    )
)

:: 4. Fallback: Query Windows Registry for R install location
for /f "tokens=2*" %%A in ('reg query "HKLM\SOFTWARE\R-core\R" /v "InstallPath" 2^>nul') do (
    if exist "%%B\bin\Rscript.exe" (
        set "RSCRIPT_PATH=%%B\bin\Rscript.exe"
        goto RUN
    )
)

:: 5. Error state
echo [ERROR] R installation could not be found on this computer.
echo Checked: System PATH, Program Files, Program Files (x86), and Registry.
pause
exit /b 1

:RUN
echo Found Rscript at: "%RSCRIPT_PATH%"
echo Launching Shiny App...

"%RSCRIPT_PATH%" launch.R

pause