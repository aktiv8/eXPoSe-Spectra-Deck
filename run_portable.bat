@echo off
REM Run from source using a Python folder carried with the program.
REM Copy a full python.org install (it includes tkinter) into a folder called
REM "python" beside this file; the first run builds .venv here, later runs are
REM fast. Without that folder this behaves like run.bat.

cd /d "%~dp0"

if exist "%~dp0python\python.exe" (
    "%~dp0python\python.exe" launch.py %*
    goto end
)

call "%~dp0run.bat" %*
goto :eof

:end
if %errorlevel% neq 0 (
    echo.
    echo The launcher exited with an error. Review the messages above.
    pause
)
