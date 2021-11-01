
SET CHIMERAX_EXE="c:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe"
FOR %%A in (%*) DO (

IF "%%A" == "release" set CHIMERAX_EXE="c:\Program Files\ChimeraX\bin\ChimeraX-console.exe"
)

for %%A in (%*) DO (
IF "%%A" == "clean" (
	%CHIMERAX_EXE% --nogui --safemode --exit --cmd "devel clean ."
	BREAK
)
)

for %%A in (%*) DO (
IF "%%A" == "app-install" (
	%CHIMERAX_EXE% --nogui --safemode --exit --cmd "devel install ."
	BREAK
)
)
