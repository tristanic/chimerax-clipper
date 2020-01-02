
SET CHIMERAX_EXE="c:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe"
FOR %%A in (%*) DO (

IF "%%A" == "release" set CHIMERAX_EXE="c:\Program Files\ChimeraX\bin\ChimeraX-console.exe"
)

for %%A in (%*) DO (
IF "%%A" == "clean" (
	%CHIMERAX_EXE% --nogui --cmd "devel clean .; exit"
	BREAK
)
)

for %%A in (%*) DO (
IF "%%A" == "app-install" (
	%CHIMERAX_EXE% --nogui --cmd "devel install .; exit"
	BREAK
)
)
