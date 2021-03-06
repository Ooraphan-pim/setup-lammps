# LAMMPS multiple-machine Makefile

SHELL = /bin/sh
#.IGNORE:

# Definitions

ROOT =	lmp
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Targets

help:
	@echo 'Type "make target" where target is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	for file in $$files; do head -1 $$file; done

clean:
	rm -r Obj_*

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE)  "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) $(INC) Makefile*; fi

# Update Makefile.lib and Makefile.list

makelib:
	@csh Make.csh Makefile.lib

makelist:
	@csh Make.csh Makefile.list

# Packages

package:
	@echo 'Available packages: class2, dpd, granular, kspace'
#	@echo '                    manybody meam molecule, poems, xtc'
	@echo '                    manybody molecule, poems, xtc'
	@echo '  "make yes-name"   include a package'
	@echo '  "make no-name"    exclude a package'
	@echo '  "make yes-all"    include all packages'
	@echo '  "make no-all"     exclude all packages'

yes-all:
#	make yes-class2 yes-dpd yes-granular yes-kspace \
#	yes-manybody yes-meam yes-molecule yes-poems yes-xtc
	make yes-class2 yes-dpd yes-granular yes-kspace \
	yes-manybody yes-molecule yes-poems yes-xtc

no-all:
	@echo 'Removing files, ignore any rm errors ...'
	@cd CLASS2; csh -f Install.csh 0
	@cd DPD; csh -f Install.csh 0
	@cd GRANULAR; csh -f Install.csh 0
	@cd KSPACE; csh -f Install.csh 0
	@cd MANYBODY; csh -f Install.csh 0
#	@cd MEAM; csh -f Install.csh 0
	@cd MOLECULE; csh -f Install.csh 0
	@cd POEMS; csh -f Install.csh 0
	@cd XTC; csh -f Install.csh 0
	@make clean

yes-class2:
	@cd CLASS2; csh -f Install.csh 1
no-class2:
	@echo 'Removing files, ignore any rm errors ...'
	@cd CLASS2; csh -f Install.csh 0
	@make clean

yes-dpd:
	@cd DPD; csh -f Install.csh 1
no-dpd:
	@echo 'Removing files, ignore any rm errors ...'
	@cd DPD; csh -f Install.csh 0
	@make clean

yes-granular:
	@cd GRANULAR; csh -f Install.csh 1
no-granular:
	@echo 'Removing files, ignore any rm errors ...'
	@cd GRANULAR; csh -f Install.csh 0
	@make clean

yes-kspace:
	@cd KSPACE; csh -f Install.csh 1
no-kspace:
	@echo 'Removing files, ignore any rm errors ...'
	@cd KSPACE; csh -f Install.csh 0
	@make clean

yes-manybody:
	@cd MANYBODY; csh -f Install.csh 1
no-manybody:
	@echo 'Removing files, ignore any rm errors ...'
	@cd MANYBODY; csh -f Install.csh 0
	@make clean

#yes-meam:
#	@cd MEAM; csh -f Install.csh 1
#no-meam:
#	@echo 'Removing files, ignore any rm errors ...'
#	@cd MEAM; csh -f Install.csh 0
#	@make clean

yes-molecule:
	@cd MOLECULE; csh -f Install.csh 1
no-molecule:
	@echo 'Removing files, ignore any rm errors ...'
	@cd MOLECULE; csh -f Install.csh 0
	@make clean

yes-poems:
	@cd POEMS; csh -f Install.csh 1
no-poems:
	@echo 'Removing files, ignore any rm errors ...'
	@cd POEMS; csh -f Install.csh 0
	@make clean

yes-xtc:
	@cd XTC; csh -f Install.csh 1
no-xtc:
	@echo 'Removing files, ignore any rm errors ...'
	@cd XTC; csh -f Install.csh 0
	@make clean

# update src files with package files

package-update:
	@csh -f Package.csh CLASS2 update
	@csh -f Package.csh DPD update
	@csh -f Package.csh GRANULAR update
	@csh -f Package.csh KSPACE update
	@csh -f Package.csh MANYBODY update
#	@csh -f Package.csh MEAM update
	@csh -f Package.csh MOLECULE update
	@csh -f Package.csh POEMS update
	@csh -f Package.csh XTC update

# overwrite package files with src files

package-overwrite:
	@csh -f Package.csh CLASS2 overwrite
	@csh -f Package.csh DPD overwrite
	@csh -f Package.csh GRANULAR overwrite
	@csh -f Package.csh KSPACE overwrite
	@csh -f Package.csh MANYBODY overwrite
#	@csh -f Package.csh MEAM overwrite
	@csh -f Package.csh MOLECULE overwrite
	@csh -f Package.csh POEMS overwrite
	@csh -f Package.csh XTC overwrite

# check differences between src and pacakge files

package-check:
	@csh -f Package.csh CLASS2 check
	@csh -f Package.csh DPD check
	@csh -f Package.csh GRANULAR check
	@csh -f Package.csh KSPACE check
	@csh -f Package.csh MANYBODY check
#	@csh -f Package.csh MEAM check
	@csh -f Package.csh MOLECULE check
	@csh -f Package.csh POEMS check
	@csh -f Package.csh XTC check
