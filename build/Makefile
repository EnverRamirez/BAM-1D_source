#
#  $Author: pkubota $
#  $Date: 2007/10/10 20:28:03 $
#  $Revision: 1.4 $
#
#  alteration of the name of the executable and creation of variables and LV EXEC
#  Enter the number of Levels Model (affects SpecDump.f90)
#
# $Id: Makefile,v 1.4 2007/10/10 20:28:03 pkubota Exp $
#
SHELL=/bin/sh
#
# Makefile for Pos-Proccessing output global model
#
# Targets are the suffixes of the system-specific makefiles in
# the makefiles/ directory.
# For example, to build PosGrib for Solaris, give the command
#
#   make solaris
#
# This builds an intermediate library in the util/ directory,
# then builds the nedit and nc executables in the source/ directory.
#
MAKE=make
#LV=28
LV=64
PATH2=../exec
SOURCES_DIR = ../source

all: 
	@echo "Please specify target:"
	@echo "(For example, type \"make linux\" for a Linux system.)"
	@(cd  makefiles && ls -C Makefile* | sed -e 's/Makefile.//g')

.DEFAULT: copy_sources_to_build
	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f ./Makefile.$@ ];\
	#cat Assimilation/SpecDump.f90.in | sed s/"NLEVS=28"/"NLEVS=$(LV)"/g > Assimilation/SpecDump.f90; \
	   then ln -fs makefiles/Makefile.$@ Makefile.$@; fi)
	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Physics/Radiation/RRTMG/ext/rrtmg_sw/Makefile.$@ ];\
	   then ln -fs ../../../../../makefiles/Makefile.$@ Physics/Radiation/RRTMG/ext/rrtmg_sw/Makefile.$@; fi)
	@- (cd .; if [ -f makefiles/Makefile.$@ -a ! -f Physics/Radiation/RRTMG/ext/rrtmg_lw/Makefile.$@ ];\
	   then ln -fs ../../../../../makefiles/Makefile.$@ Physics/Radiation/RRTMG/ext/rrtmg_lw/Makefile.$@; fi)

	#(/opt/modules/3.2.6.6/bin/modulecmd bash load craype-hugepages64M)
	(cd .;cd Physics/Radiation/RRTMG/ext/rrtmg_sw;$(MAKE) -f Makefile.$@)
	(cd .;cp -fp Physics/Radiation/RRTMG/ext/rrtmg_sw/*.o .)
	(cd .;echo ll.mod > Physics/Radiation/RRTMG/ext/rrtmg_sw/ll.mod ;cp -fp Physics/Radiation/RRTMG/ext/rrtmg_sw/*.mod .)
	(cd .;cd Physics/Radiation/RRTMG/ext/rrtmg_lw;$(MAKE) -f Makefile.$@)
	(cd .;cp -fp Physics/Radiation/RRTMG/ext/rrtmg_lw/*.o .)
	(cd .;echo ll.mod > Physics/Radiation/RRTMG/ext/rrtmg_lw/ll.mod ;cp -fp Physics/Radiation/RRTMG/ext/rrtmg_lw/*.mod .)
	(cd .;$(MAKE) -f Makefile.$@)
	(mv ParModel_MPI $(PATH2)/ParModel_MPI)
	(rm -f Makefile.$@)

# We need a "dev-all" target that builds the docs plus binaries, but
# that doesn't work since we require the user to specify the target.  More
# thought is needed


copy_sources_to_build:
	@echo " =====> copying sources to build ..."
	-rsync -Cravzp $(SOURCES_DIR)/* .


clean:
	@echo "$(MAKE)"

	(cd .;$(MAKE) -f Makefile.common clean)
	(cd Physics/Radiation/RRTMG/ext/rrtmg_sw;$(MAKE) -f Makefile.common clean)
	(cd Physics/Radiation/RRTMG/ext/rrtmg_lw;$(MAKE) -f Makefile.common clean)
