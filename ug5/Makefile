#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*             This file is part of the program and software framework       *
#*                  UG --- Ubquity Generator Framework                       *
#*                                                                           *
#*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  UG is distributed under the terms of the ZIB Academic Licence.           *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with UG; see the file COPYING. If not email to scip@zib.de.        *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#@file    Makefile
#@brief   Makefile for a parallel steiner tree problem 
#@author  Yuji Shinano

UG_VERSION 	= 0.8.5
#-----------------------------------------------------------------------------
# COMM 
#-----------------------------------------------------------------------------
COMM		= pth
#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------
UGDIR		=	../../../scipoptsuite-5.0.0/ug
SCIPDIR         =       $(UGDIR)/lib/scip
STPDIR		= 	../scip5
#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project
#-----------------------------------------------------------------------------
# set the UGLIB 
#-----------------------------------------------------------------------------
UGLIB		=	ugscip-$(UG_VERSION).$(BASE).$(COMM)
#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------
TEST            =      short
TIME            =      3600
SETTINGS        =      default
VERBOSE         =      false
OPT             =      opt
COMP            =      gnu
LPS             =      spx2
DET             =      false
ZLIB            =      true
THREADS         =      4
#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

ifeq ($(COMM), pth)
MAINNAME	=	fscip
endif
ifeq ($(COMM), mpi)
MAINNAME	=	parascip
endif
UGOBJ		=	logreg_plugins.o
#MAINOBJ	=	reader_logreg.o probdata_logreg.o mysolpool.o \
#				convenient_tool.o call_cblas.o get_lineardependence.o \
#				relax_newton.o logistic_regression.o stepsize.o

MAINOBJXX	=	reader_logreg.o probdata_logreg.o mysolpool.o \
					convenient_tool.o call_cblas.o get_lineardependence.o \
					relax_newton.o logistic_regression.o stepsize.o \
					heur_backward.o

USRLDFLAGS= -llapack -lblas
#USRLDFLAGS=  -L/usr/lib64/libopenblas.so -lopenblas -llapack -lpthread

UGSRC		=	$(addprefix $(SRCDIR)/,$(UGOBJ:.o=.cpp))
MAINSRC		=	$(addprefix $(STPDIR)/src/,$(MAINOBJ:.o=.c))
MAINSRC         +=      $(addprefix $(STPDIR)/src/,$(MAINOBJXX:.o=.cpp))
MAINSRC         +=      $(addprefix $(SRCDIR)/,$(MAINCXXOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
UGOBJFILES	=	$(addprefix $(OBJDIR)/,$(UGOBJ))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))
MAINOBJFILES		+=	$(addprefix $(OBJDIR)/,$(MAINOBJXX))
MAINCXXOBJFILES =	$(addprefix $(OBJDIR)/,$(MAINCXXOBJ))

ifeq ($(CPOPT),true)
FLAGS           +=      -DWITH_CPOPTIMIZER -DIL_STD -I$(LIBDIR)/cpinc -I$(LIBDIR)/cpxinc
FLAGS           +=      -I$(LIBDIR)/ilconcertinc
LDFLAGS         +=       $(LIBDIR)/libcp.a $(LIBDIR)/libconcert.a  $(LIBDIR)/libcplex.a -lpthread
endif

# LASTSETTINGS    =       $(OBJDIR)/make.lastsettings

# maybe setting the include path to SCIP's application directory is enough
CFLAGS		+=	-I$(SCIPDIR)/src -DWITH_UG -I$(STPDIR)/src 
ifeq ($(COMM), pth)
CXX             =       g++
LINKCXX         =       g++
CXXFLAGS	+=	-I$(UGDIR)/src -D_COMM_PTH -I$(STPDIR)/src
endif
ifeq ($(COMM), mpi)
CC              =       mpicc
CXX             =       mpicxx
LINKCXX         =       mpicxx
CXXFLAGS	+=	-I$(UGDIR)/src -D_COMM_MPI_WORLD -DMPICH_IGNORE_CXX_SEEK -I$(STPDIR)/src
endif

LDFLAGS		+=	-lpthread

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(UGOBJFILES) $(MAINOBJFILES) $(MAINCXXOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) $(SCIPDIR)/lint/scip.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

.PHONY: doc
doc:
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scip.css);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/scippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/miniscippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scipfooter.html footer.html);
		cd doc; $(DOXY) scheduler.dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(CHECKDIR):
		@-mkdir -p $(CHECKDIR)

.PHONY: test
test:		$(STPDIR)
		@-(ln -fs $(STPDIR)/data);
		@-(cd check && ln -fs ../$(STPDIR)/check/testset);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
ifeq ($(COMM),pth)
		cd check; \
		$(SHELL) ./check_fscip.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(MEM) $(CONTINUE) $(LOCK) $(LPS) $(THREADS);
endif
ifeq ($(COMM),mpi)
		cd check; \
		$(SHELL) ./check_parascip.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(MEM) $(CONTINUE) $(LOCK) $(LPS) $(THREADS);
endif

# include local targets
-include $(SCIPDIR)/make/local/make.targets


.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'

-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:  $(CPOPTDEP)
ifneq ($(CPOPT),$(LAST_CPOPT))
		@-touch $(CPOPTSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_CPOPT=$(CPOPT)" >> $(LASTSETTINGS)

-include	$(MAINDEP)

# $(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES) $(MAINCXXOBJFILES)
$(MAINFILE):	$(BINDIR) $(OBJDIR) $(UGOBJFILES) $(MAINOBJFILES) $(MAINCXXOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(UGOBJFILES) $(MAINOBJFILES) $(MAINCXXOBJFILES) $(LINKCXXSCIPALL) \
		$(LINKCXX_L)$(UGDIR)/lib $(LINKCXX_l)$ $(UGLIB) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(STPDIR)/src/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(STPDIR)/src/%.cpp
		@echo "-> compiling $@"
		#$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c -g $< $(CXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
