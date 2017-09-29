SUBDIRS = functions 

jagsmod_LTLIBRARIES = nuclear.la

nuclear_la_SOURCES = nuclear.cc

nuclear_la_CPPFLAGS = -I$(includedir)/JAGS

nuclear_la_LIBADD = functions/nuclear.la
if WINDOWS
nuclear_la_LIBADD += -ljags-$(JAGS_MAJOR) -ljrmath-0
else
nuclear_la_LIBADD += -ljags -ljrmath
endif

nuclear_la_LDFLAGS = -module -avoid-version
if WINDOWS
nuclear_la_LDFLAGS += -no-undefined
endif
