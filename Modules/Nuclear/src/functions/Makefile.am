noinst_LTLIBRARIES = nuclear.la

nuclear_la_CPPFLAGS = -I$(top_srcdir)/src \
												 -I$(includedir)/JAGS

nuclear_la_LDFLAGS = -no-undefined -module -avoid-version

nuclear_la_SOURCES = sfactorTdn_fast.cc sfactorTdn.cc

noinst_HEADERS = sfactorTdn_fast.h  sfactorTdn.h 
