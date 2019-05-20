noinst_LTLIBRARIES = nuclear.la

nuclear_la_CPPFLAGS = -I$(top_srcdir)/src \
												 -I$(includedir)/JAGS

nuclear_la_LDFLAGS = -no-undefined -module -avoid-version

nuclear_la_SOURCES = sfactorTdn.cc sfactor3Hedp.cc   \
    sfactor3Hedpx.cc \
    sfactorTdn2.cc sfactor3Hdnx.cc sigma7Benpx.cc sfactor7Lipnx.cc

noinst_HEADERS = sfactorTdn.h sfactor3Hedp.h      \
    sfactor3Hedpx.h \
    sfactorTdn2.h sfactor3Hdnx.h sigma7Benpx.h sfactor7Lipnx.h
