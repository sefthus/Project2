TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    the_jacobi_method.cpp

HEADERS += \
    catch.hpp \
    the_jacobi_method.hpp

LIBS += -LC:\Armadillo\examples\lib_win64 -lblas_win64_MT
LIBS += -LC:\Armadillo\examples\lib_win64 -llapack_win64_MT

INCLUDEPATH += C:\Armadillo\include
DEPENDPATH += C:\Armadillo\include
