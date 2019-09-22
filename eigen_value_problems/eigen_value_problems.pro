TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/Cellar
LIBS += -larmadillo -llapack -lblas



INCLUDEPATH += D:\armadillo-9.700.2\include
DEPENDPATH += D:\armadillo-9.700.2\include
LIBS += \
    -LD:\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT
