TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

release {
    DEFINES += ARMA_NO_DEBUG
}


INCLUDEPATH += /usr/local/include
#LIBS += -L/usr/local/lib -framework Accelerate
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo
QMAKE_MAC_SDK = macosx10.11 # Something wrong with SDK... (idk lol)

SOURCES += main.cpp \
    #arma_solve.cpp \
    #num_solve.cpp \
    #lu_decomposition.cpp \
    numerical.cpp \
    lu_decomp.cpp

#RESOURCES += qml.qrc

# Default rules for deployment.
#include(deployment.pri)

#DISTFILES += \
#    plotting_results.py

HEADERS += \
    #arma_solve.h \
    #num_solve.h \
    #lu_decomposition.h \
    numerical.h \
    lu_decomp.h

