TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
release {
    DEFINES += ARMA_NO_DEBUG
}

SOURCES += main.cpp

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo
